
import argparse
import logging
import hashlib
import sqlite3
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate UPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--nrp', action='store', help='Path to NCBI NRP fasta file.')
parser.add_argument('--pcla-proteins', action='store', dest='pcla_proteins', help='Path to NCBI NRP PCLA proteins file.')
parser.add_argument('--pcla-clusters', action='store', dest='pcla_clusters', help='Path to NCBI NRP PCLA clusters file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
nrp_path = Path(args.nrp).resolve()
pcla_proteins_path = Path(args.pcla_proteins).resolve()
pcla_clusters_path = Path(args.pcla_clusters).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-NRP - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger_ups = logging.getLogger('UPS')
logger_psc = logging.getLogger('PSC')


print('import PCLA proteins and cluster information...')
pcla_cluster_genes = {}
with pcla_clusters_path.open() as fh:
    for line in fh:
        (id_, id, product, proteins, organisms, conserved_in_organism, conserved_in_taxid, gene, hmm) = line.split('\t')
        gene = gene.strip()
        if(gene != ''):
            pcla_cluster_genes[id] = gene

nrp_genes = {}
with pcla_proteins_path.open() as fh:
    for line in fh:
        (cluster_id, id, definition, organism, taxid, length) = line.split('\t')
        gene = pcla_cluster_genes.get(cluster_id, None)
        if(gene is not None):
            nrp_genes[id] = gene
del pcla_cluster_genes
print("\tfound %d NRP / gene annotations" % len(nrp_genes))


print('lookup UPS by NRP hash & update WP_* / gene annotations...')
nrps_processed = 0
nrps_not_found = 0
nrps_lenth_mm = 0
nrps_updated = 0
nrpcs_updated = 0
ncbi_nrp_path = Path(args.nrp).resolve()
with ncbi_nrp_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()
    conn.row_factory = sqlite3.Row
    c_ups = conn.cursor()
    for record in SeqIO.parse(fh, 'fasta'):
        nrps_processed += 1
        nrp_id = record.id
        seq = str(record.seq).upper()
        hash = hashlib.md5(seq.encode()).hexdigest()
        c_ups.execute('SELECT length, uniref100_id, uniref90_id FROM ups WHERE hash=?', (hash,))
        r = c_ups.fetchone()
        if(r is None):
            print("UPS hash not found! nrp-id=%s, hash=%s" % (nrp_id, hash))
            nrps_not_found += 1
            continue
        elif(r['length'] != len(seq)):
            print("UPS length mismatch! nrp-id=%s, hash=%s, ups-length=%d, nrp-length=%d" % (nrp_id, hash, r['length'], len(seq)))
            nrps_lenth_mm += 1
            continue
        else:
            gene = nrp_genes.get(nrp_id, None)
            logger_ups.debug(
                'ncbi-nrp-id=%s, hash=%s, uniref100=%s, uniref90=%s, gene=%s', 
                nrp_id, hash, r['uniref100_id'], r['uniref90_id'], gene
            )
            c_ups.execute('UPDATE ups SET ncbi_nrp_id=? WHERE hash=?', (nrp_id[3:], hash))  # annotate UPS with NCBI nrp id (WP_*)
            logger_ups.info('UPDATE ups SET ncbi_nrp_id=%s WHERE hash=%s', nrp_id[3:], hash)
            nrps_updated += 1
            if(gene is not None):  # annotate UniRef90 with NCBI PCLA cluster's gene if present
                conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, r['uniref90_id']))
                logger_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', gene, r['uniref90_id'])
                nrpcs_updated += 1
        if((nrps_processed % 1000000) == 0):
            conn.commit()
            print("\t... %d" % nrps_processed)
    c_ups.close()
    conn.commit()

print('\n')
print("NRPs processed: %d" % nrps_processed)
logger_ups.debug('summary: # UPS with annotated NRP IDs=%d', nrps_updated)
print("NRPs with annotated WP_* id: %d" % nrps_updated)
logger_psc.debug('summary: # PSC with annotated genes=%d', nrpcs_updated)
print("NRPCs with annotated gene names: %d" % nrpcs_updated)
print("NRPs not found: %d" % nrps_not_found)
print("NRPs with length mismatches: %d" % nrps_lenth_mm)

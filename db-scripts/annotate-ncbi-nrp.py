
import argparse
import logging
import hashlib
import sqlite3
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
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
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


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
print("\tfound %i NRP / gene annotations" % len(nrp_genes))


print('lookup UPS by NRP hash & update WP_* / gene annotations...')
nrps_processed = 0
nrps_not_found = 0
nrps_wo_ips = 0
ups_updated = 0
ips_updated = 0
psc_updated = 0
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
    for record in SeqIO.parse(fh, 'fasta'):
        nrps_processed += 1
        nrp_id = record.id
        seq = str(record.seq).upper()
        seq_hash = hashlib.md5(seq.encode())
        seq_hash_digest = seq_hash.digest()
        seq_hash_hexdigest = seq_hash.hexdigest()
        rec_ups = conn.execute('SELECT length, uniref100_id FROM ups WHERE hash=?', (seq_hash_digest,)).fetchone()
        if(rec_ups is not None):
            assert rec_ups['length'] == len(seq), "Detected hash duplicate with different seq length! hash=%s, NCRBI-NRP-id=%s, UniParc-id=%s, NCRBI-NRP-length=%s, db-length=%s" % (seq_hash_hexdigest, nrp_id, rec_ups['uniparc_id'], length, rec_ups['length'])
            if(rec_ups['uniref100_id']):
                rec_ips = conn.execute('SELECT * FROM ips WHERE uniref100_id=?', (rec_ups['uniref100_id'],)).fetchone()
                if(rec_ips is not None):
                    gene = nrp_genes.get(nrp_id, None)
                    log_ups.debug(
                        'ncbi-nrp-id=%s, hash=%s, uniref100=%s, uniref90=%s, gene=%s', 
                        nrp_id, seq_hash.hexdigest(), rec_ips['uniref100_id'], rec_ips['uniref90_id'], gene
                    )
                    conn.execute('UPDATE ups SET ncbi_nrp_id=? WHERE hash=?', (nrp_id[3:], seq_hash_digest))  # annotate IPS with NCBI nrp id (WP_*)
                    log_ups.info('UPDATE ups SET ncbi_nrp_id=%s WHERE hash=%s', nrp_id[3:], seq_hash_hexdigest)
                    ups_updated += 1
                    if(gene is not None):  # annotate IPS or PSC with NCBI PCLA cluster's gene if present
                        if(rec_ips['uniref90_id'] is None):
                            conn.execute('UPDATE ips SET gene=? WHERE uniref100_id=?', (gene, rec_ups['uniref100_id']))
                            log_ips.info('UPDATE ips SET gene=%s WHERE uniref100_id=%s', gene, rec_ups['uniref100_id'])
                            ips_updated += 1
                        else:
                            conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, rec_ips['uniref90_id']))
                            log_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', gene, rec_ips['uniref90_id'])
                            psc_updated += 1
                else:
                    nrps_wo_ips += 1
            else:
                nrps_wo_ips += 1
        else:
            nrps_not_found += 1
        if((nrps_processed % 1000000) == 0):
            conn.commit()
            print("\t... %i" % nrps_processed)
    conn.commit()

print('\n')
print("NRPs processed: %i" % nrps_processed)

log_ups.debug('summary: # UPS with annotated NRP IDs=%i', ups_updated)
print("UPSs with annotated WP_* id: %i" % ups_updated)

log_ips.debug('summary: # IPSs with annotated genes=%i', ips_updated)
print("IPSs with annotated gene names: %i" % ips_updated)

log_psc.debug('summary: # PSC with annotated genes=%i', psc_updated)
print("PSCs with annotated gene names: %i" % psc_updated)

print("NRPs not found: %i" % nrps_not_found)
print("NRPs w/o IPS: %i" % nrps_wo_ips)


import argparse
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
    for record in SeqIO.parse(fh, 'fasta'):
        nrps_processed += 1
        nrp_id = record.id
        seq = str(record.seq).upper()
        hash = hashlib.md5(seq.encode()).hexdigest()
        try:
            conn.row_factory = sqlite3.Row
            c_ups = conn.cursor()
            c_ups.execute('SELECT length, uniref90_id FROM ups WHERE hash=?', (hash,))
            r = c_ups.fetchone()
            if(r is None):
                # print("WARNING: could not find UPS! nrp-id=%s, hash=%s" % (nrp_id, hash))
                nrps_not_found += 1
                continue
            elif(r['length'] != len(seq)):
                # print("WARNING: length mismatch! nrp-id=%s, hash=%s, ups-length=%d, nrp-length=%d" % (nrp_id, hash, r['length'], len(seq)))
                nrps_lenth_mm += 1
                continue
            else:
                c_ups.execute('UPDATE ups SET ncbi_nrp_id=? WHERE hash=?', (nrp_id[3:], hash))  # annotate UPS with NCBI nrp id (WP_*)
                c_ups.close()
                nrps_updated += 1
                gene = nrp_genes.get(nrp_id, None)
                if(gene is not None):  # annotate UniRef90 with NCBI PCLA cluster's gene if present
                    c_psc = conn.cursor()
                    c_psc.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, r['uniref90_id']))
                    c_psc.close()
                    nrpcs_updated += 1
        except sqlite3.Error as e:
            print("SQL ERROR: could not set annotation for nrp-id=%s" % nrp_id)
            print(e)
        if((nrps_processed % 1000000) == 0):
            conn.commit()
            print("\t... %d" % nrps_processed)
    conn.commit()

print('\n')
print("NRPs processed: %d" % nrps_processed)
print("NRPs with annotated WP_* id: %d" % nrps_updated)
print("NRPCs with annotated gene names: %d" % nrpcs_updated)
print("NRPs not found: %d" % nrps_not_found)
print("NRPs with length mismatches: %d" % nrps_lenth_mm)

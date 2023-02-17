
import argparse
import logging
import hashlib
import sqlite3

from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO


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
pcla_cluster_annotations = {}
with pcla_clusters_path.open() as fh, alive_bar() as bar:
    for line in fh:
        (id_, id, product, proteins, organisms, conserved_in_organism, conserved_in_taxid, gene, hmm) = line.split('\t')
        gene = gene.strip()
        product = product.strip().replace('"', '')
        if(gene == ''):
            gene = None
        if(product == '' or product.lower() == 'hypothetical protein' or product.lower() == 'uncharacterized protein' or product.lower() == 'uncharacterized conserved protein'):
            product = None
        pcla_cluster_annotations[id] = (gene, product)
    bar()

nrp_annotations = {}
with pcla_proteins_path.open() as fh, alive_bar() as bar:
    for line in fh:
        (cluster_id, id, definition, organism, taxid, length) = line.split('\t')
        if(cluster_id in pcla_cluster_annotations):
            nrp_annotations[id] = pcla_cluster_annotations[cluster_id]
        bar()
del pcla_cluster_annotations
print(f'found {len(nrp_annotations)} NRP / gene annotations')


print('lookup UPS by NRP hash & update WP_* / gene annotations...')
nrps_processed = 0
nrps_not_found = 0
nrps_wo_ips = 0
ups_updated = 0
psc_updated_gene = 0
psc_updated_product = 0
ncbi_nrp_path = Path(args.nrp).resolve()
with ncbi_nrp_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, alive_bar() as bar:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
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
            assert rec_ups['length'] == len(seq), f"Detected hash duplicate with different seq length! hash={seq_hash_hexdigest}, NCRBI-NRP-id={nrp_id}, UniParc-id={rec_ups['uniparc_id']}, NCRBI-NRP-length={length}, db-length={rec_ups['length']}"
            if(rec_ups['uniref100_id']):
                rec_ips = conn.execute('SELECT * FROM ips WHERE uniref100_id=?', (rec_ups['uniref100_id'],)).fetchone()
                if(rec_ips is not None):
                    (gene, product) = nrp_annotations.get(nrp_id, (None, None))
                    log_ups.debug(
                        'ncbi-nrp-id=%s, hash=%s, uniref100=%s, uniref90=%s, gene=%s, product=%s',
                        nrp_id, seq_hash.hexdigest(), rec_ips['uniref100_id'], rec_ips['uniref90_id'], gene, product
                    )
                    conn.execute('UPDATE ups SET ncbi_nrp_id=? WHERE hash=?', (nrp_id[3:], seq_hash_digest))  # annotate IPS with NCBI nrp id (WP_*)
                    log_ups.info('UPDATE ups SET ncbi_nrp_id=%s WHERE hash=%s', nrp_id[3:], seq_hash_hexdigest)
                    ups_updated += 1
                    if(rec_ips['uniref90_id'] is not None):
                        if(gene is not None):
                            conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (gene, rec_ips['uniref90_id']))
                            log_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', gene, rec_ips['uniref90_id'])
                            psc_updated_gene += 1
                        if(product is not None):
                            rec_psc = conn.execute('SELECT * FROM psc WHERE uniref90_id=?', (rec_ips['uniref90_id'],)).fetchone()
                            if(rec_psc is not None and rec_psc['product'] is None):
                                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, rec_ips['uniref90_id']))
                                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, rec_ips['uniref90_id'])
                                psc_updated_product += 1
                else:
                    nrps_wo_ips += 1
            else:
                nrps_wo_ips += 1
        else:
            nrps_not_found += 1
        if((nrps_processed % 1000000) == 0):
            conn.commit()
        bar()
    conn.commit()
print(f'NRPs processed: {nrps_processed}')
log_ups.debug('summary: # UPS with annotated NRP IDs=%i', ups_updated)
print(f'UPSs with annotated WP_* id: {ups_updated}')
log_psc.debug('summary: # PSC with annotated genes=%i', psc_updated_gene)
print(f'PSCs with annotated gene names: {psc_updated_gene}')
log_psc.debug('summary: # PSC with annotated products=%i', psc_updated_product)
print(f'PSCs with annotated gene products: {psc_updated_product}')
print(f'NRPs not found: {nrps_not_found}')
print(f'NRPs w/o IPS: {nrps_wo_ips}')

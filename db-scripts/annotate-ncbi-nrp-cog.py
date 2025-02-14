
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
parser.add_argument('--cog-ids', action='store', dest='cog_ids', help='Path to NCBI COG id / functional class file.')
parser.add_argument('--nrp-cog-mapping', action='store', dest='nrp_cog_mapping', help='Path to GI / COG mapping file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
nrp_path = Path(args.nrp).resolve()
cog_ids_path = Path(args.cog_ids).resolve()
nrp_cog_mapping_path = Path(args.nrp_cog_mapping).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-NRP-COG - %(levelname)s - %(message)s',
    level=logging.INFO
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


print('import NCBI COG id / function class information...')
cogs_by_id = {}
with cog_ids_path.open() as fh:
# with cog_ids_path.open(encoding='windows-1252') as fh:
    for line in fh:
        if(line[0] != '#'):
            (id, cat, product, gene, pathways, pubmed, pdb) = line.split('\t')
            if(product.lower() == 'hypothetical protein' or product.lower() == 'uncharacterized protein' or product.lower() == 'uncharacterized conserved protein' or product.lower() == 'Uncharacterized membrane protein'):
                product = None
            if(gene == ''):
                gene = None
            if(cat == ''):
                cat = None
            cogs_by_id[id] = {
                'id': id,
                'cat': cat,
                'product': product,
                'gene': gene
            }
print(f'\tstored COGs: {len(cogs_by_id)}')

print('import NCBI GI / COG mapping information...')
nrp_to_cog = {}
with nrp_cog_mapping_path.open() as fh:
    for line in fh:
        cols = line.strip().split(',')
        # (gene_id, assembly_id, gb_id, protein_length, footprint_coordinates_protein, cog_footprint_length, cog_id, reserved_, cog_membership_class, bitscore, evalue, profile_length, footprint_coordinates_profile)
        cog_id = cols[6]
        nrp_id = cols[2]
        if cog_id != ''  and  nrp_id.startswith('WP_') and cog_id in cogs_by_id:
            nrp_to_cog[nrp_id] = cog_id
print(f'\tmapped COG IDs: {len(nrp_to_cog)}')


print('lookup UPS by NRP hash & update WP_* / gene annotations...')
nrps_processed = 0
nrps_not_found = 0
nrps_wo_ips = 0
ups_updated = 0
psc_updated_gene = 0
psc_updated_product = 0
ncbi_nrp_path = Path(args.nrp).resolve()
uniref90_with_cog_annotations = set()
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
            assert rec_ups['length'] == len(seq), f"Detected hash duplicate with different seq length! hash={seq_hash_hexdigest}, NCRBI-NRP-id={nrp_id}, UniParc-id={rec_ups['uniparc_id']}, NCRBI-NRP-length={rec_ups['length']}, db-length={rec_ups['length']}"
            conn.execute('UPDATE ups SET ncbi_nrp_id=? WHERE hash=?', (nrp_id[3:], seq_hash_digest))  # annotate UPS with NCBI nrp id (WP_*)
            log_ups.info('UPDATE ups SET ncbi_nrp_id=%s WHERE hash=%s', nrp_id[3:], seq_hash_hexdigest)
            ups_updated += 1
            if(rec_ups['uniref100_id']):
                rec_ips = conn.execute('SELECT * FROM ips WHERE uniref100_id=?', (rec_ups['uniref100_id'],)).fetchone()
                if(rec_ips is not None):
                    uniref90_id = rec_ips['uniref90_id']
                    cog_id = nrp_to_cog.get(nrp_id, None)
                    if(uniref90_id is not None  and  cog_id is not None  and  uniref90_id not in uniref90_with_cog_annotations):  # skip known PSC/COG annotations
                        uniref90_with_cog_annotations.add(uniref90_id)
                        cog = cogs_by_id[cog_id]  # assured in line 67
                        conn.execute('UPDATE psc SET cog_id=?, cog_category=? WHERE uniref90_id=?', (cog['id'][3:], cog['cat'], uniref90_id))
                        conn.execute('UPDATE psc SET cog_id=?, cog_category=? WHERE uniref90_id=?', (cog['id'][3:], cog['cat'], uniref90_id))
                        log_psc.info('UPDATE psc SET cog_id=%s, cog_category=%s WHERE uniref90_id=%s', cog['id'], cog['cat'], uniref90_id)
                        
                        rec_psc = conn.execute('SELECT gene, product FROM psc WHERE uniref90_id=?', (uniref90_id,)).fetchone()  # get current PSC gene/product annotation
                        if(rec_psc is not None):
                            if(cog['gene'] is not None  and  rec_psc['gene'] is None):  # update PSC gene annotation if non-existent
                                conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (cog['gene'], uniref90_id))
                                log_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', cog['gene'], uniref90_id)
                                psc_updated_gene += 1
                            if(cog['product'] is not None  and  rec_psc['product'] is None):  # update PSC product annotation if non-existent
                                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (cog['product'], uniref90_id))
                                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', cog['product'], uniref90_id)
                                psc_updated_product += 1
                else:
                    nrps_wo_ips += 1
            else:
                nrps_wo_ips += 1
        else:
            nrps_not_found += 1
        if((nrps_processed % 1_000_000) == 0):
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

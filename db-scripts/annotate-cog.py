import argparse
import logging
import sqlite3

from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--alignments', action='store', help='Path to diamond alignment file.')
parser.add_argument('--cog-ids', action='store', dest='cog_ids', help='Path to NCBI COG id / functional class file.')
parser.add_argument('--gi-cog-mapping', action='store', dest='gi_cog_mapping', help='Path to GI / COG mapping file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
alignments_path = Path(args.alignments).resolve()
cog_ids_path = Path(args.cog_ids).resolve()
gi_cog_mapping_path = Path(args.gi_cog_mapping).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - COG - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


print('import NCBI COG id / function class information...')
cog_id_fclass = {}
with cog_ids_path.open(encoding='windows-1252') as fh:
    for line in fh:
        if(line[0] != '#'):
            (id, cat, product, gene, pathways, pubmed, pdb) = line.split('\t')
            if(product.lower() == 'hypothetical protein' or product.lower() == 'uncharacterized protein' or product.lower() == 'uncharacterized conserved protein'):
                product = None
            if(gene == ''):
                gene = None
            cog_id_fclass[id] = {
                'id': id[3:],
                'cat': cat,
                'product': product,
                'gene': gene
            }
print(f'\tstored COG ids: {len(cog_id_fclass)}')

print('import NCBI GI / COG mapping information...')
gb_id_cog = {}
with gi_cog_mapping_path.open() as fh:
    for line in fh:
        cols = line.strip().split(',')
        # (gene_id, assembly_id, gb_id, protein_length, footprint_coordinates_protein, cog_footprint_length, cog_id, reserved_, cog_membership_class, bitscore, evalue, profile_length, footprint_coordinates_profile)
        cog_id = cols[6]
        gb_id = cols[2]
        cog = cog_id_fclass.get(cog_id, None)
        if(cog is not None):
            gb_id_cog[gb_id] = cog
del cog_id_fclass
print(f'\tmapped COG ids: {len(gb_id_cog)}')

print('parse PSC alignments and add NCBI COG annotation...')
psc_processed = 0
psc_updated = 0
with alignments_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    conn.row_factory = sqlite3.Row
    for line in fh:
        (uniref90_id, sseqid, stitle, length, pident, qlen, slen, evalue) = line.strip().split('\t')
        length = int(length)
        qcov = length / int(qlen)
        scov = length / int(slen)
        if(qcov >= 0.8 and scov >= 0.8 and float(pident) >= 90. and float(evalue) < 1e-6):
            cog = gb_id_cog.get(sseqid, None)
            if(cog is not None):
                conn.execute('UPDATE psc SET cog_id=?, cog_category=? WHERE uniref90_id=?', (cog['id'], cog['cat'], uniref90_id))
                log.info('UPDATE psc SET cog_id=%s, cog_category=%s WHERE uniref90_id=%s', cog['id'], cog['cat'], uniref90_id)
                psc_updated += 1
                rec_psc = conn.execute('SELECT * FROM psc WHERE uniref90_id=?', (uniref90_id,)).fetchone()
                if(rec_psc is not None):
                    if(cog['product'] is not None):
                        conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (cog['product'], uniref90_id))
                        log.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', cog['product'], uniref90_id)
                    if(cog['gene'] is not None):
                        conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (cog['gene'], uniref90_id))
                        log.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', cog['gene'], uniref90_id)
        psc_processed += 1
        if((psc_processed % 100000) == 0):
            conn.commit()
            print(f'\t... {psc_processed}')
    conn.commit()

print('\n')
print(f'PSCs processed: {psc_processed}')
print(f'PSCs with annotated COG id/category: {psc_updated}')
log.debug('summary: PSC annotated=%i', psc_updated)

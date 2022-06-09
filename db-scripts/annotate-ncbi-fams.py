import argparse
import logging
import re
import sqlite3

from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate PSCs by NCBIfams.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--hmms', action='store', help='Path to NCBIfams HMM description file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to NCBIfams HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
hmms_path = Path(args.hmms).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-AMR - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_psc = logging.getLogger('PSC')


psc_annotated = 0
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    print('parse NCBIfam AMR HMMs...')
    amr_hmms = {}
    with hmms_path.open() as fh:
        for line in fh:
            (
                hmm_accession, tc_1, tc_2, gene_symbol, length, hmm_description, scope, hmm_type, subtype, clazz, subclass, db_version, hierarchy_node
            ) = line.split('\t')
            hmm = {
                'acc': hmm_accession,
                'gene': gene_symbol,
                'product': hmm_description.strip()
            }
            if(hmm_accession != '-'):
                amr_hmms[hmm_accession] = hmm
    print(f'read {len(amr_hmms)} HMMs')

    print('parse NCBI AMR hits...')
    hit_per_psc = {}
    with hmm_result_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (psc_id, _, hmm_name, hmm_id, evalue, bitscore, _) = re.split(r'\s+', line.strip(), maxsplit=6)
                hit = {
                    'psc_id': psc_id,
                    'hmm_id': hmm_id,
                    'bitscore': float(bitscore)
                }
                if(psc_id not in hit_per_psc):
                    hit_per_psc[psc_id] = hit
                else:
                    existing_hit = hit_per_psc[psc_id]
                    if(hit['bitscore'] > existing_hit['bitscore']):
                        hit_per_psc[psc_id] = hit
    print(f'read {len(hit_per_psc)} hits')

    for psc_id, hit in hit_per_psc.items():
        hmm = amr_hmms.get(hit['hmm_id'], None)
        if(hmm is not None):
            if(hmm['gene'] == ''):
                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (hmm['product'], psc_id))  # annotate PSC
                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', hmm['product'], psc_id)
                psc_annotated += 1
            else:
                conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (hmm['gene'], hmm['product'], psc_id))  # annotate PSC
                log_psc.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', hmm['gene'], hmm['product'], psc_id)
                psc_annotated += 1

print('\n')

print(f'PSCs with annotated AMR gene / product: {psc_annotated}')
log_psc.debug('summary: PSC annotated=%d', psc_annotated)

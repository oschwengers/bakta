import argparse
import logging
import re
import sqlite3

from pathlib import Path

from alive_progress import alive_bar

parser = argparse.ArgumentParser(
    description='Annotate PSCs by PFAM family descriptions.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--hmms', dest='hmms', action='store', help='Path to PFAM HMM file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to PFAM HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
hmms_path = Path(args.hmms).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - PFAM - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_psc = logging.getLogger('PSC')


pfam_descriptions = {}
PFAM_ACC_PATTERN = re.compile(r'^ACC\s{3}(PF.+)\n')
PFAM_DESC_PATTERN = re.compile(r'^DESC\s{2}(.+)\n')
print('parse PFAM entries...')
with hmms_path.open() as fh, alive_bar() as bar:
    acc = None
    for line in fh:
        if(acc is None):
            m = PFAM_ACC_PATTERN.fullmatch(line)
            if(m):
                acc = m.group(1)
        else:
            m = PFAM_DESC_PATTERN.fullmatch(line)
            if(m):
                description = m.group(1)
                pfam_descriptions[acc] = description
                acc = None
        bar()


print('parse Pfam hits...')
psc_ids = set()
best_hits = {}
with hmm_result_path.open() as fh, alive_bar() as bar:
    for line in fh:
        if(line[0] != '#'):
            cols = re.split(r'\s+', line.strip(), maxsplit=18)
            psc_id = cols[0]
            psc_ids.add(psc_id)
            pfam_acc = cols[3]
            bitscore = float(cols[5])
            best_hit = best_hits.get(psc_id, None)
            if(best_hit is None or bitscore > best_hit['bitscore']):
                best_hits[psc_id] = {
                    'psc_id': psc_id,
                    'pfam_acc': pfam_acc,
                    'bitscore': bitscore
                }
        bar()
print(f'PSC ids: {len(psc_ids)}')
print(f'best PFAM hits: {len(best_hits)}')
print('\n')

psc_annotated = 0
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, alive_bar(total=len(best_hits)) as bar:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()
    conn.row_factory = sqlite3.Row

    print('store best Pfam hits...')
    for psc_id, hit in best_hits.items():
        product = pfam_descriptions[hit['pfam_acc']]
        product_terms = product.split()
        if(len(product_terms) > 1):
            if(product_terms[-1] == 'domain'):
                product = f'{product}-containing protein'
        conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, psc_id))  # annotate PSC with Pfam family hit
        log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, psc_id)
        psc_annotated += 1
        if((psc_annotated % 1000000) == 0):
            conn.commit()
        bar()
print(f'PSCs annotated by Pfam family: {psc_annotated}')
log_psc.debug('summary: PSC annotated=%i', psc_annotated)

import argparse
import logging
import re
import sqlite3
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to NCBI AMR HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - PFAM - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
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

    conn.row_factory = sqlite3.Row

    print('parse Pfam hits...')
    with hmm_result_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = re.split('\s+', line.strip(), maxsplit=18)
                psc_id = cols[2]
                product = cols[18]
                product_terms = product.split()
                if(len(product_terms) > 1):
                    if(product_terms[-1] == 'domain'):
                        product = f'{product}-containing protein'
                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, psc_id))  # annotate PSC with Pfam family hit
                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, psc_id)
                psc_annotated += 1
            if((psc_annotated % 1000000) == 0):
                conn.commit()
                print(f'\t... {psc_annotated}')

print('\n')
print(f'PSCs annotated by Pfam family: {psc_annotated}')
log_ips.debug('summary: PSCs annotated=%d', psc_annotated)

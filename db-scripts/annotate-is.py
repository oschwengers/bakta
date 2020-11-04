import argparse
import logging
import sqlite3
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate PSCs by ISfinder transposon protein sequences.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--alignments', action='store', help='Path to diamond alignment file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
alignments_path = Path(args.alignments).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - IS - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


print('parse PSC alignments and add IS annotation...')
ups_processed = 0
ups_updated = 0
with alignments_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    for line in fh:
        ups_processed += 1
        (uniref90_id, sseqid, stitle, length, pident, qlen, slen, evalue) = line.strip().split('\t')
        length = int(length)
        qcov = length / int(qlen)
        scov = length / int(slen)
        if(qcov >= 0.99 and scov >= 0.99 and float(pident) >= 98. and float(evalue) < 1e-6):
            descriptions = stitle.split(' ', 1)[1][3:-3].split('~~~')
            if(descriptions[1].lower() == 'transposase'):
                descriptions = descriptions[0].split('_')
                is_name = descriptions[0]
                is_family = descriptions[2]
                product = f"{is_family}-like element {is_name} family transposase"
                gene = f'tnp-{is_name}'
                conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (gene, product, uniref90_id))
                log.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', gene, product, uniref90_id)
                ups_updated += 1
        if((ups_processed % 1000) == 0):
            conn.commit()
            print(f'\t... {ups_processed}')
    conn.commit()

print('\n')
print(f'PSCs processed: {ups_processed}')
print(f'PSCs with annotated IS transposons: {ups_updated}')
log.debug('summary: PSC annotated=%i', ups_updated)

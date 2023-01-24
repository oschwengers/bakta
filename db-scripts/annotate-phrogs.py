import argparse
import logging
import sqlite3

from pathlib import Path

from alive_progress import alive_bar

parser = argparse.ArgumentParser(
    description='Annotate PSCs by PHROG protein sequences.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--annotation', action='store', help='Path to PHROG annotation file.')
parser.add_argument('--psc-alignments', action='store', dest='psc_alignments', help='Path to PSC diamond alignment file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
annotation_path = Path(args.annotation).resolve()
psc_alignments_path = Path(args.psc_alignments).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - PHROG - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log = logging.getLogger('PSC')


phrogs = {}
with annotation_path.open() as fh_in:
    for line in fh_in:
        if not line.startswith('phrog'):
            (id, color, product, category) = line.strip().split('\t')
            if product != 'NA':
                phrogs[f'phrog_{id}'] = product


print('parse PSC alignments and add PHROG annotations...')
psc_processed = 0
psc_updated = 0
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    with psc_alignments_path.open() as fh, alive_bar() as bar:
        for line in fh:
            (uniref90_id, sseqid, stitle, length, pident, qlen, slen, evalue) = line.strip().split('\t')
            length = int(length)
            qcov = length / int(qlen)
            scov = length / int(slen)
            if(qcov >= 0.80 and scov >= 0.80 and float(pident) >= 90. and float(evalue) < 1e-6):
                if sseqid in phrogs:
                    product = phrogs[sseqid]
                    conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (product, uniref90_id))
                    log.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', product, uniref90_id)
                    psc_updated += 1
            psc_processed += 1
            if((psc_processed % 10000) == 0):
                conn.commit()
            bar()
        conn.commit()
    print(f'PSCs processed: {psc_processed}')
    print(f'PSCs annotated by PHROGs: {psc_updated}')
    log.debug('summary: PSC annotated=%i', psc_updated)
    

import argparse
import sqlite3
from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Extract hypothetical proteins for subsequent annotation.'
)
parser.add_argument('--psc', action='store', help='Path to PSC Fasta file.')
parser.add_argument('--db', action='store', help='Path to Bakta sqlite3 db file.')
parser.add_argument('--hypotheticals', action='store', help='Path to hypothetical protein Fasta file.')
args = parser.parse_args()


db_path = Path(args.db)
psc_path = Path(args.psc).resolve()
hypotheticals_path = Path(args.hypotheticals).resolve()

hypothetical_ids = set()

print('fetch hypothetical proteins (product is null) PSC IDs...')
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, alive_bar() as bar:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()
    conn.row_factory = sqlite3.Row
    for rec in conn.execute('SELECT uniref90_id FROM psc WHERE product IS NULL'):
        hypothetical_ids.add(rec['uniref90_id'])
        bar()
print(f'fetched {len(hypothetical_ids)} hypothetical IDs.')
print('\n')

print('extract sequences from psc.faa...')
hypotheticals = 0
with psc_path.open(mode='r') as fh_psc, hypotheticals_path.open(mode='w') as fh_hypothetical, alive_bar() as bar:
    for record in SeqIO.parse(fh_psc, 'fasta'):
        id = record.id
        seq = record.seq
        if(id in hypothetical_ids):
            fh_hypothetical.write(f'>{id}\n{seq}\n')
            hypotheticals += 1
        bar()
print(f'hypothetical proteins extracted: {hypotheticals}')

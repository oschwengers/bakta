
import argparse
import sqlite3
from pathlib import Path

parser = argparse.ArgumentParser(
    description='Create bakta sql db.'
)
parser.add_argument('--db', action='store', help='Path to bakta sqlite3 db file.')
args = parser.parse_args()

db_path = Path(args.db)

with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')

    print('create SQL table UPS...')
    conn.execute('DROP TABLE IF EXISTS ups;')
    conn.execute('''CREATE TABLE ups (
        hash TEXT PRIMARY KEY,
        length INTEGER NOT NULL,
        uniref100_id TEXT NOT NULL,
        uniref90_id TEXT,
        uniparc_id TEXT,
        ncbi_nrp_id TEXT,
        uniprotkb_acc TEXT,
        gene TEXT,
        product TEXT
        ) WITHOUT ROWID;''')
    conn.commit()
    print('\t...done')

    print('create SQL table PSC...')
    conn.execute('DROP TABLE IF EXISTS psc;')
    conn.execute('''CREATE TABLE psc (
        uniref90_id TEXT PRIMARY KEY,
        uniref50_id TEXT,
        gene TEXT,
        product TEXT,
        ec_id TEXT,
        cog_id TEXT,
        cog_category TEXT,
        go_ids TEXT
        ) WITHOUT ROWID;''')
    conn.commit()
    print('\t...done')

print("\nSQLite bakta db successfully created: %s" % db_path)

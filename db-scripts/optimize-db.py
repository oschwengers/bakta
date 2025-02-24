import argparse
import os
import sqlite3

from pathlib import Path

parser = argparse.ArgumentParser(
    description='Optimize fragmented DB.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file')
parser.add_argument('--tmp', action='store', default='.', help='Path to tmp directory')
args = parser.parse_args()


tmp_path = Path(args.tmp).resolve()
os.environ['SQLITE_TMPDIR'] = str(tmp_path)  # set SQLITE tmp dir env var
print(f"SQLITE tmp-dir: {os.environ['SQLITE_TMPDIR']}")

db_path = Path(args.db).resolve()
db_file_size_before = db_path.stat().st_size
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.execute(f"PRAGMA temp_store_directory = '{tmp_path}';")
    conn.commit()

    conn.execute('VACUUM')
    conn.commit()
db_file_size_after = db_path.stat().st_size
print('successfully optimized DB')
print(f'\tsize before: {db_file_size_before}')
print(f'\tsize after: {db_file_size_after}')

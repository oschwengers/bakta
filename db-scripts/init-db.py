
import argparse
import logging
import sqlite3
from pathlib import Path

parser = argparse.ArgumentParser(
    description='Create bakta sql db.'
)
parser.add_argument('--db', action='store', help='Path to bakta sqlite3 db file.')
args = parser.parse_args()


db_path = Path(args.db)


logging.basicConfig(
    filename='bakta.db.log',
    filemode='w',
    format='%(name)s - INIT - %(levelname)s - %(message)s',
    level=logging.INFO
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


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
    log_ups.info('DROP TABLE IF EXISTS ups;')
    stmt = '''CREATE TABLE ups (
        hash BLOB PRIMARY KEY,
        length INTEGER NOT NULL,
        uniparc_id TEXT,
        ncbi_nrp_id TEXT,
        uniref100_id TEXT
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_ups.info(stmt)
    conn.commit()
    print('\t...done')

    print('create SQL table IPS...')
    conn.execute('DROP TABLE IF EXISTS ips;')
    log_ips.info('DROP TABLE IF EXISTS ips;')
    stmt = '''CREATE TABLE ips (
        uniref100_id TEXT PRIMARY KEY,
        uniref90_id TEXT,
        gene TEXT,
        product TEXT,
        ec_ids TEXT,
        go_ids TEXT
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_ips.info(stmt)
    conn.commit()
    print('\t...done')

    print('create SQL table PSC...')
    log_psc.info('DROP TABLE IF EXISTS psc;')
    conn.execute('DROP TABLE IF EXISTS psc;')
    stmt = '''CREATE TABLE psc (
        uniref90_id TEXT PRIMARY KEY,
        uniref50_id TEXT,
        gene TEXT,
        product TEXT,
        ec_ids TEXT,
        cog_id TEXT,
        cog_category TEXT,
        go_ids TEXT
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_psc.info(stmt)
    conn.commit()
    print('\t...done')

print("\nSQLite bakta db successfully created: %s" % db_path)

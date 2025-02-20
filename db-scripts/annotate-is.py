import argparse
import logging
import sqlite3

from pathlib import Path

from alive_progress import alive_bar

parser = argparse.ArgumentParser(
    description='Annotate PSCs by ISfinder transposon protein sequences.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--ips-alignments', action='store', dest='ips_alignments', help='Path to IPS diamond alignment file.')
parser.add_argument('--psc-alignments', action='store', dest='psc_alignments', help='Path to PSC diamond alignment file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
ips_alignments_path = Path(args.ips_alignments).resolve()
psc_alignments_path = Path(args.psc_alignments).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - IS - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


print('parse PSC alignments and add IS annotation...')
ips_processed = 0
ips_updated = 0
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

    with ips_alignments_path.open() as fh, alive_bar(enrich_print=False) as bar:
        for line in fh:
            (uniref100_id, sseqid, stitle, length, pident, qlen, slen, evalue) = line.strip().split('\t')
            length = int(length)
            qcov = length / int(qlen)
            scov = length / int(slen)
            if(qcov >= 0.90 and scov >= 0.90 and float(pident) >= 95. and float(evalue) < 1e-6):
                (is_name, is_group, is_family, is_orf) = sseqid.split('_')
                gene = 'tnp'
                product = f"{is_family} family {is_name} transposase"
                if(is_orf == 'ORFA'):
                    product = f"{product} ORF A"
                elif(is_orf == 'ORFB'):
                    product = f"{product} ORF B"
                conn.execute('UPDATE ips SET gene=?, product=? WHERE uniref100_id=?', (gene, product, uniref100_id))
                log_ips.info('UPDATE ips SET gene=%s, product=%s WHERE uniref100_id=%s', gene, product, uniref100_id)
                ips_updated += 1
            ips_processed += 1
            if((ips_processed % 10000) == 0):
                conn.commit()
            bar()
        conn.commit()
    print(f'IPSs processed: {ips_processed}')
    print(f'IPSs with annotated IS transposases: {ips_updated}')
    print('\n')
    log_ips.debug('summary: IPS annotated=%i', ips_updated)


    with psc_alignments_path.open() as fh, alive_bar(enrich_print=False) as bar:
        for line in fh:
            (uniref90_id, sseqid, stitle, length, pident, qlen, slen, evalue) = line.strip().split('\t')
            length = int(length)
            qcov = length / int(qlen)
            scov = length / int(slen)
            if(qcov >= 0.80 and scov >= 0.80 and float(pident) >= 90. and float(evalue) < 1e-6):
                (is_name, is_group, is_family, is_orf) = sseqid.split('_')
                gene = 'tnp'
                product = f"{is_family} family transposase"
                if(is_orf == 'ORFA'):
                    product = f"{product} ORF A"
                elif(is_orf == 'ORFB'):
                    product = f"{product} ORF B"
                conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (gene, product, uniref90_id))
                log_psc.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', gene, product, uniref90_id)
                psc_updated += 1
            psc_processed += 1
            if((psc_processed % 10000) == 0):
                conn.commit()
            bar()
        conn.commit()
    print(f'PSCs processed: {psc_processed}')
    print(f'PSCs with annotated IS transposases: {psc_updated}')
    log_psc.debug('summary: PSC annotated=%i', psc_updated)
    

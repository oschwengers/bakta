import argparse
import logging
import hashlib
import sqlite3
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--amr', action='store', help='Path to NCBI Pathogen reference gene catalog file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
amr_path = Path(args.amr).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-AMR - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


print('lookup IPS by AMR NRP id (WP_*) and update gene/product annotations...')
nrps_processed = 0
ips_updated = 0
with amr_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    conn.execute('CREATE INDEX nrp ON ups ( ncbi_nrp_id )')
    conn.commit()
    log_ups.info('CREATE INDEX nrp ON ups ( ncbi_nrp_id )')

    conn.row_factory = sqlite3.Row
    for line in fh:
        nrps_processed += 1
        (allele, gene_family, whitelisted_taxa, product, scope, amr_type,
            subtype, clazz, subclazz, refseq_protein_accession,
            refseq_nucleotide_accession, curated_refseq_start,
            genbank_protein_accession, genbank_nucleotide_accession,
            genbank_strand_orientation, genbank_cds_start, genbank_cds_stop,
            pubmed_reference, blacklisted_taxa, db_version) = line.split('\t')
        if(refseq_protein_accession != '' and 'WP_' in refseq_protein_accession):
            refseq_protein_accession = refseq_protein_accession[3:]  # remove 'WP_' in NCBI NRP IDs
            if(scope == 'core'):
                if(amr_type == 'AMR' and subtype == 'AMR'):
                    rec_ups = conn.execute('SELECT * FROM ups WHERE ncbi_nrp_id=?', (refseq_protein_accession)).fetchone()
                    if(rec_ups is None):
                        log_ups.debug('no NRP hit: NRP-id=%s', refseq_protein_accession)
                        continue
                    log_ips.debug(
                        'nrp-id=%s, allele=%s, gene-family=%s, product=%s', 
                        refseq_protein_accession, allele, gene_family, product
                    )
                    gene = allele if allele != '' else gene_family
                    uniref100_id = rec_ups['uniref100_id']
                    if(gene == ''):
                        conn.execute('UPDATE ips SET product=? WHERE uniref100_id=?', (product, uniref100_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100_id
                        log_ips.info('UPDATE ips SET product=%s WHERE uniref100_id=%s', product, uniref100_id)
                        ips_updated += 1
                    else:
                        conn.execute('UPDATE ips SET gene=?, product=? WHERE uniref100_id=?', (gene, product, uniref100_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100
                        log_ips.info('UPDATE ips SET gene=%s, product=%s WHERE uniref100_id=%s', gene, product, uniref100_id)
                        ips_updated += 1
            elif(scope == 'plus'):
                pass
        if((nrps_processed % 100) == 0):
            conn.commit()
            print("\t... %d" % nrps_processed)
    conn.commit()

    conn.execute('DROP INDEX nrp')
    conn.commit()
    log_ups.info('DROP INDEX nrp')

print('\n')
print("NRPs processed: %d" % nrps_processed)
print("IPSs with annotated AMR gene / product: %d" % ips_updated)
log_ips.debug('summary: IPS annotated=%d', ips_updated)

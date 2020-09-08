import argparse
import logging
import hashlib
import sqlite3
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate UPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
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
log = logging.getLogger('UPS')


print('lookup UPS by AMR NRP id (WP_*) and update gene/product annotations...')
nrps_processed = 0
ups_updated = 0
with amr_path.open() as fh, sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute("PRAGMA mmap_size = %i;" % (20 * 1024 * 1024 * 1024))
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    for line in fh:
        nrps_processed += 1
        (allele, gene_family, whitelisted_taxa, product, scope, amr_type,
            subtype, clazz, subclazz, refseq_protein_accession,
            refseq_nucleotide_accession, curated_refseq_start,
            genbank_protein_accession, genbank_nucleotide_accession,
            genbank_strand_orientation, genbank_cds_start, genbank_cds_stop,
            pubmed_reference, blacklisted_taxa, db_version) = line.split('\t')
        try:
            if(refseq_protein_accession != '' and 'WP_' in refseq_protein_accession):
                refseq_protein_accession = refseq_protein_accession[3:]  # remove 'WP_' in NCBI NRP IDs
                if(scope == 'core'):
                    if(amr_type == 'AMR' and subtype == 'AMR'):
                        log.debug(
                            'nrp-id=%s, allele=%s, gene-family=%s, product=%s', 
                            refseq_protein_accession, allele, gene_family, product
                        )
                        gene = allele if allele != '' else gene_family
                        if(gene == ''):
                            conn.execute('UPDATE ups SET product=? WHERE ncbi_nrp_id=?', (product, refseq_protein_accession))  # annotate UPS with NCBI nrp id (WP_*)
                            log.info('UPDATE ups SET product=%s WHERE ncbi_nrp_id=%s', product, refseq_protein_accession)
                            ups_updated += 1
                        else:
                            conn.execute('UPDATE ups SET gene=?, product=? WHERE ncbi_nrp_id=?', (gene, product, refseq_protein_accession))  # annotate UPS with NCBI nrp id (WP_*)
                            log.info('UPDATE ups SET gene=%s, product=%s WHERE ncbi_nrp_id=%s', gene, product, refseq_protein_accession)
                            ups_updated += 1
                elif(scope == 'plus'):
                    pass
        except sqlite3.Error as e:
            print("SQL ERROR: could not set annotation for nrp-id=WP_%s" % refseq_protein_accession)
            print(e)
        if((nrps_processed % 100) == 0):
            conn.commit()
            print("\t... %d" % nrps_processed)
    conn.commit()

print('\n')
print("NRPs processed: %d" % nrps_processed)
print("UPSs with annotated AMR gene / product: %d" % ups_updated)
log.debug('summary: UPS annotated=%d', ups_updated)

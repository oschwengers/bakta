import argparse
import logging
import sqlite3

from pathlib import Path

from alive_progress import alive_bar


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--genes', action='store', help='Path to NCBI Pathogen reference gene catalog file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
genes_path = Path(args.genes).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-AMR - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')


nrps_processed = 0
ips_annotated = 0
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()

    print('create UPS index on NCBI NRP ids...')
    conn.execute('CREATE INDEX nrp ON ups ( ncbi_nrp_id )')
    conn.commit()
    log_ups.info('CREATE INDEX nrp ON ups ( ncbi_nrp_id )')

    print('lookup IPS by AMR NRP id (WP_*) and update gene/product annotations...')
    conn.row_factory = sqlite3.Row
    with genes_path.open() as fh, alive_bar() as bar:
        for line in fh:
            nrps_processed += 1
            (
                allele, gene_family, whitelisted_taxa, product, scope, amr_type,
                subtype, clazz, subclazz, refseq_protein_accession,
                refseq_nucleotide_accession, curated_refseq_start,
                genbank_protein_accession, genbank_nucleotide_accession,
                genbank_strand_orientation, genbank_cds_start, genbank_cds_stop,
                pubmed_reference, blacklisted_taxa, synonyms, hierarchy_node, db_version
            ) = line.split('\t')
            if('WP_' in refseq_protein_accession):
                refseq_protein_accession = refseq_protein_accession[3:]  # remove 'WP_' in NCBI NRP IDs
                rec_ups = conn.execute('SELECT * FROM ups WHERE ncbi_nrp_id=?', (refseq_protein_accession,)).fetchone()
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
                    conn.execute('UPDATE ips SET product=? WHERE uniref100_id=?', (product, uniref100_id))
                    log_ips.info('UPDATE ips SET product=%s WHERE uniref100_id=%s', product, uniref100_id)
                    ips_annotated += 1
                else:
                    conn.execute('UPDATE ips SET gene=?, product=? WHERE uniref100_id=?', (gene, product, uniref100_id))
                    log_ips.info('UPDATE ips SET gene=%s, product=%s WHERE uniref100_id=%s', gene, product, uniref100_id)
                    ips_annotated += 1
            if((nrps_processed % 1000) == 0):
                conn.commit()
            bar()
    conn.commit()

    print('drop UPS index on NCBI NRP ids...')
    conn.execute('DROP INDEX nrp')
    conn.commit()
    log_ups.info('DROP INDEX nrp')

print(f'NRPs processed: {nrps_processed}')
print(f'IPSs with annotated AMR gene / product: {ips_annotated}')
log_ips.debug('summary: IPS annotated=%d', ips_annotated)

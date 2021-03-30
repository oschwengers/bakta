import argparse
import logging
import re
import sqlite3
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Annotate IPSs by NCBI nrp IDs and PSCs by nrp cluster gene labels.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--genes', action='store', help='Path to NCBI Pathogen reference gene catalog file.')
parser.add_argument('--hmms', action='store', help='Path to NCBI AMR HMM description file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to NCBI AMR HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
genes_path = Path(args.genes).resolve()
hmms_path = Path(args.hmms).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-AMR - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_ups = logging.getLogger('UPS')
log_ips = logging.getLogger('IPS')
log_psc = logging.getLogger('PSC')


nrps_processed = 0
ips_updated = 0
psc_annotated = set()
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
    with genes_path.open() as fh:
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
                print(f'\t... {nrps_processed}')
    conn.commit()

    print('drop UPS index on NCBI NRP ids...')
    conn.execute('DROP INDEX nrp')
    conn.commit()
    log_ups.info('DROP INDEX nrp')

    print('parse NCBI AMR HMMs...')
    amr_hmms = {}
    with hmms_path.open() as fh:
        for line in fh:
            (accession, short_name, TC1, TC2, gene_symbol, length, family_type, curated, product_name,
                AMRFinder_reportable, db_type, subtype, clazz, subclass) = line.split('\t')
            hmm = {
                'acc': accession,
                'gene': gene_symbol,
                'type': family_type,
                'product': product_name
            }
            amr_hmms[accession] = hmm

    print('parse NCBI AMR hits...')
    exception_hits_per_psc = set()
    equivalog_hits_per_psc = set()
    with hmm_result_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (psc_id, _, hmm_name, hmm_id, evalue, bitscore, _) = re.split('\s+', line.strip(), maxsplit=6)
                hit = {
                    'psc_id': psc_id,
                    'hmm_id': hmm_id,
                    'evalue': float(evalue),
                    'bitscore': float(bitscore)
                }
                hmm = amr_hmms.get(hit['hmm_id'], None)
                if(hmm is not None):
                    if(hmm['type'] == 'exception'):
                        if(hmm['gene'] == ''):
                            conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100_id
                            log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', hmm['product'], psc_id)
                            psc_annotated.add(psc_id)
                        else:
                            conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (hmm['gene'], hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100
                            log_psc.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', hmm['gene'], hmm['product'], psc_id)
                            psc_annotated.add(psc_id)
                        exception_hits_per_psc.add(psc_id)
                    elif(hmm['type'] == 'equivalog'):
                        if(psc_id not in exception_hits_per_psc):
                            if(hmm['gene'] == ''):
                                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100_id
                                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', hmm['product'], psc_id)
                                psc_annotated.add(psc_id)
                            else:
                                conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (hmm['gene'], hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100
                                log_psc.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', hmm['gene'], hmm['product'], psc_id)
                                psc_annotated.add(psc_id)
                            equivalog_hits_per_psc.add(psc_id)
                    elif(hmm['type'] == 'subfamily'):
                        if(psc_id not in exception_hits_per_psc and psc_id not in equivalog_hits_per_psc):
                            if(hmm['gene'] == ''):
                                conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100_id
                                log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', hmm['product'], psc_id)
                                psc_annotated.add(psc_id)
                            else:
                                conn.execute('UPDATE psc SET gene=?, product=? WHERE uniref90_id=?', (hmm['gene'], hmm['product'], psc_id))  # annotate IPS with NCBI nrp id (WP_*) -> UniRef100
                                log_psc.info('UPDATE psc SET gene=%s, product=%s WHERE uniref90_id=%s', hmm['gene'], hmm['product'], psc_id)
                                psc_annotated.add(psc_id)

print('\n')
print(f'NRPs processed: {nrps_processed}')
print(f'IPSs with annotated AMR gene / product: {ips_updated}')
log_ips.debug('summary: IPS annotated=%d', ips_updated)

psc_updated = len(psc_annotated)
print(f'PSCs with annotated AMR gene / product: {psc_updated}')
log_ips.debug('summary: PSCs annotated=%d', psc_updated)

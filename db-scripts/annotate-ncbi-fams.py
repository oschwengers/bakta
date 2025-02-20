import argparse
import logging
import re
import sqlite3

from pathlib import Path

from alive_progress import alive_bar


parser = argparse.ArgumentParser(
    description='Annotate PSCs by NCBIfams.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--hmms', action='store', help='Path to NCBIfams HMM description file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to NCBIfams HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
hmms_path = Path(args.hmms).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - NCBI-AMR - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_psc = logging.getLogger('PSC')

family_type_ranks = {
    'exception': 77,
    'equivalog': 70
}

print('parse NCBIfam HMMs...')
hmms = {}
with hmms_path.open() as fh, alive_bar(enrich_print=False) as bar:
    for line in fh:
        (
            ncbi_accession,
            source_identifier,
            label,
            sequence_cutoff,
            domain_cutoff,
            hmm_length,
            family_type,
            for_structural_annotation,
            for_naming,
            for_AMRFinder,
            product_name,
            gene_symbol,
            gene_synonyms,
            ec_numbers,
            go_terms,
            pmids,
            taxonomic_range,
            taxonomic_range_name,
            taxonomic_rank_name,
            n_refseq_protein_hits,
            source,
            name_orig,
            hmm_name,
            comment
        ) = line.split('\t')
        if(family_type in family_type_ranks.keys() and for_naming == 'Y'):  # only accept exception/equivalog HMMs eligible for naming proteins
            hmm = {
                'acc': ncbi_accession,
                'type': family_type,
                'product': product_name.strip()
            }
            if(gene_symbol != ''):
                hmm['gene'] = gene_symbol
            if(ec_numbers != ''):
                hmm['ecs'] = ec_numbers.split(',')
            if(go_terms != ''):
                hmm['gos'] = [term.split(':')[1] for term in go_terms.split(',')]
            if(ncbi_accession != '-'):
                hmms[ncbi_accession] = hmm
        bar()
print(f'read {len(hmms)} HMMs')
print('\n')

print('parse NCBIfam hits...')
hit_per_psc = {}
with hmm_result_path.open() as fh, alive_bar(enrich_print=False) as bar:
    for line in fh:
        if(line[0] != '#'):
            (psc_id, _, hmm_name, hmm_id, evalue, bitscore, _) = re.split(r'\s+', line.strip(), maxsplit=6)
            hit = {
                'psc_id': psc_id,
                'hmm_id': hmm_id,
                'bitscore': float(bitscore)
            }
            if(psc_id not in hit_per_psc):
                hit_per_psc[psc_id] = hit
            else:
                hit_hmm = hmms.get( hit['hmm_id'], {} )
                hit_rank = family_type_ranks.get( hit_hmm.get('type', ''), -1 )
                existing_hit = hit_per_psc[psc_id]
                existing_hit_hmm = hmms.get( existing_hit['hmm_id'], {} )
                existing_hit_rank = family_type_ranks.get( existing_hit_hmm.get('type', ''), -1 )
                if(hit_rank < 0):
                    continue
                elif(hit_rank > existing_hit_rank):  # give precedence to HMMs of higher ranking family types
                    hit_per_psc[psc_id] = hit
                elif(hit_rank == existing_hit_rank):
                    if(hit['bitscore'] > existing_hit['bitscore']):
                        hit_per_psc[psc_id] = hit
        bar()
print(f'read {len(hit_per_psc)} hits')
print('\n')

psc_annotated = 0
psc_with_gene = 0
psc_with_ec = 0
psc_with_go = 0
with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, alive_bar(total=len(hit_per_psc), enrich_print=False) as bar:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()
    conn.row_factory = sqlite3.Row
    for psc_id, hit in hit_per_psc.items():
        hmm = hmms.get(hit['hmm_id'], None)
        if(hmm is not None):
            conn.execute('UPDATE psc SET product=? WHERE uniref90_id=?', (hmm['product'], psc_id))  # annotate PSC
            log_psc.info('UPDATE psc SET product=%s WHERE uniref90_id=%s', hmm['product'], psc_id)
            if('gene' in hmm):
                conn.execute('UPDATE psc SET gene=? WHERE uniref90_id=?', (hmm['gene'], psc_id))  # annotate PSC
                log_psc.info('UPDATE psc SET gene=%s WHERE uniref90_id=%s', hmm['gene'], psc_id)
                psc_with_gene += 1
            if('ecs' in hmm):
                rec_psc = conn.execute('SELECT ec_ids FROM psc WHERE uniref90_id=?', (psc_id,)).fetchone()
                if(rec_psc is None):
                    continue
                if(rec_psc['ec_ids'] is not None):
                    ecs = set(rec_psc['ec_ids'].split(','))
                    for ec in hmm['ecs']:
                        ecs.add(ec)
                else:
                    ecs = hmm['ecs']
                ec_list = ','.join(sorted(ecs))
                conn.execute('UPDATE psc SET ec_ids=? WHERE uniref90_id=?', (ec_list, psc_id))  # annotate PSC
                log_psc.info('UPDATE psc SET ec_ids=%s WHERE uniref90_id=%s', ec_list, psc_id)
                psc_with_ec += 1
            if('gos' in hmm):
                rec_psc = conn.execute('SELECT go_ids FROM psc WHERE uniref90_id=?', (psc_id,)).fetchone()
                if(rec_psc is None):
                    continue
                if(rec_psc['go_ids'] is not None):
                    gos = set(rec_psc['go_ids'].split(','))
                    for go in hmm['gos']:
                        gos.add(go)
                else:
                    gos = hmm['gos']
                go_list = ','.join(sorted(gos))
                conn.execute('UPDATE psc SET go_ids=? WHERE uniref90_id=?', (go_list, psc_id))  # annotate PSC
                log_psc.info('UPDATE psc SET go_ids=%s WHERE uniref90_id=%s', go_list, psc_id)
                psc_with_go += 1
            psc_annotated += 1
        bar()
print(f'PSCs with annotated product: {psc_annotated}')
log_psc.debug('summary: PSC annotated=%d', psc_annotated)
print(f'PSCs with annotated gene: {psc_with_gene}')
log_psc.debug('summary: PSC gene annotated=%d', psc_with_gene)
print(f'PSCs with annotated EC: {psc_with_ec}')
log_psc.debug('summary: PSC EC annotated=%d', psc_with_ec)
print(f'PSCs with annotated GO: {psc_with_go}')
log_psc.debug('summary: PSC GO annotated=%d', psc_with_go)

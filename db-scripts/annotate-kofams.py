import argparse
import logging
import re
import sqlite3

from pathlib import Path

from alive_progress import alive_bar

RE_EC = re.compile(r'(?:\[EC:([0-9\.\-\s]+)\])')

parser = argparse.ArgumentParser(
    description='Annotate PSCs by kofams.'
)
parser.add_argument('--db', action='store', help='Path to Bakta db file.')
parser.add_argument('--hmms', action='store', help='Path to KEGG kofams HMM description file.')
parser.add_argument('--hmm-results', dest='hmm_results', action='store', help='Path to kofams HMM output file.')
args = parser.parse_args()


db_path = Path(args.db).resolve()
hmms_path = Path(args.hmms).resolve()
hmm_result_path = Path(args.hmm_results).resolve()


logging.basicConfig(
    filename='bakta.db.log',
    filemode='a',
    format='%(name)s - KOFAMS - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_psc = logging.getLogger('PSC')


print('parse kofams HMMs...')
hmms = {}
with hmms_path.open() as fh:
    skip = True
    for line in fh:
        if(skip):
            skip = False
            continue
        (
            knum, threshold, score_type, profile_type, f_measure, nseq, nseq_used, alen, mlen, eff_nseq, re_pos, definition
        ) = line.split('\t')
        knum = knum
        hmm = {
            'ko': knum[1:],
            'threshold': float(threshold),
            'score_type': score_type,
            'profile_type': profile_type,
            'f_measure': float(f_measure),
            'product': definition.strip()
        }
        m = RE_EC.search(definition)
        if(m):
            hmm['ecs'] = m[1].split()
        if(knum != '-' and hmm['f_measure'] > 0.77):  # discard the lower 10th percentile
            hmms[knum] = hmm
print(f'read {len(hmms)} HMMs')

print('parse Kofam hits...')
hit_per_psc = {}
with hmm_result_path.open() as fh, alive_bar(enrich_print=False) as bar:
    for line in fh:
        if(line[0] != '#'):
            (psc_id, _, hmm_name, hmm_id, evalue, bitscore, _) = re.split(r'\s+', line.strip(), maxsplit=6)
            hit = {
                'psc_id': psc_id,
                'hmm_id': hmm_name,
                'bitscore': float(bitscore)
            }
            hmm = hmms.get(hmm_name, None)
            if(hmm is not None and hit['bitscore'] > hmm['threshold']):
                if(psc_id not in hit_per_psc):
                    hit_per_psc[psc_id] = hit
                else:
                    existing_hit = hit_per_psc[psc_id]
                    if(hit['bitscore'] > existing_hit['bitscore']):
                        hit_per_psc[psc_id] = hit
        bar()
print(f'parsed and selected {len(hit_per_psc)} valid hits')
print('\n')

psc_annotated = 0
ecs_added = 0
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
            conn.execute('UPDATE psc SET kegg_orthology_id=? WHERE uniref90_id=?', (hmm['ko'], psc_id))  # annotate PSC
            log_psc.info('UPDATE psc SET kegg_orthology_id=%s WHERE uniref90_id=%s', hmm['ko'], psc_id)
            psc_annotated += 1
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
                ecs_added += 1
        bar()
print(f'PSCs with annotated kofam: {psc_annotated}')
log_psc.debug('summary: PSC annotated=%d', psc_annotated)
print(f'Added ECs: {ecs_added}')

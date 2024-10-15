import logging

from typing import Sequence

import pyhmmer

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf


log = logging.getLogger('EXPERT_AA_HMM')


def search(cdss: Sequence[dict], user_hmms_path):
    """Detect Pfam-A entries"""
    cds_found = set()
    orf_by_aa_digest = orf.get_orf_dictionary(cdss)
    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    proteins: list[pyhmmer.easel.DigitalSequence] = [ pyhmmer.easel.TextSequence(sequence=cds['aa'], name=bytes(orf.get_orf_key(cds), 'UTF-8')).digitize(alphabet) for cds in cdss ]
    with pyhmmer.plan7.HMMFile(user_hmms_path) as hmms_fh:
        hmms = list(hmms_fh)
        for hmm_query_hits in pyhmmer.hmmsearch(hmms, proteins, bit_cutoffs='trusted', cpus=cfg.threads):
            hmm_id = hmm_query_hits.query.accession.decode()
            hmm_length = hmm_query_hits.query.M
            hmm_description = hmm_query_hits.query.description.decode()
            hmms_description_fields = hmm_description.split('~~~')
            for hmm_query_hit in hmm_query_hits.reported:
                aa_identifier = hmm_query_hit.name.decode()
                cds = orf_by_aa_digest[aa_identifier]
                if hmm_query_hit.evalue > bc.MIN_HMM_EVALUE:
                    log.debug(
                        'discard low evalue: seq=%s, start=%i, stop=%i, strand=%s, id=%s, evalue=%1.1e, bitscore=%f',
                        cds['sequence'], cds['start'], cds['stop'], cds['strand'], hmm_id, hmm_query_hit.evalue, hmm_query_hit.score
                    )
                else:
                    hit_domain_lengths_sum = sum([len(dom.alignment.hmm_sequence) for dom in hmm_query_hit.domains.included])
                    hit = {
                        'type': 'user_hmms',
                        'rank': 100,
                        'id': hmm_id,
                        'length': hit_domain_lengths_sum,
                        'aa_cov': hit_domain_lengths_sum / len(cds['aa']),
                        'hmm_cov': hit_domain_lengths_sum / hmm_length,
                        'evalue': hmm_query_hit.evalue,
                        'score': hmm_query_hit.score,
                        'start': hmm_query_hit.best_domain.alignment.target_from,
                        'stop': hmm_query_hit.best_domain.alignment.target_to,
                        'db_xrefs': [f'UserHMM:{hmm_id}'],
                        'gene': None,
                        'product': None
                    }
                    # read user-provided gene, product, db_xrefs via HMM's description field provided as <gene>~~~<product>~~~<db_xrefs>
                    if(len(hmms_description_fields) == 1):
                        hit['product'] = hmm_description if hmm_description != '' else None
                    elif(len(hmms_description_fields) == 3):
                        (gene, product, dbxrefs) = hmms_description_fields
                        hit['gene'] = gene
                        hit['product'] = product
                        hit['db_xrefs'].extend(set(dbxrefs.split(',')))
                    else:
                        log.warning('wrong HMM description format: id=%s, desc=%s', hmm_id, hmm_description)

                    cds.setdefault('expert', [])
                    cds['expert'].append(hit)
                    log.debug(
                        'hit: source=UserHMMs, rank=99, seq=%s, start=%i, stop=%i, strand=%s, query-cov=%0.3f, model-cov=%0.3f, hmm-id=%s, gene=%s, product=%s, evalue=%1.1e, bitscore=%f',
                        cds['sequence'], cds['start'], cds['stop'], cds['strand'], hit['aa_cov'], hit['hmm_cov'], hmm_id, hit['gene'], hit['product'], hit['evalue'], hit['score']
                    )
                    cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found

import logging

from typing import Sequence

import pyhmmer

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf


log = logging.getLogger('EXPERT_AA_HMM')


def search(cdss: Sequence[dict], user_hmms):
    """Detect Pfam-A entries"""
    cds_found = set()
    orf_by_aa_digest = orf.get_orf_dictionary(cdss)
    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    proteins: list[pyhmmer.easel.DigitalSequence] = [ pyhmmer.easel.TextSequence(sequence=cds['aa'], name=bytes(orf.get_orf_key(cds), 'UTF-8')).digitize(alphabet) for cds in cdss ]
    with pyhmmer.plan7.HMMFile(user_hmms) as hmm:
        for hmm_query_hits in pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs='trusted', cpus=cfg.threads):
            print(f'TopHits: {hmm_query_hits}')
            hmm_name = hmm_query_hits.query_name.decode()
            print(f'query_name: {hmm_name}')
            hmm_id = hmm_query_hits.query_accession.decode()
            print(f'query_accession: {hmm_id}')
            hmm_length = hmm_query_hits.query_length
            print(f'query_length {hmm_length}')
            # print(f'description: {hmm_query_hits.description.decode()}')
            
            for hmm_query_hit in hmm_query_hits.reported:
                print(f'Hit: {hmm_query_hit}')
                # print(f'\taccession: {hmm_query_hit.accession.decode()}')
                # print(f'\tdescription: {hmm_query_hit.description.decode()}')
                print(f'\tname: {hmm_query_hit.name.decode()}')
                print(f'\tevalue: {hmm_query_hit.evalue}')
                print(f'\tlength: {hmm_query_hit.length}')
                print(f'\tscore: {hmm_query_hit.score}')

                aa_identifier = hmm_query_hit.name.decode()
                cds = orf_by_aa_digest[aa_identifier]
                if hmm_query_hit.evalue > bc.MIN_HMM_EVALUE:
                    log.debug(
                        'discard low evalue: contig=%s, start=%i, stop=%i, strand=%s, id=%s, evalue=%1.1e, bitscore=%f',
                        cds['contig'], cds['start'], cds['stop'], cds['strand'], hmm_id, hmm_query_hit.evalue, hmm_query_hit.score
                    )
                else:
                    hit_domain_lengths_sum = sum([len(dom.alignment.hmm_sequence) for dom in hmm_query_hit.domains.included])
                    print(f'hit domain length sum: {hit_domain_lengths_sum}')
                    hmm_cov = hit_domain_lengths_sum / hmm_length
                    print(f'domain cov: {hmm_cov}')
                    aa_cov = hit_domain_lengths_sum / hmm_length
                    aa_cov_best_domain = (hmm_query_hit.best_domain.alignment.target_to - hmm_query_hit.best_domain.alignment.target_from + 1) / len(cds['aa'])
                    assert aa_cov == aa_cov_best_domain, f'ERROR: aa_cov ({aa_cov}) != aa_cov_best_domain ({aa_cov_best_domain})'
                    print(f'aa cov: {aa_cov}')

                    hit = {
                        'type': 'user_hmms',
                        'source': 'UserHMMs',
                        'rank': 100,
                        'id': hmm_id,
                        'length': hit_domain_lengths_sum,
                        'aa_cov': aa_cov,
                        'hmm_cov': hmm_cov,
                        'evalue': hmm_query_hit.evalue,
                        'score': hmm_query_hit.score,
                        'start': hmm_query_hit.best_domain.alignment.target_from,
                        'stop': hmm_query_hit.best_domain.alignment.target_to,
                        'db_xrefs': [f'UserHMM:{hmm_id}'],
                        'gene': None,
                        'product': hmm_name  # None
                    }

                    # ToDo: implement user-provided gene, product, db_xrefs via HMM's description field provided as <gene>~~~<product>~~~<db_xrefs>, waiting for https://github.com/althonos/pyhmmer/issues/76
                    # hmm_description = hmm_query_hits.description.decode()
                    # cols = hmm_description.split('~~~')
                    # if(len(cols) == 1):
                    #     hit['product'] = hmm_description
                    # else:
                    #     (gene, product, dbxrefs) = cols
                    #     hit['gene'] = gene
                    #     hit['product'] = product
                    #     hit['db_xrefs'].extend(set(dbxrefs.split(',')))

                    cds.setdefault('expert', [])
                    cds['expert'].append(hit)
                    log.debug(
                        'hit: source=UserHMMs, rank=99, contig=%s, start=%i, stop=%i, strand=%s, query-cov=%0.3f, model-cov=%0.3f, gene=%s, product=%s, evalue=%1.1e, bitscore=%f',
                        cds['contig'], cds['start'], cds['stop'], cds['strand'], hit['aa_cov'], hit['hmm_cov'], hit['gene'], hit['product'], hit['evalue'], hit['score']
                    )
                    cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found

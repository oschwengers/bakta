import logging
import subprocess as sp

from pathlib import Path
from typing import Sequence

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf


log = logging.getLogger('EXPERT-AMRFINDER')


def search(cdss: Sequence[dict], cds_fasta_path: Path):
    """Conduct expert CDS analysis with AMRFinderPlus."""
    amrfinder_output_path = cfg.tmp_path.joinpath('amrfinder.tsv')
    amrfinderplus_db_path = cfg.db_path.joinpath('amrfinderplus-db')
    amrfinderplus_db_latest_path = amrfinderplus_db_path.joinpath('latest')

    amrfinderplus_tmp_path = cfg.tmp_path.joinpath('amrfinderplus')
    amrfinderplus_tmp_path.mkdir()
    env = cfg.env.copy()
    env['TMPDIR'] = str(amrfinderplus_tmp_path)

    cmd = [
        'amrfinder',
        '--database', str(amrfinderplus_db_latest_path),
        '--protein', str(cds_fasta_path),
        '--plus',
        '--translation_table', str(cfg.translation_table),
        '--output', str(amrfinder_output_path),
        '--threads', str(cfg.threads)
    ]
    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('AMR expert system failed! amrfinder-error-code=%d', proc.returncode)
        raise Exception(f"amrfinder error! error code: {proc.returncode}. Please, try 'amrfinder_update --force_update --database {amrfinderplus_db_path}' to update AMRFinderPlus's internal database.")

    cds_found = set()
    cds_by_hexdigest = orf.get_orf_dictionary(cdss)
    with amrfinder_output_path.open() as fh:
        for line in fh:
            if(line[:7] != 'Protein'):
                (
                    aa_identifier, gene, product, scope, element_type, element_subtype, clazz, subclass, method, target_length, reference_sequence_length,
                    cov_ref_seq, ident_ref_seq, alignment_length, accession_closest_seq, name_closest_seq, hmm_id, hmm_description
                ) = line.split('\t')
                cds = cds_by_hexdigest[aa_identifier]
                hit = {
                    'type': 'amrfinder',
                    'rank': 95,
                    'gene': gene if gene != '' else None,
                    'product': product,
                    'method': method
                }
                if(method.lower() != 'hmm'):
                    hit['query_cov'] = int(alignment_length) / len(cds['aa'])
                    model_cov = float(cov_ref_seq) / 100
                    hit['model_cov'] = model_cov
                    identity = float(ident_ref_seq) / 100
                    hit['identity'] = identity
                    hit['id'] = accession_closest_seq
                    hit['db_xrefs'] = [f'{bc.DB_XREF_NCBI_PROTEIN}:{accession_closest_seq}']
                else:
                    model_cov = 0
                    identity = 0
                    hit['id'] = hmm_id
                    hit['db_xrefs'] = [f'{bc.DB_XREF_NCBI_FAMILIES}:{hmm_id}']

                cds.setdefault('expert', [])
                cds['expert'].append(hit)
                log.debug(
                    'hit: gene=%s, product=%s, method=%s, target-cov=%0.3f, identity=%0.3f, seq=%s, start=%i, stop=%i, strand=%s',
                    gene, product, method, model_cov, identity, cds['sequence'], cds['start'], cds['stop'], cds['strand']
                )
                cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found

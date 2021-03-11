import logging
from concurrent.futures import ThreadPoolExecutor
import subprocess as sp
import sqlite3

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('EXPERT-AMRFINDER')

def search(cdss):
    """Conduct expert CDS analysis with AMRFinderPlus."""
    cds_fasta_path = cfg.tmp_path.joinpath('cds.faa')
    amrfinder_output_path = cfg.tmp_path.joinpath('amrfinder.tsv')
    cmd = [
        'amrfinder',
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
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('AMR expert system failed! amrfinder-error-code=%d', proc.returncode)
        raise Exception(f'amrfinder error! error code: {proc.returncode}')

    cds_by_hexdigest = {cds['aa_hexdigest']: cds for cds in cdss}
    with amrfinder_output_path.open() as fh:
        for line in fh:
            if(line[:7] != 'Protein'):
                (aa_hexdigest, gene, product, scope, element_type, element_subtype, clazz, subclass, method, target_length, reference_sequence_length,
                cov_ref_seq, ident_ref_seq, alignment_length, accession_closest_seq, name_closest_seq, hmm_id, hmm_description) = line.split('\t')
                target_cov = float(cov_ref_seq) / 100
                identity = float(ident_ref_seq) / 100
                cds = cds_by_hexdigest[aa_hexdigest]
                if('expert' not in cds):
                    cds['expert'] = {}
                cds['expert']['amrfinder'] = {
                    'gene': gene,
                    'product': product,
                    'target-cov': target_cov,
                    'identity': identity
                }
                log.debug(
                    'hit: gene=%s, product=%s, target-cov=%0.3f, identity=%0.3f, contig=%s, start=%i, stop=%i, strand=%s',
                    gene, product, target_cov, identity, cds['contig'], cds['start'], cds['stop'], cds['strand']
                )

    pscs_found = []
    pscs_not_found = []
    for cds in cdss:
        if('psc' in cds):
            pscs_found.append(cds)
        else:
            pscs_not_found.append(cds)
    log.info('found=%i', len(pscs_found))
    return pscs_found, pscs_not_found

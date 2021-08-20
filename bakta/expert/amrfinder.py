import logging
import subprocess as sp

import bakta.config as cfg


log = logging.getLogger('EXPERT-AMRFINDER')


def search(cdss, cds_fasta_path):
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
    cds_by_hexdigest = {f"{cds['aa_hexdigest']}-{cds['contig']}-{cds['start']}": cds for cds in cdss}
    with amrfinder_output_path.open() as fh:
        for line in fh:
            if(line[:7] != 'Protein'):
                (
                    aa_identifier, gene, product, scope, element_type, element_subtype, clazz, subclass, method, target_length, reference_sequence_length,
                    cov_ref_seq, ident_ref_seq, alignment_length, accession_closest_seq, name_closest_seq, hmm_id, hmm_description
                ) = line.split('\t')
                cds = cds_by_hexdigest[aa_identifier]
                query_cov = int(alignment_length) / len(cds['aa'])
                model_cov = float(cov_ref_seq) / 100
                identity = float(ident_ref_seq) / 100
                expert_hit = {
                    'system': 'amrfinder',
                    'rank': 95,
                    'gene': gene if gene != '' else None,
                    'product': product,
                    'query_cov': query_cov,
                    'model_cov': model_cov,
                    'identity': identity
                }
                if('expert' not in cds):
                    cds['expert'] = {}
                cds['expert']['amrfinder'] = expert_hit
                log.debug(
                    'hit: gene=%s, product=%s, target-cov=%0.3f, identity=%0.3f, contig=%s, start=%i, stop=%i, strand=%s',
                    gene, product, model_cov, identity, cds['contig'], cds['start'], cds['stop'], cds['strand']
                )
                cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found

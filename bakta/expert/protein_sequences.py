import logging
import subprocess as sp

import bakta.config as cfg


log = logging.getLogger('EXPERT_AA_SEQ')


def search(cdss, cds_fasta_path):
    """Conduct homology search of CDSs against PCS db."""
    diamond_db_path = cfg.db_path.joinpath('expert-protein-sequences.dmnd')
    diamond_output_path = cfg.tmp_path.joinpath('diamond.cds.expert.tsv')
    cmd = [
        'diamond',
        'blastp',
        '--query', str(cds_fasta_path),
        '--db', str(diamond_db_path),
        '--out', str(diamond_output_path),
        '--id', str(50),  # '50',
        '--query-cover', str(80),  # '80'
        '--subject-cover', str(80),  # '80'
        '--max-target-seqs', '5',  # single best output
        '--outfmt', '6', 'qseqid', 'sseqid', 'slen', 'length', 'pident', 'evalue', 'bitscore', 'stitle',
        '--threads', str(cfg.threads),
        '--tmpdir', str(cfg.tmp_path),  # use tmp folder
        '--block-size', '4',  # increase block size for faster executions
        '--index-chunks', '1'  # set index chunks to 1 for faster executions
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
        log.warning('diamond failed! diamond-error-code=%d', proc.returncode)
        raise Exception(f'diamond error! error code: {proc.returncode}')

    cds_found = set()
    cds_by_hexdigest = {f"{cds['aa_hexdigest']}-{cds['contig']}-{cds['start']}": cds for cds in cdss}
    with diamond_output_path.open() as fh:
        for line in fh:
            (aa_identifier, model_id, model_length, alignment_length, identity, evalue, bitscore, model_title) = line.strip().split('\t')
            cds = cds_by_hexdigest[aa_identifier]
            query_cov = int(alignment_length) / len(cds['aa'])
            model_cov = int(alignment_length) / int(model_length)
            identity = float(identity) / 100
            evalue = float(evalue)
            bitscore = float(bitscore)
            (expert_system, rank, min_identity, min_query_cov, min_model_cov, gene, product, dbxrefs) = model_title.split(' ', 1)[1].split('~~~')
            rank = int(rank)
            min_identity = float(min_identity) / 100
            min_query_cov = float(min_query_cov) / 100
            min_model_cov = float(min_model_cov) / 100
            if(query_cov >= min_query_cov and model_cov >= min_model_cov and identity >= min_identity):
                expert_hit = {
                    'system': expert_system,
                    'rank': rank,
                    'gene': gene if gene != '' else None,
                    'product': product,
                    'query_cov': query_cov,
                    'model_cov': model_cov,
                    'identity': identity,
                    'evalue': evalue,
                    'bitscore': bitscore
                }
                dbxrefs = None if dbxrefs == '' else dbxrefs.split(',')
                if(dbxrefs is not None):
                    expert_hit['db_xrefs'] = dbxrefs
                if('expert' not in cds):
                    cds['expert'] = {}
                cds['expert']['aa_seq'] = expert_hit
                log.debug(
                    'hit: contig=%s, start=%i, stop=%i, strand=%s, system=%s, rank=%i, query-cov=%0.3f, model-cov=%0.3f, identity=%0.3f, gene=%s, product=%s, evalue=%1.1e, bitscore=%f',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], expert_system, rank, query_cov, model_cov, identity, gene, product, evalue, bitscore
                )
                cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found

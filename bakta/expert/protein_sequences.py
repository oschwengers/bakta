import logging
import subprocess as sp
import sys

import bakta.config as cfg

from Bio import SeqIO
from xopen import xopen


log = logging.getLogger('EXPERT_AA_SEQ')


def search(cdss, cds_fasta_path, expert_system, db_path):
    """Conduct homology search of CDSs against PCS db."""
    diamond_output_path = cfg.tmp_path.joinpath('diamond.cds.expert.tsv')
    cmd = [
        'diamond',
        'blastp',
        '--query', str(cds_fasta_path),
        '--db', str(db_path),
        '--out', str(diamond_output_path),
        '--id', str(50),  # '50',
        '--query-cover', str(80),  # '80'
        '--subject-cover', str(80),  # '80'
        '--max-target-seqs', '1',  # single best output
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
            (source, rank, min_identity, min_query_cov, min_model_cov, gene, product, dbxrefs) = model_title.split(' ', 1)[1].split('~~~')
            rank = int(rank)
            min_identity = float(min_identity) / 100
            min_query_cov = float(min_query_cov) / 100
            min_model_cov = float(min_model_cov) / 100
            if(query_cov >= min_query_cov and model_cov >= min_model_cov and identity >= min_identity):
                expert_hit = {
                    'source': source,
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
                cds['expert'][expert_system] = expert_hit
                log.debug(
                    'hit: contig=%s, start=%i, stop=%i, strand=%s, source=%sm, rank=%i, query-cov=%0.3f, model-cov=%0.3f, identity=%0.3f, gene=%s, product=%s, evalue=%1.1e, bitscore=%f',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], source, rank, query_cov, model_cov, identity, gene, product, evalue, bitscore
                )
                cds_found.add(aa_identifier)

    log.info('found=%i', len(cds_found))
    return cds_found


def write_user_protein_sequences(aa_fasta_path):
    try:
        with xopen(str(cfg.user_proteins), threads=0) as fh_in, aa_fasta_path.open('w') as fh_out:
            for record in SeqIO.parse(fh_in, 'fasta'):
                seq = str(record.seq)
                model_id = record.id
                cols = record.description.split(' ', 1)[1].split('~~~')
                min_id = 50
                min_query_cov = 80
                min_model_cov = 80
                dbxrefs = []
                gene = ''
                if(len(cols) == 3):
                    gene = cols[0]
                    product = cols[1]
                    if(cols[2] != ''):
                        dbxrefs = cols[2].split(',')
                elif(len(cols) == 6):
                    min_id = float(cols[0])
                    min_query_cov = float(cols[1])
                    min_model_cov = float(cols[2])
                    gene = cols[3]
                    product = cols[4]
                    if(cols[5] != ''):
                        dbxrefs = cols[5].split(',')
                else:
                    log.error(
                        'wrong description format in user aa! description=%s',
                        record.description
                    )        
                    raise ValueError(f'Wrong description format! description={record.description}')
                gene = gene.strip()
                product = product.strip()
                for dbxref in dbxrefs:
                    if(':' not in dbxref):
                        log.error(
                            'wrong dbxref format in user aa! id=%s, min-id=%f, min-query-cov=%f, min-model-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                            model_id, min_id, min_query_cov, min_model_cov, gene, product, dbxrefs
                        )        
                        raise ValueError(f'Wrong dbxref format! dbxref={dbxref}')
                if(product == ''):
                    raise ValueError('Protein product must not be empty!')

                fh_out.write(f">{model_id} UserProteins~~~{100}~~~{min_id}~~~{min_query_cov}~~~{min_model_cov}~~~{gene}~~~{product}~~~{','.join(dbxrefs)}\n")
                fh_out.write(f"{seq}\n")
                log.debug(
                    'imported user aa: id=%s, length=%i, min-id=%f, min-query-cov=%f, min-model-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                    model_id, len(seq), min_id, min_query_cov, min_model_cov, gene, product, dbxrefs
                )
    except Exception as e:
        log.error('provided user proteins file format not valid!', exc_info=True)
        sys.exit(f'ERROR: User proteins file format not valid!')

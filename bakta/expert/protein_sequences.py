import logging
import subprocess as sp
import sys

from pathlib import Path
from typing import Sequence

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf

from Bio import SeqIO
from xopen import xopen


log = logging.getLogger('EXPERT_AA_SEQ')


def search(cdss: Sequence[dict], cds_fasta_path: Path, expert_system: str, db_path: Path):
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
    cds_by_hexdigest = orf.get_orf_dictionary(cdss)
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


def write_user_protein_sequences(aa_fasta_path: Path):
    user_proteins = []
    try:
        with xopen(str(cfg.user_proteins), threads=0) as fh_in:
            for record in SeqIO.parse(fh_in, 'fasta'):
                user_proteins.append(parse_user_protein_sequences_fasta(record))
    except Exception as e:
        log.error('provided user proteins file Fasta format not valid!', exc_info=True)
        sys.exit(f'ERROR: User proteins file Fasta format not valid!')

    try:
        with xopen(str(cfg.user_proteins), threads=0) as fh_in:
            for record in SeqIO.parse(fh_in, 'genbank'):
                for feature in record.features:
                    if(feature.type.lower() == 'cds'  and  'pseudo' not in feature.qualifiers):
                        user_proteins.append(parse_user_protein_sequences_genbank(feature))
    except Exception as e:
        log.error('provided user proteins file GenBank format not valid!', exc_info=True)
        sys.exit(f'ERROR: User proteins file GenBank format not valid!')

    if(len(user_proteins) > 0):
        try:
            with aa_fasta_path.open('w') as fh_out:
                for user_protein in user_proteins:
                    (model_id, min_id, min_query_cov, min_model_cov, gene, product, dbxrefs, seq) = user_protein
                    fh_out.write(f">{model_id} UserProteins~~~{100}~~~{min_id}~~~{min_query_cov}~~~{min_model_cov}~~~{gene}~~~{product}~~~{','.join(dbxrefs)}\n{seq}\n")
                    log.debug(
                        'imported user aa: id=%s, length=%i, min-id=%f, min-query-cov=%f, min-model-cov=%f, gene=%s, product=%s, dbxrefs=%s',
                        model_id, len(seq), min_id, min_query_cov, min_model_cov, gene, product, dbxrefs
                    )
        except Exception as e:
            log.error('cannot write user protein file!', exc_info=True)
            sys.exit(f'ERROR: Cannot write user protein file!')
    else:
        log.error('no user proteins detected!', exc_info=True)
        sys.exit(f'ERROR: No user proteins detected in file!')


def parse_user_protein_sequences_fasta(record):
    model_id = record.id
    min_id = bc.MIN_PSC_IDENTITY * 100
    min_query_cov = bc.MIN_PSC_COVERAGE * 100
    min_model_cov = bc.MIN_PSC_COVERAGE * 100

    cols = record.description.split(' ', 1)[1].split('~~~')
    db_xrefs = []
    gene = ''
    if(len(cols) == 3):
        gene = cols[0]
        product = cols[1]
        if(cols[2] != ''):
            db_xrefs = cols[2].split(',')
    elif(len(cols) == 6):
        min_id = float(cols[0])
        min_query_cov = float(cols[1])
        min_model_cov = float(cols[2])
        gene = cols[3]
        product = cols[4]
        if(cols[5] != ''):
            db_xrefs = cols[5].split(',')
    else:
        log.error(
            'wrong description format in Fasta user protein file! description=%s',
            record.description
        )
        raise ValueError(f'Wrong description format in Fasta user protein file! description={record.description}')

    gene = gene.strip()
    product = product.strip()
    if(product == ''):
        log.error('missing product in Fasta user protein file!')
        raise ValueError('Missing product in Fasta user protein file!')

    for db_xref in db_xrefs:
        if(':' not in db_xref):
            log.error('wrong dbxref format in Fasta user protein file!')
            raise ValueError(f'Wrong dbxref format in Fasta user protein file! dbxref={db_xref}')
    seq = str(record.seq).upper()
    return (model_id, min_id, min_query_cov, min_model_cov, gene, product, db_xrefs, seq)


def parse_user_protein_sequences_genbank(feature: dict):
    min_id = bc.MIN_PSC_IDENTITY * 100
    min_query_cov = bc.MIN_PSC_COVERAGE * 100
    min_model_cov = bc.MIN_PSC_COVERAGE * 100

    q = feature.qualifiers
    model_ids = q.get('locus_tag', [])
    model_id = model_ids[0] if len(model_ids) > 0 else None
    if(model_id is None):
        log.error('missing locus_tag in GenBank user protein file!')
        raise ValueError('Missing locus_tag in GenBank user protein file!')

    genes = q.get('gene', [])
    gene = genes[0].strip() if len(genes) > 0 else ''
    
    products = q.get('product', [])
    product = products[0].strip() if len(products) > 0 else None
    if(product is None or product == ''):
        log.error('missing product in GenBank user protein file!')
        raise ValueError('Missing product in GenBank user protein file!')

    db_xrefs = q.get('db_xref', [])
    for db_xref in db_xrefs:
        if(':' not in db_xref):
            log.error('wrong dbxref format in Fasta user protein file!')
            raise ValueError(f'Wrong dbxref format in Fasta user protein file! dbxref={db_xref}')
    ecs = q.get('EC_number', [])
    for ec in ecs:
        db_xrefs.append(f'EC:{ec}')

    seqs = q.get('translation', [])
    seq = str(seqs[0]).upper() if len(seqs) > 0 else None
    if(seq is None):
        log.error('missing sequence translation in GenBank user protein file!')
        raise ValueError('Missing sequence translation in GenBank user protein file!')

    return (model_id, min_id, min_query_cov, min_model_cov, gene, product, db_xrefs, seq)

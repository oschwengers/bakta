
import logging
import subprocess as sp

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:r_rna')


def predict_r_rnas(data, contigs_path):
    """Search for ribosomal RNA sequences."""

    output_path = cfg.tmp_path.joinpath('rrna.tsv')

    cmd = [
        'cmsearch',
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path),
        str(cfg.db_path.joinpath('rRNA')),
        str(contigs_path)
    ]
    if(data['genome_size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(data['genome_size'] // 1000000))
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning('rRNAs failed! cmscan-error-code=%d', proc.returncode)
        log.debug(
            'rRNAs: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("cmsearch error! error code: %i" % proc.returncode)

    rrnas = []
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, accession, subject, subject_id, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = line.strip().split()
                
                if(strand == '-'):
                    (start, stop) = (stop, start)
                
                db_xrefs = ['GO:0005840', 'GO:0003735']
                if(subject_id == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += ['RFAM:RF00001', 'SO:0000652']
                elif(subject_id == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += ['RFAM:RF00177', 'SO:0001000']
                elif(subject_id == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += ['RFAM:RF02541', 'SO:0001001']

                rrna = {
                    'type': bc.FEATURE_R_RNA,
                    'gene': "%s_rrna" % rrna_tag,
                    'product': "%s ribosomal RNA" % rrna_tag,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'score': float(score),
                    'evalue': float(evalue),
                    'db_xrefs': db_xrefs
                }
                rrnas.append(rrna)
                log.debug(
                    'rRNA: contig=%s, gene=%s, start=%i, stop=%i, strand=%s',
                    rrna['contig'], rrna['gene'], rrna['start'], rrna['stop'], rrna['strand']
                )
    log.info('rRNAs: # %i', len(rrnas))
    return rrnas

import logging
import re
import subprocess as sp

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:r_rna')


def predict_r_rnas(data, contigs_path):
    """Search for ribosomal RNA sequences."""

    output_path = cfg.tmp_path.joinpath('rrna.tsv')
    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '-g',  # activate glocal mode
        '--nohmmonly',  # strictly use CM models
        '--rfam',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path)
    ]
    if(data['genome_size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(2 * data['genome_size'] // 1000000))
    cmd.append(str(cfg.db_path.joinpath('rRNA')))
    cmd.append(str(contigs_path))
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
        log.warning('rRNAs failed! cmscan-error-code=%d', proc.returncode)
        raise Exception("cmscan error! error code: %i" % proc.returncode)

    rrnas = []
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (subject, accession, contig_id, contig_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = re.split('\s+', line.strip(), maxsplit=17)
                
                if(strand == '-'):
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                length = stop - start + 1
                partial = trunc != 'no'

                db_xrefs = ['GO:0005840', 'GO:0003735']
                if(accession == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += ['RFAM:RF00001', 'SO:0000652']
                    consensus_length = 119
                elif(accession == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += ['RFAM:RF00177', 'SO:0001000']
                    consensus_length = 1533
                elif(accession == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += ['RFAM:RF02541', 'SO:0001001']
                    consensus_length = 2925
                
                coverage = length / consensus_length
                if( coverage < 0.8):
                    partial = True
                
                if(coverage < 0.3):
                    log.debug(
                        'discard low coverage: contig=%s, rRNA=%s, start=%i, stop=%i, strand=%s, length=%i, coverage=%0.3f',
                        contig_id, rrna_tag, start, stop, strand, length, coverage
                    )
                else:
                    rrna = {
                        'type': bc.FEATURE_R_RNA,
                        'gene': "%s_rrna" % rrna_tag,
                        'product': "(partial) %s ribosomal RNA" % rrna_tag if partial else "%s ribosomal RNA" % rrna_tag,
                        'contig': contig_id,
                        'start': start,
                        'stop': stop,
                        'strand': strand,
                        'partial': partial,
                        'coverage': coverage,
                        'score': float(score),
                        'evalue': float(evalue),
                        'db_xrefs': db_xrefs
                    }
                    if('5' in trunc):
                        rrna['trunc_5'] = True
                    if('3' in trunc):
                        rrna['trunc_3'] = True
                    rrnas.append(rrna)
                    log.info(
                        'contig=%s, gene=%s, start=%i, stop=%i, strand=%s, partial=%s, length=%i, coverage=%0.3f',
                        rrna['contig'], rrna['gene'], rrna['start'], rrna['stop'], rrna['strand'], partial, length, coverage
                    )

    log.info('# %i', len(rrnas))
    return rrnas
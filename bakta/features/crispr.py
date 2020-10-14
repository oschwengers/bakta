
import logging
import subprocess as sp

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:crispr')


def predict_crispr(data, contigs_path):
    """Predict CRISPR arrays with PILER-CR."""
    # SO:0001459 <- Sequence Ontology

    output_path = cfg.tmp_path.joinpath('crispr.txt')
    cmd = [
        'pilercr',
        '-in', str(contigs_path),
        '-out', str(output_path),
        '-noinfo',  # omit help in output
        '-quiet'  # silent mode
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
        log.warning('CRISPRs failed! pilercr-error-code=%d', proc.returncode)
        raise Exception("PILER-CR error! error code: %i" % proc.returncode)

    # parse orfs
    crispr_arrays = []
    with output_path.open() as fh:
        contig_id = None
        skip_lines = True
        for line in fh:
            line = line.strip()
            if(line == 'SUMMARY BY POSITION'):
                skip_lines = False
            elif(skip_lines is False):
                if(len(line) > 0):
                    if(line[0] == '>'):
                        contig_id = line[1:]
                    elif(line[0] != 'A' and line[0] != '='):
                        cols = line.split()
                        if(len(cols) == 8):
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, repeat_consensus) = cols
                        else:
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, distance, repeat_consensus) = cols
                        crispr = {
                            'type': bc.FEATURE_CRISPR,
                            'contig': contig_id,
                            'start': int(position),
                            'stop': int(position) + int(length) - 1,
                            'strand': '+',
                            'product': "CRISPR array with %s repeats of length %s, consensus sequence %s and spacer length %s" % (copies, repeat_length, repeat_consensus, spacer_length),
                            'spacer_length': int(spacer_length),
                            'repeat_length': int(repeat_length),
                            'repeats': int(copies),
                            'repeat_consensus': repeat_consensus
                        }
                        crispr_arrays.append(crispr)
                        log.debug(
                            'contig=%s, start=%i, stop=%i, spacer-length=%i, repeat-length=%i, # repeats=%i, repeat-consensus=%s',
                            crispr['contig'], crispr['start'], crispr['stop'], crispr['spacer_length'], crispr['repeat_length'], crispr['repeats'], crispr['repeat_consensus']
                        )
    log.info('# %i', len(crispr_arrays))
    return crispr_arrays
    


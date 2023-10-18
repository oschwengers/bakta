import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('CRISPR')


def predict_crispr(genome: dict, contigs_path: Path):
    """Predict CRISPR arrays with PILER-CR."""

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
        raise Exception(f'PILER-CR error! error code: {proc.returncode}')

    # parse crispr arrays
    crispr_arrays = {}
    contigs = {c['id']: c for c in genome['contigs']}
    with output_path.open() as fh:
        output_section = None
        contig_id = None
        array_id = None
        skip_lines = True
        crispr_array = None
        gap_count = 0
        for line in fh:
            line = line.strip()
            if(line == ''):
                continue
            if(line == 'DETAIL REPORT'):
                output_section = 'DETAIL'
                skip_lines = False
            elif(line == 'SUMMARY BY POSITION'):
                output_section = 'POSITION'
                skip_lines = False
            elif(line == 'SUMMARY BY SIMILARITY'):
                output_section = 'SIMILARITY'
                skip_lines = False
            elif(skip_lines is False):
                if(output_section == 'DETAIL'):
                    if(line[0:5] == 'Array'):
                        gap_count = 0
                        array_id = line.split()[1]
                        crispr_array = OrderedDict()
                        crispr_array['type'] = bc.FEATURE_CRISPR
                        crispr_array['strand'] = bc.STRAND_UNKNOWN
                        crispr_array['repeats'] = []
                        crispr_array['spacers'] = []
                        crispr_arrays[array_id] = crispr_array
                    elif(line[0] == '>'):
                        contig_id = line[1:]
                        crispr_array['contig'] = contig_id
                    elif(line[0] != '='):
                        cols = line.split()
                        if(len(cols) == 7  and  cols[0] != 'Pos'):
                            (position, repeat_length, id, spacer_length, left_flank, repeat_seq, spacer_seq) = cols
                            position, repeat_length, spacer_length = int(position), int(repeat_length), int(spacer_length)
                            spacer_seq = spacer_seq.upper()
                            crispr_repeat = OrderedDict()
                            crispr_repeat['strand'] = bc.STRAND_UNKNOWN
                            crispr_repeat['start'] = position - gap_count
                            crispr_repeat['stop'] = position + repeat_length - 1 - gap_count
                            crispr_array['repeats'].append(crispr_repeat)
                            gap_count += repeat_seq.count('-')  # correct wrong PILER-CR detail positions by gaps
                            crispr_spacer = OrderedDict()
                            crispr_spacer['strand'] = bc.STRAND_UNKNOWN
                            crispr_spacer['start'] = position + repeat_length  - gap_count
                            crispr_spacer['stop'] = position + repeat_length + spacer_length - 1 - gap_count
                            crispr_spacer['sequence'] = spacer_seq
                            crispr_array['spacers'].append(crispr_spacer)
                            spacer_genome_seq = bu.extract_feature_sequence(crispr_spacer, contigs[contig_id])
                            assert spacer_seq == spacer_genome_seq  # assure PILER-CR spacer sequence equal extraction from genome
                        elif(len(cols) == 6  and  cols[0] != 'Pos'):  # last line in array without spacer
                            (position, repeat_length, id, left_flank, repeat_seq, spacer_seq) = cols
                            position, repeat_length, spacer_length = int(position), int(repeat_length), int(spacer_length)
                            crispr_repeat = OrderedDict()
                            crispr_repeat['strand'] = bc.STRAND_UNKNOWN
                            crispr_repeat['start'] = position - gap_count
                            crispr_repeat['stop'] = position + repeat_length - 1 - gap_count
                            crispr_array['repeats'].append(crispr_repeat)
                elif(output_section == 'POSITION'):
                    if(line[0] == '>'):
                        contig_id = line[1:]
                    elif(line[0] != 'A' and line[0] != '='):
                        cols = line.split()
                        if(len(cols) == 8):
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, repeat_consensus) = cols
                        else:
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, distance, repeat_consensus) = cols
                        crispr_array = crispr_arrays[array_id]
                        crispr_array['start'] = int(position)
                        crispr_array['stop'] = int(position) + int(length) - 1
                        crispr_array['product'] = f'CRISPR array with {copies} repeats of length {repeat_length}, consensus sequence {repeat_consensus} and spacer length {spacer_length}'
                        crispr_array['spacer_length'] = int(spacer_length)
                        crispr_array['repeat_length'] = int(repeat_length)
                        assert len(crispr_array['repeats']) == int(copies), print(f"len(reps)={len(crispr_array['repeats'])}, int(copies)={int(copies)}")
                        crispr_array['repeat_consensus'] = repeat_consensus
                        crispr_array['db_xrefs'] = [so.SO_CRISPR.id]

                        nt = bu.extract_feature_sequence(crispr_array, contigs[contig_id])  # extract nt sequences
                        crispr_array['nt'] = nt
                        log.info(
                            'contig=%s, start=%i, stop=%i, spacer-length=%i, repeat-length=%i, # repeats=%i, repeat-consensus=%s, nt=[%s..%s]',
                            crispr_array['contig'], crispr_array['start'], crispr_array['stop'], crispr_array['spacer_length'], crispr_array['repeat_length'], len(crispr_array['repeats']), crispr_array['repeat_consensus'], nt[:10], nt[-10:]
                        )
    crispr_arrays = crispr_arrays.values()                        
    log.info('predicted=%i', len(crispr_arrays))
    return crispr_arrays

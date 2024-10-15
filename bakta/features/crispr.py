import logging
import re
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


RE_CRISPR = re.compile(r'(\d{1,8})\s+(\d{2})\s+(\d{1,3}\.\d)\s+(?:(\d{1,2})\s+)?([ATGCN]+)?\s+([ATGCN\.-]+)\s*(?:([ATGCN]+))?')


log = logging.getLogger('CRISPR')


def predict_crispr(data: dict, sequences_path: Path):
    """Predict CRISPR arrays with PILER-CR."""

    output_path = cfg.tmp_path.joinpath('crispr.txt')
    cmd = [
        'pilercr',
        '-in', str(sequences_path),
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
    sequences = {seq['id']: seq for seq in data['sequences']}
    with output_path.open() as fh:
        output_section = None
        sequence_id = None
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
                        sequence_id = line[1:]
                        crispr_array['sequence'] = sequence_id
                    elif(line[0] != '='):
                        m = RE_CRISPR.fullmatch(line)
                        if(m is not None):
                            position = int(m.group(1))
                            repeat_length = int(m.group(2))
                            repeat_seq = m.group(6)
                            spacer_seq = m.group(7)
                            crispr_repeat = OrderedDict()
                            crispr_repeat['strand'] = bc.STRAND_UNKNOWN
                            crispr_repeat['start'] = position - gap_count
                            crispr_repeat['stop'] = position + repeat_length - 1 - gap_count
                            crispr_array['repeats'].append(crispr_repeat)
                            log.debug('repeat: array-id=%s, start=%i, stop=%i', array_id, crispr_repeat['start'], crispr_repeat['stop'])
                            gap_count += repeat_seq.count('-')  # correct wrong PILER-CR detail positions by gaps
                            if(spacer_seq is not None):
                                spacer_seq = spacer_seq.upper()
                                spacer_length = len(spacer_seq)
                                crispr_spacer = OrderedDict()
                                crispr_spacer['strand'] = bc.STRAND_UNKNOWN
                                crispr_spacer['start'] = position + repeat_length - gap_count
                                crispr_spacer['stop'] = position + repeat_length + spacer_length - 1 - gap_count
                                crispr_spacer['sequence'] = spacer_seq
                                crispr_array['spacers'].append(crispr_spacer)
                                spacer_genome_seq = bu.extract_feature_sequence(crispr_spacer, sequences[sequence_id])
                                log.debug('spacer: array-id=%s, start=%i, stop=%i, genome-seq=%s, spacer-seq=%s', array_id, crispr_spacer['start'], crispr_spacer['stop'], spacer_genome_seq, spacer_seq)
                                assert spacer_seq == spacer_genome_seq  # assure PILER-CR provided sequence equals sequence extracted from genome
                elif(output_section == 'POSITION'):
                    if(line[0] == '>'):
                        sequence_id = line[1:]
                    elif(line[0] != 'A' and line[0] != '='):
                        cols = line.split()
                        if(len(cols) == 8):
                            (array_id, sequence, position, length, copies, repeat_length, spacer_length, repeat_consensus) = cols
                        else:
                            (array_id, sequence, position, length, copies, repeat_length, spacer_length, distance, repeat_consensus) = cols
                        crispr_array = crispr_arrays[array_id]
                        positions = [seq['start'] for seq in crispr_array['spacers']] + [seq['stop'] for seq in crispr_array['spacers']] + [seq['start'] for seq in crispr_array['repeats']] + [seq['stop'] for seq in crispr_array['repeats']]
                        crispr_array['start'] = min(positions)
                        crispr_array['stop'] = max(positions)
                        crispr_array['product'] = f'CRISPR array with {copies} repeats of length {repeat_length}, consensus sequence {repeat_consensus} and spacer length {spacer_length}'
                        crispr_array['spacer_length'] = int(spacer_length)
                        crispr_array['repeat_length'] = int(repeat_length)
                        assert (len(crispr_array['repeats']) - int(copies)) <= 1, print(f"len(reps)={len(crispr_array['repeats'])}, int(copies)={int(copies)}")
                        crispr_array['repeat_consensus'] = repeat_consensus
                        crispr_array['db_xrefs'] = [so.SO_CRISPR.id]

                        nt = bu.extract_feature_sequence(crispr_array, sequences[sequence_id])  # extract nt sequences
                        crispr_array['nt'] = nt
                        log.info(
                            'seq=%s, start=%i, stop=%i, spacer-length=%i, repeat-length=%i, # repeats=%i, repeat-consensus=%s, nt=[%s..%s]',
                            crispr_array['sequence'], crispr_array['start'], crispr_array['stop'], crispr_array['spacer_length'], crispr_array['repeat_length'], len(crispr_array['repeats']), crispr_array['repeat_consensus'], nt[:10], nt[-10:]
                        )
    crispr_arrays = crispr_arrays.values()                        
    log.info('predicted=%i', len(crispr_arrays))
    return crispr_arrays

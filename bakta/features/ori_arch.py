import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path
from typing import Sequence

import bakta.config as cfg
import bakta.constants as bc
import bakta.io.fasta as fasta
import bakta.utils as bu


log = logging.getLogger('ORI_ARCH')


def predict_orics(data: dict, sequences_path: Path) -> Sequence[dict]:
    """Predict archaeal replication origins with Ori-Finder-Arch.

    Ori-Finder-Arch emits GFF3 records with type ``rep_origin``.  Its own RIP
    CDS records are deliberately not imported: Bakta's CDS caller remains the
    sole source of CDS features.
    """
    sequences = {seq['id']: seq for seq in data['sequences']}
    orics = []
    # Ori-Finder-Arch accepts multi-record FASTA, but topology and assembly
    # level are global options. Invoke it per replicon so mixed assemblies are
    # never silently downgraded to a single linear/draft interpretation.
    for index, sequence in enumerate(data['sequences'], start=1):
        input_path = cfg.tmp_path.joinpath(f'ori-arch.{index}.fna')
        output_path = cfg.tmp_path.joinpath(f'ori-arch.{index}.gff3')
        fasta.export_sequences([sequence], input_path)
        topology = sequence['topology']
        level = 'complete' if sequence.get('complete', data['genome']['complete']) else 'draft'
        _run_ori_finder(input_path, output_path, topology, level)
        orics.extend(_parse_orics(output_path, sequences))
    log.info('predicted=%i', len(orics))
    return orics


def _run_ori_finder(input_path: Path, output_path: Path, topology: str, level: str):
    cmd = ['OriFinderArch', '-i', str(input_path), '-o', str(output_path), '-t', topology, '-l', level]
    log.debug('cmd=%s', cmd)
    proc = sp.run(cmd, cwd=str(cfg.tmp_path), env=cfg.env, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    if(proc.returncode != 0):
        log.debug('stdout=%s, stderr=%s', proc.stdout, proc.stderr)
        raise Exception(f'Ori-Finder-Arch error! exit_code={proc.returncode}; stderr={proc.stderr.strip()}')


def _parse_orics(output_path: Path, sequences: dict) -> Sequence[dict]:
    if(not output_path.is_file()):
        raise Exception('Ori-Finder-Arch completed without producing a GFF3 result')
    orics = []
    with output_path.open() as fh:
        for line in fh:
            if(line.startswith('#')):
                continue
            fields = line.rstrip('\n').split('\t')
            if(len(fields) != 9):
                log.warning('skip malformed Ori-Finder-Arch GFF3 line: %s', line.rstrip())
                continue
            sequence_id, source, feature_type, start, stop, score, strand, phase, attributes = fields
            if(feature_type != 'rep_origin' or sequence_id not in sequences):
                continue
            ori = OrderedDict()
            ori['type'] = bc.FEATURE_ORIC
            ori['sequence'] = sequence_id
            ori['start'] = int(start)
            ori['stop'] = int(stop)
            ori['strand'] = bc.STRAND_UNKNOWN
            ori['product'] = 'origin of replication'
            ori['inference'] = 'ab initio prediction:Ori-Finder-Arch:1.0.0'
            ori['attributes'] = _parse_attributes(attributes)
            if(score != '.'):
                ori['score'] = float(score)
            ori['nt'] = bu.extract_feature_sequence(ori, sequences[sequence_id])
            orics.append(ori)
            log.info('oriC: seq=%s, start=%i, stop=%i', sequence_id, ori['start'], ori['stop'])
    return orics


def _parse_attributes(raw: str) -> dict:
    return {key: value.strip("'") for key, value in (field.split('=', 1) for field in raw.split(';') if '=' in field)}

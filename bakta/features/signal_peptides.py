import logging
import sys
import subprocess as sp

from typing import Sequence, Tuple
from pathlib import Path

import bakta.constants as bc
import bakta.config as cfg

log = logging.getLogger('SIGNAL_PEPTIDES')

def search(orfs: Sequence[dict], orf_aa_path: Path):
    """Search for signal peptides with DeepSig."""

    deepsig_output_path = cfg.tmp_path.joinpath('deepsig.gff3')
    gram = 'gramn' if cfg.gram == '-' else 'gramp'
    cmd = [
        'deepsig',
        '--fasta', str(orf_aa_path),
        '--organism', gram,
        '--outf', str(deepsig_output_path),
        '--outfmt', 'gff3',
        '--threads', str(cfg.threads)
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
        log.warning('signal peptides failed! deep-sig-error-code=%d', proc.returncode)
        raise Exception(f'deepsig error! error code: {proc.returncode}')

    sequences_by_hexdigest = {f"{orf['aa_hexdigest']}": orf for orf in orfs}

    sig_peps = []
    with deepsig_output_path.open() as fh:
        for line in fh:
            (identifier, tool, feature_type, start_aa, stop_aa, score, placeholder, placeholder2, description) = line.split('\t')
            if feature_type == 'Signal peptide':
                start_aa = int(start_aa)
                stop_aa = int(stop_aa)
                score = float(score)
                aa_identifier = identifier.split('-')[0]
                if (aa_identifier in sequences_by_hexdigest):
                    orf = sequences_by_hexdigest[aa_identifier]
                    start_nt, stop_nt = orf_nt_start_stop(orf, start_aa, stop_aa)

                    sig_pep = {
                        'type': bc.FEATURE_SIGNAL_PEPTIDE,
                        'start': start_nt,
                        'stop': stop_nt,
                        'score': score
                    }
                    if (bc.FEATURE_SIGNAL_PEPTIDE not in orf):
                        orf[bc.FEATURE_SIGNAL_PEPTIDE] = {}
                    orf[bc.FEATURE_SIGNAL_PEPTIDE] = sig_pep
                    log.debug(
                        'hit: contig=%s, nt-start=%i, nt-stop=%i, aa-start=%i, aa-stop=%i, score=%0.2f',
                        orf['contig'], start_nt, stop_nt, start_aa, stop_aa, score
                    )
                    sig_peps.append(sig_pep)
                else:
                    log.error('signal peptide: unknown aa_identifier=%s', aa_identifier)
                    sys.exit('ERROR: signal peptide found for unknown aa_identifier=%s!', aa_identifier)
    return sig_peps


def orf_nt_start_stop(orf: dict, start_aa: int, stop_aa: int) -> Tuple[int, int]:
    """Determin signal peptide nucleotide position.
        orf: dictonary of ORF, either CDS or sORF
        start_aa: sig pep start position within aa seq
        stop_aa: sig pep stop position within aa seq
    """
    if(orf['strand'] == bc.STRAND_FORWARD):
        start_nt = orf['start'] + ((start_aa-1) * 3)
        stop_nt = orf['start'] + ((stop_aa-1) * 3) + 2
    else:
        start_nt = orf['stop'] - ((stop_aa-1) * 3 + 2)
        stop_nt = orf['stop'] - ((start_aa-1) * 3)
    return start_nt, stop_nt

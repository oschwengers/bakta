import logging
import sys
import subprocess as sp
from collections import OrderedDict

import bakta
import bakta.constants as bc
import bakta.config as cfg

log = logging.getLogger('SIGNAL_PEPTIDES')

def execute_deepsig(orfs, orf_fasta_path):
    """Search for Signal Peptides using DeepSig tool."""
    deepsig_output_path = cfg.tmp_path.joinpath('deepsig.gff3')
    gram = 'gramn' if cfg.gram == '-' else 'gramp'
    cmd = [
        'deepsig.py',
        '-f', str(orf_fasta_path),
        '-f', str(seqs_fasta_path),
        '-k', gram,
        '-o', str(deepsig_output_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        log.error('signal peptide: deepsig failed to run through: command=%s, stdout=%s, stderr=%s, error-code=%s', cmd, proc.stdout, proc.stderr, proc.returncode)
        sys.exit('ERROR: deepsig failed to run through! error code: %s, error output: %s', proc.returncode, proc.stderr)

    sequences_by_hexdigest = {f"{orf['aa_hexdigest']}": orf for orf in orfs}
    sig_pep = OrderedDict()
    sig_peps = []
    with deepsig_output_path.open() as fh:
        for line in fh:
            feature_type = line.split("\t")[2]
            if feature_type == "Signal peptide":
                identifier, tool, feature_type, start, stop, feature_annotation_score, placeholder, placeholder2, description = line.split("\t")
                aa_identifier = identifier.split("-")[0]
                orf = sequences_by_hexdigest[aa_identifier]

                start_nucleotides, stop_nucleotides = start_stop_orf(orf, start, stop)

                sig_pep['type'] = bc.FEATURE_SIGNAL_PEPTIDE
                sig_pep['start'] = start_nucleotides
                sig_pep['stop'] = stop_nucleotides
                sig_pep['feature_annotation_score'] = feature_annotation_score
                if (bc.FEATURE_SIGNAL_PEPTIDE not in orf):
                    orf[bc.FEATURE_SIGNAL_PEPTIDE]={}
                orf[bc.FEATURE_SIGNAL_PEPTIDE]=sig_pep
                log.debug(
                    'hit: identifier=%s, start=%s, stop=%s, feature annotation score=%s',
                    aa_identifier, start_nucleotides, stop_nucleotides, feature_annotation_score
                )
                sig_peps.append(sig_pep)
    return sig_peps

def start_stop_orf(orf, start, stop):
    """Method for determining correct position of the signal peptide on nucleotide sequence."""
    start = int(start)
    stop = int(stop)
    if orf['strand'] == '-':
        start_nucleotides = orf['stop']-((stop-1)*3+2)
        stop_nucleotides = orf['stop']-((start-1)*3)
    else:
        start_nucleotides = orf['start']+((start-1)*3)
        stop_nucleotides = orf['start']+((stop-1)*3)+2
    return start_nucleotides,stop_nucleotides

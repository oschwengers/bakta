import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path
from typing import Dict, Sequence

import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('ORF')


def detect_spurious(orfs: Sequence[dict], orf_aa_path: Path):
    """Detect spurious ORFs with AntiFam"""
    output_path = cfg.tmp_path.joinpath('cds.spurious.hmm.tsv')
    cmd = [
        'hmmsearch',
        '--noali',
        '--cut_ga',  # use gathering cutoff
        '--tblout', str(output_path),
        '--cpu', str(cfg.threads if cfg.threads <= 4 else 4),
        str(cfg.db_path.joinpath('antifam')),
        str(orf_aa_path)
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
        log.warning('spurious ORF detection failed! hmmsearch-error-code=%d', proc.returncode)
        raise Exception(f'hmmsearch error! error code: {proc.returncode}')

    discarded_orfs = []
    orf_by_aa_digest = get_orf_dictionary(orfs)
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (aa_identifier, _, subject_name, subject_id, evalue, bitscore, _) = line.strip().split(maxsplit=6)
                orf = orf_by_aa_digest[aa_identifier]
                evalue = float(evalue)
                bitscore = float(bitscore)
                if(evalue > 1E-5):
                    log.debug(
                        'discard low spurious E value: contig=%s, start=%i, stop=%i, strand=%s, subject=%s, evalue=%1.1e, bitscore=%f',
                        orf['contig'], orf['start'], orf['stop'], orf['strand'], subject_name, evalue, bitscore
                    )
                else:
                    discard = OrderedDict()
                    discard['type'] = bc.DISCARD_TYPE_SPURIOUS
                    discard['description'] = f'(partial) homology to spurious sequence HMM (AntiFam:{subject_id})'
                    discard['score'] = bitscore
                    discard['evalue'] = evalue

                    orf['discarded'] = discard
                    discarded_orfs.append(orf)
                    log.info(
                        'discard spurious: contig=%s, start=%i, stop=%i, strand=%s, homology=%s, evalue=%1.1e, bitscore=%f',
                        orf['contig'], orf['start'], orf['stop'], orf['strand'], subject_name, evalue, bitscore
                    )
    log.info('discarded=%i', len(discarded_orfs))
    return discarded_orfs


def get_orf_key(orf: dict) -> str:
    """Generate a standardized and unique ORF-like feature key for internal store/analyze/parse/retrieval cycles."""
    return f"{orf['aa_hexdigest']}-{orf['contig']}-{orf['start']}-{orf['stop']}-{orf['strand']}"


def get_orf_dictionary(orfs: Sequence[dict]) -> Dict[str, dict]:
    """create a standardized ORF-like feature dict for internal store/analyze/parse/retrieval cycles."""
    return {get_orf_key(orf): orf for orf in orfs}


def write_internal_faa(features: Sequence[dict], faa_path: Path):
    """Write aa sequences to internal temporary Fasta file."""
    log.info('write internal aa seqs: # seqs=%i, path=%s', len(features), faa_path)
    with faa_path.open(mode='wt') as fh:
        for orf in features:
            fh.write(f">{get_orf_key(orf)}\n{orf['aa']}\n")
    
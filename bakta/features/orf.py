import logging

from collections import OrderedDict
from pathlib import Path
from typing import Dict, Sequence

import pyhmmer

import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('ORF')


def detect_spurious(orfs: Sequence[dict]):
    """Detect spurious ORFs with AntiFam"""
    discarded_orfs = []
    orf_by_aa_digest = get_orf_dictionary(orfs)
    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    proteins: list[pyhmmer.easel.DigitalSequence] = [pyhmmer.easel.TextSequence(sequence=orf['aa'], name=bytes(get_orf_key(orf), 'UTF-8')).digitize(alphabet) for orf in orfs]
    with pyhmmer.plan7.HMMFile(cfg.db_path.joinpath('antifam')) as hmm:
        for top_hits in pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs='gathering', cpus=cfg.threads):
            for hit in top_hits:
                orf = orf_by_aa_digest[hit.name.decode()]
                if hit.evalue > bc.MIN_HMM_EVALUE:
                    log.debug(
                        'discard low spurious E value: seq=%s, start=%i, stop=%i, strand=%s, subject=%s, evalue=%1.1e, bitscore=%f',
                        orf['sequence'], orf['start'], orf['stop'], orf['strand'], hit.best_domain.alignment.hmm_name.decode(), hit.evalue, hit.score
                    )
                else:
                    discard = OrderedDict()
                    discard['type'] = bc.DISCARD_TYPE_SPURIOUS
                    discard['description'] = f'(partial) homology to spurious sequence HMM (AntiFam:{hit.best_domain.alignment.hmm_accession.decode()})'
                    discard['score'] = hit.score
                    discard['evalue'] = hit.evalue

                    orf['discarded'] = discard
                    discarded_orfs.append(orf)
                    log.info(
                        'discard spurious: seq=%s, start=%i, stop=%i, strand=%s, homology=%s, evalue=%1.1e, bitscore=%f',
                        orf['sequence'], orf['start'], orf['stop'], orf['strand'], hit.best_domain.alignment.hmm_name.decode(), hit.evalue, hit.score
                    )
    log.info('discarded=%i', len(discarded_orfs))
    return discarded_orfs


def get_orf_key(orf: dict) -> str:
    """Generate a standardized and unique ORF-like feature key for internal store/analyze/parse/retrieval cycles."""
    return f"{orf['aa_hexdigest']}-{orf['sequence']}-{orf['start']}-{orf['stop']}-{orf['strand']}-{orf.get('source', 'internal')}"


def get_orf_dictionary(orfs: Sequence[dict]) -> Dict[str, dict]:
    """create a standardized ORF-like feature dict for internal store/analyze/parse/retrieval cycles."""
    return {get_orf_key(orf): orf for orf in orfs}


def write_internal_faa(features: Sequence[dict], faa_path: Path):
    """Write aa sequences to internal temporary Fasta file."""
    log.info('write internal aa seqs: # seqs=%i, path=%s', len(features), faa_path)
    with faa_path.open(mode='wt') as fh:
        for orf in features:
            fh.write(f">{get_orf_key(orf)}\n{orf['aa']}\n")
    
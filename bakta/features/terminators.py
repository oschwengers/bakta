import logging
import re
import shutil
import subprocess as sp

from pathlib import Path
from typing import Sequence

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu
import bakta.so as so
import bakta.features.orf as orf


HIT_SCORE_THRESHOLD = 70  # minimum confidence score for a terminator prediction


log = logging.getLogger('TERMINATOR')


TRANSTERMHP_BAG_FORMAT = re.compile(r'([A-Z0-9_\-+]+)\s+(\d{1,7})\s+..\s+(\d{1,7})\s([+-])\s+-[0-9.]+\s+-[0-9.]+\s+[ATCG-]+\s+[ATCG-]+\s+[ATCG-]+\s+[ATCG-]+\s+[ATCG-]+\s+(\d+)\s\d+', re.IGNORECASE)


def export_gene_positions(features: Sequence[dict], output_path: Path):
    """Export gene positions to a tab-delimited file."""
    with output_path.open('w') as f:
        for feature in features:
            start, stop = (feature['start'], feature['stop']) if feature['strand'] == bc.STRAND_FORWARD else (feature['stop'], feature['start'])
            f.write(f"{orf.get_orf_key(feature)}\t{start}\t{stop}\t{feature['sequence']}\n")


def predict_terminators(data: dict, sequences_path: Path) -> Sequence[dict]:
    """Predict rho-independent terminators using TransTermHP."""
    log.info('predicting rho-independent terminators')
    cds_sorf_features = [feat for feat in data['features'] if feat['type'] in [bc.FEATURE_CDS, bc.FEATURE_SORF] and 'discarded' not in feat]

    cds_sorf_dict = orf.get_orf_dictionary(cds_sorf_features)  # build CDS/sORF dictionary for TransTermHP

    coordinates_path = cfg.tmp_path.joinpath('genes.coords')
    export_gene_positions(cds_sorf_features, coordinates_path)
    expterm_path = Path(shutil.which('transterm')).resolve().parent.parent.joinpath('data', 'expterm.dat')
    output_bag_path = cfg.tmp_path.joinpath('bag.tsv')
    cmd = [
        'transterm',
        f'--min-conf={HIT_SCORE_THRESHOLD}',  # minimum confidence score for a terminator prediction to be reported
        '-p', str(expterm_path),  # experm.dat path
        '--bag-output', str(output_bag_path),  # output Best After Gene (BAG) terminators
        str(sequences_path),
        str(coordinates_path)
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
        log.warning('terminators failed! transtermhp-error-code=%d', proc.returncode)
        raise Exception(f'transtermhp error! error code: {proc.returncode}')

    sequences_per_id = {seq['id']: seq for seq in data['sequences']}
    with output_bag_path.open() as f:
        terminators = []
        for line in f:
            if 'NONE' not in line:
                match = TRANSTERMHP_BAG_FORMAT.fullmatch(line.strip())
                if(match):
                    orf_id = match.group(1)
                    seq_id = cds_sorf_dict[orf_id]['sequence']
                    start = int(match.group(2))
                    stop = int(match.group(3))
                    strand = bc.STRAND_FORWARD if match.group(4) == '+' else bc.STRAND_REVERSE
                    score = int(match.group(5))
                    if (score >= HIT_SCORE_THRESHOLD):
                        terminator = {
                            'type': bc.FEATURE_TERMINATOR,
                            'sequence': seq_id,
                            'start': start,
                            'stop': stop,
                            'strand': strand,
                            'product': 'Rho-independent terminator',
                            'score': score,
                            'class': so.SO_CIS_REG_TERMINATOR,
                            'db_xrefs': [so.SO_CIS_REG_TERMINATOR.id]
                        }
                        nt = bu.extract_feature_sequence(terminator, sequences_per_id[seq_id])  # extract nt sequences
                        terminator['nt'] = nt

                        terminators.append(terminator)
                        log.info(
                            'seq=%s, start=%i, stop=%i, strand=%s, length=%i, score=%d, nt=[%s..%s]',
                            terminator['sequence'], terminator['start'], terminator['stop'], terminator['strand'], len(nt), score, nt[:10], nt[-10:]
                        )
                    else:
                        log.debug(
                            'skipped terminator with low score: seq=%s, start=%i, stop=%i, strand=%s, length=%i, score=%d',
                            seq_id, start, stop, strand, stop - start + 1, score
                        )
                else:
                    log.warning('could not parse line: %s', line.strip())
    log.info('predicted=%i', len(terminators))
    return terminators

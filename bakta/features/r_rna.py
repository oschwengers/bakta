import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_COVERAGE = 0.3
HIT_COVERAGE_TRUNCATED = 0.8


log = logging.getLogger('R_RNA')


def predict_r_rnas(data: dict, sequences_path: Path):
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
    if(data['stats']['size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(2 * data['stats']['size'] // 1000000))
    cmd.append(str(cfg.db_path.joinpath('rRNA')))
    cmd.append(str(sequences_path))
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
        raise Exception(f'cmscan error! error code: {proc.returncode}')

    rrnas = []
    sequences = {seq['id']: seq for seq in data['sequences']}
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (
                    subject, accession, sequence_id, sequence_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description
                ) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

                if(strand == '-'):
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                score = float(score)
                length = stop - start + 1
                if(trunc == "5'"):
                    truncated = bc.FEATURE_END_5_PRIME
                elif(trunc == "3'"):
                    truncated = bc.FEATURE_END_3_PRIME
                else:
                    truncated = None

                db_xrefs = [f'{bc.DB_XREF_GO}:0005840', f'{bc.DB_XREF_GO}:0003735']
                if(accession == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00001', f'{bc.DB_XREF_KOFAM}:K01985', so.SO_RRNA_5S.id]
                    consensus_length = 119
                elif(accession == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00177', f'{bc.DB_XREF_KOFAM}:K01977', so.SO_RRNA_16S.id]
                    consensus_length = 1533
                elif(accession == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF02541', f'{bc.DB_XREF_KOFAM}:K01980', so.SO_RRNA_23S.id]
                    consensus_length = 2925
                else:
                    log.warning(
                        'unknown rRNA detected! accession=%s, seq=%s, start=%i, stop=%i, strand=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        accession, sequence_id, start, stop, strand, length, truncated, score, evalue
                    )
                    continue

                coverage = length / consensus_length
                if(coverage < HIT_COVERAGE_TRUNCATED):
                    truncated = bc.FEATURE_END_UNKNOWN

                if(coverage < HIT_COVERAGE):
                    log.debug(
                        'discard low coverage: seq=%s, rRNA=%s, start=%i, stop=%i, strand=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e',
                        sequence_id, rrna_tag, start, stop, strand, length, coverage, truncated, score, evalue
                    )
                else:
                    rrna = OrderedDict()
                    rrna['type'] = bc.FEATURE_R_RNA
                    rrna['sequence'] = sequence_id
                    rrna['start'] = start
                    rrna['stop'] = stop
                    rrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    if(accession == 'RF00001'):
                        rrna['gene'] = 'rrf'
                    elif(accession == 'RF00177'):
                        rrna['gene'] = 'rrs'
                    elif(accession == 'RF02541'):
                        rrna['gene'] = 'rrl'
                    rrna['product'] = f'{rrna_tag} ribosomal RNA'

                    if(truncated):
                        rrna['truncated'] = truncated

                    rrna['coverage'] = coverage
                    rrna['score'] = score
                    rrna['evalue'] = evalue
                    rrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(rrna, sequences[sequence_id])  # extract nt sequences
                    rrna['nt'] = nt

                    rrnas.append(rrna)
                    log.info(
                        'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        rrna['sequence'], rrna['start'], rrna['stop'], rrna['strand'], rrna['gene'], rrna['product'], length, coverage, truncated, score, evalue, nt[:10], nt[-10:]
                    )

    log.info('predicted=%i', len(rrnas))
    return rrnas

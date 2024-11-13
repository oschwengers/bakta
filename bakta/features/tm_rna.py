import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('TM_RNA')


def predict_tm_rnas(data: dict, sequences_path: Path):
    """Search for tmRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    cmd = [
        'aragorn',
        '-m',  # detect tmRNAs
        f'-gc{cfg.translation_table}',
        '-w',  # batch mode
        '-o', str(txt_output_path),
        str(sequences_path)
    ]
    if(cfg.complete):
        cmd.append('-c')  # complete circular sequence(s)
    else:
        cmd.append('-l')  # linear sequence(s)

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
        log.warning('tmRNAs failed! aragorn-error-code=%d', proc.returncode)
        raise Exception(f'aragorn error! error code: {proc.returncode}')

    tmrnas = []
    sequences = {seq['id']: seq for seq in data['sequences']}
    with txt_output_path.open() as fh:
        sequence_id = None
        for line in fh:
            line = line.strip()
            cols = line.split()
            if(line[0] == '>'):
                sequence_id = cols[0][1:]
            elif(len(cols) == 5):
                (nr, type, location, tag_location, tag_aa) = line.split()
                strand = bc.STRAND_FORWARD
                if(location[0] == 'c'):
                    strand = bc.STRAND_REVERSE
                    location = location[1:]
                (start, stop) = location[1:-1].split(',')
                start = int(start)
                stop = int(stop)
                tag_start, tag_stop = [int(pos) for pos in tag_location.split(',')]

                if(start > 0 and stop > 0):  # prevent edge tmRNA on linear sequences
                    tmrna = OrderedDict()
                    tmrna['type'] = bc.FEATURE_TM_RNA
                    tmrna['sequence'] = sequence_id
                    tmrna['start'] = start
                    tmrna['stop'] = stop
                    tmrna['strand'] = strand
                    tmrna['gene'] = 'ssrA'
                    tmrna['product'] = 'transfer-messenger RNA, SsrA'
                    tmrna['db_xrefs'] = [so.SO_TMRNA.id]
                    tmrna['tag'] = {
                        'start': start + tag_start - 1,
                        'stop': start + tag_stop - 1,
                        'aa': tag_aa.replace('*', '')
                    }

                    nt = bu.extract_feature_sequence(tmrna, sequences[sequence_id])  # extract nt sequences
                    tmrna['nt'] = nt

                    tag = tmrna['tag']
                    tag_nt = bu.extract_feature_sequence({'start': tag['start'], 'stop': tag['stop'], 'strand': strand}, sequences[sequence_id])  # extract nt sequences
                    tag['nt'] = tag_nt

                    if(start > stop):
                        tmrna['edge'] = True  # mark tmRNA as edge feature

                    tmrnas.append(tmrna)
                    log.info(
                        'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s]',
                        tmrna['sequence'], tmrna['start'], tmrna['stop'], tmrna['strand'], tmrna['gene'], tmrna['product'], nt[:10], nt[-10:]
                    )
    log.info('predicted=%i', len(tmrnas))
    return tmrnas

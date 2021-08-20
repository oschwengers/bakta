import logging
import subprocess as sp

from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('TM_RNA')


def predict_tm_rnas(genome, contigs_path):
    """Search for tmRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    cmd = [
        'aragorn',
        '-m',  # detect tmRNAs
        f'-gc{cfg.translation_table}',
        '-w',  # batch mode
        '-o', str(txt_output_path),
        str(contigs_path)
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
    contigs = {c['id']: c for c in genome['contigs']}
    with txt_output_path.open() as fh:
        contig_id = None
        for line in fh:
            line = line.strip()
            cols = line.split()
            if(line[0] == '>'):
                contig_id = cols[0][1:]
            elif(len(cols) == 5):
                (nr, type, location, tag_location, tag_aa) = line.split()
                strand = bc.STRAND_FORWARD
                if(location[0] == 'c'):
                    strand = bc.STRAND_REVERSE
                    location = location[1:]
                (start, stop) = location[1:-1].split(',')
                start = int(start)
                stop = int(stop)

                tmrna = OrderedDict()
                tmrna['type'] = bc.FEATURE_TM_RNA
                tmrna['contig'] = contig_id
                tmrna['start'] = start
                tmrna['stop'] = stop
                tmrna['strand'] = strand
                tmrna['gene'] = 'ssrA'
                tmrna['product'] = 'transfer-messenger RNA, SsrA'
                tmrna['db_xrefs'] = [so.SO_TMRNA.id]

                nt = bu.extract_feature_sequence(tmrna, contigs[contig_id])  # extract nt sequences
                tmrna['nt'] = nt

                if(start > stop):
                    tmrna['edge'] = True  # mark tmRNA as edge feature

                tmrnas.append(tmrna)
                log.info(
                    'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s]',
                    tmrna['contig'], tmrna['start'], tmrna['stop'], tmrna['strand'], tmrna['gene'], tmrna['product'], nt[:10], nt[-10:]
                )
    log.info('predicted=%i', len(tmrnas))
    return tmrnas

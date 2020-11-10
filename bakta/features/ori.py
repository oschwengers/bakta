
import logging
import subprocess as sp
from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('ORI')


def predict_oris(genome, contigs_path, ori_type):
    """Search for oriT/C sequences."""

    database = 'oric.fna' if ori_type == bc.FEATURE_ORIC else 'orit.fna'
    contigs = {c['id']: c for c in genome['contigs']}
    output_path = cfg.tmp_path.joinpath('ori.blastn.tsv')
    cmd = [
        'blastn',
        '-query', str(cfg.db_path.joinpath(database)),
        '-subject', str(contigs_path),
        '-culling_limit', '1',
        '-evalue', '1E-5',
        '-num_threads', str(cfg.threads),
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send length nident sstrand',
        '-out', str(output_path)
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
        log.warning('%s failed! blastn-error-code=%d', ori_type, proc.returncode)
        raise Exception(f'blastn error! error code: {proc.returncode}')

    # parse raw hits
    hits = {}
    with output_path.open() as fh:
        for line in fh:
            cols = line.strip().split('\t')
            hit = {
                'ori_id': cols[0],
                'ori_start': int(cols[1]),
                'ori_end': int(cols[2]),
                'ori_length': int(cols[3]),
                'contig': cols[4],
                'contig_start': int(cols[5]),
                'contig_stop': int(cols[6]),
                'strand': '+' if cols[9] == 'plus' else '-',
                'coverage': int(cols[7]) / int(cols[3]),
                'identity': int(cols[8]) / int(cols[7])
            }
            if(hit['strand'] == '-'):
                hit['contig_start'], hit['contig_stop'] = hit['contig_stop'], hit['contig_start']
            if(hit['coverage'] >= 0.8 and hit['identity'] >= 0.8):
                contig_hits = hits.get(hit['contig'], [])
                contig_hits.append(hit)
                if(len(contig_hits) == 1):
                    hits[hit['contig']] = contig_hits
                log.debug(
                    'raw hit: type=%s, contig=%s, start=%i, stop=%i, strand=%s, coverage=%f, identity=%f',
                    ori_type, hit['contig'], hit['contig_start'], hit['contig_stop'], hit['strand'], hit['coverage'], hit['identity']
                )
    
    # combine overlapping hits (simple 1D array peak detection)
    oris = []
    for contig in genome['contigs']:
        contig_hits = hits.get(contig['id'], None)
        if(contig_hits):
            region_hits = [0] * (contig['length'] + 1)  # init with extra leading slot (start at 1)
            for hit in contig_hits:
                for i in range(hit['contig_start'], hit['contig_stop'] + 1):
                    region_hits[i] += 1
            start = -1
            stop = -1
            for i, hit_count in enumerate(region_hits):  # detect regions of overlapping hits to combine overlapping hits
                if(hit_count == 0):
                    if(start != -1):  # new stop
                        stop = i - 1
                        if(ori_type == bc.FEATURE_ORIC and contig['type'] == bc.REPLICON_PLASMID):
                            ori_type = bc.FEATURE_ORIV
                        ori = OrderedDict()
                        ori['type'] = ori_type
                        ori['contig'] = contig['id']
                        ori['start'] = start
                        ori['stop'] = stop
                        ori['strand'] = bc.STRAND_UNKNOWN
                        oris.append(ori)
                        (refined_start, refined_stop) = refine_ori_region(region_hits, ori)
                        ori['start'] = refined_start
                        ori['stop'] = refined_stop
                        log.info(
                            'type=%s, contig=%s, start=%i, stop=%i',
                            ori_type, ori['contig'], ori['start'], ori['stop']
                        )
                        start = -1
                        stop = -1
                else:
                    if(start == -1):  # new start
                        start = i
    
    log.info('predicted=%i', len(oris))
    return oris

def refine_ori_region(region_hits, ori):
    log.debug('refine ori: [%i-%i]', ori['start'], ori['stop'])
    hit_region = region_hits[ori['start']:ori['stop'] + 1]
    hit_min = min(hit_region)
    hit_max = max(hit_region)
    hit_mean = sum(hit_region) / len(hit_region)
    log.debug('ori hit peak height: min=%i, max=%i, mean=%f', hit_min, hit_max, hit_mean)
    start = -1
    stop = -1
    for i in range(ori['start'], ori['stop'] + 2):  # extend range by extra bp to detect stop at 3' flank
        hit_value = region_hits[i]
        if(hit_value < hit_max / 2):
            if(start != -1):  # new stop
                stop = i - 1
                log.debug('new stop=%i', stop)
                return (start, stop)
        else:
            if(start == -1):  # new start
                start = i
                log.debug('new start=%i', start)
    return (start, ori['stop'])

import logging
import re

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:assembly_gap')

re_assembly_gap = re.compile(r'N{1,}', flags=0)

def detect_assembly_gaps(contigs):
    assembly_gaps = []
    for contig in contigs:
        m = re_assembly_gap.search(contig['sequence'])
        while m:
            start, end = m.span()
            assembly_gap = {
                'type': bc.FEATURE_GAP,
                'contig': contig['id'],
                'start': start,
                'stop': end + 1,
                'strand': '+',
                'length': end - start
            }
            assembly_gaps.append(assembly_gap)
            m = re_assembly_gap.search(contig['sequence'], end + 1)
    return assembly_gaps

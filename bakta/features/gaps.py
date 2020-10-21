
import logging
import re
from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:assembly_gap')

re_assembly_gap = re.compile(r'N{1,}', flags=0)

def detect_assembly_gaps(contigs):
    gaps = []
    for contig in contigs:
        m = re_assembly_gap.search(contig['sequence'])
        while m:
            start, end = m.span()
            
            gap = OrderedDict()
            gap['type'] = bc.FEATURE_GAP
            gap['contig'] = contig['id']
            gap['start'] = start
            gap['stop'] = end + 1
            gap['strand'] = '+'
            gap['length'] = end - start
            
            gaps.append(gap)
            m = re_assembly_gap.search(contig['sequence'], end + 1)
    return gaps

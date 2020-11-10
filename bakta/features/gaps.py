
import logging
import re
from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('GAP')

re_assembly_gap = re.compile(r'N{1,}', flags=0)

def detect_assembly_gaps(genome):
    gaps = []
    for contig in genome['contigs']:
        m = re_assembly_gap.search(contig['sequence'])
        while m:
            start, end = m.span()
            
            gap = OrderedDict()
            gap['type'] = bc.FEATURE_GAP
            gap['contig'] = contig['id']
            gap['start'] = start
            gap['stop'] = end + 1
            gap['strand'] = bc.STRAND_NA
            gap['length'] = end - start
            
            gaps.append(gap)
            log.info(
                'contig=%s, start=%i, stop=%i, length=%s',
                gap['contig'], gap['start'], gap['stop'], gap['length']
            )
            m = re_assembly_gap.search(contig['sequence'], end + 1)
    return gaps

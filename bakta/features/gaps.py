import logging
import re

from collections import OrderedDict

import bakta.constants as bc


log = logging.getLogger('GAP')
RE_ASSEMBLY_GAP = re.compile(r'N{1,}', flags=0)


def detect_assembly_gaps(genome):
    gaps = []
    for contig in genome['contigs']:
        m = RE_ASSEMBLY_GAP.search(contig['sequence'])
        while m:
            start, end = m.span()

            gap = OrderedDict()
            gap['type'] = bc.FEATURE_GAP
            gap['contig'] = contig['id']
            gap['start'] = start + 1
            gap['stop'] = end
            gap['strand'] = bc.STRAND_NA
            gap['length'] = end - start

            gaps.append(gap)
            log.info(
                'contig=%s, start=%i, stop=%i, length=%s',
                gap['contig'], gap['start'], gap['stop'], gap['length']
            )
            m = RE_ASSEMBLY_GAP.search(contig['sequence'], end + 1)
    return gaps

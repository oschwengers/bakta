import logging
import re

from collections import OrderedDict
from typing import Sequence

import bakta.constants as bc


log = logging.getLogger('GAP')
RE_ASSEMBLY_GAP = re.compile(r'N{1,}', flags=0)


def detect_assembly_gaps(data: dict) -> Sequence[dict]:
    gaps = []
    for seq in data['sequences']:
        m = RE_ASSEMBLY_GAP.search(seq['nt'])
        while m:
            start, end = m.span()

            gap = OrderedDict()
            gap['type'] = bc.FEATURE_GAP
            gap['sequence'] = seq['id']
            gap['start'] = start + 1
            gap['stop'] = end
            gap['strand'] = bc.STRAND_NA
            gap['length'] = end - start

            gaps.append(gap)
            log.info(
                'seq=%s, start=%i, stop=%i, length=%s',
                gap['sequence'], gap['start'], gap['stop'], gap['length']
            )
            m = RE_ASSEMBLY_GAP.search(seq['nt'], end + 1)
    return gaps

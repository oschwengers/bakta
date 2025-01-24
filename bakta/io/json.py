import json
import logging

from collections import OrderedDict
from pathlib import Path
from typing import Sequence

import bakta
import bakta.constants as bc
import bakta.config as cfg


log = logging.getLogger('JSON')


def write_json(data: dict, features: Sequence[dict], json_path: Path):
    log.info('write JSON: path=%s', json_path)

    # clean feature attributes
    for feat in features:
        if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):
            feat.pop('aa_digest')  # remove binary aa digest before JSON serialization

            # remove redundant IPS Dbxrefs
            ips = feat.get('ips', None)
            if(ips):
                ips.pop('db_xrefs')

            # remove redundant PSC Dbxrefs
            psc = feat.get('psc', None)
            if(psc):
                psc.pop('db_xrefs')

    version = OrderedDict()
    version['bakta'] = cfg.version
    version['db'] = {
        'version': f"{cfg.db_info['major']}.{cfg.db_info['minor']}",
        'type': cfg.db_info['type']
    }
    data['version'] = version

    with json_path.open('wt') as fh:
        json.dump(data, fh, indent=4)

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

    # replace features type dict by sorted feature list
    output = OrderedDict()
    if data is not None:
        ordered_genome = OrderedDict()
        ordered_genome['genus'] = data['genus']
        ordered_genome['species'] = data['species']
        ordered_genome['strain'] = data['strain']
        if('plasmid' in data):
            ordered_genome['plasmid'] = data['plasmid']
        ordered_genome['complete'] = data['complete']
        ordered_genome['gram'] = data['gram']
        ordered_genome['translation_table'] = data['translation_table']
        output['genome'] = ordered_genome

        stats = OrderedDict()
        stats['no_sequences'] = len(data['sequences'])
        stats['size'] = data['size']
        stats['gc'] = data['gc']
        stats['n_ratio'] = data['n_ratio']
        stats['n50'] = data['n50']
        stats['coding_ratio'] = data['coding_ratio']
        output['stats'] = stats

    output['features'] = features
    if data is not None:
        output['sequences'] = data['sequences']

    run = OrderedDict()
    run['start'] = cfg.run_start.strftime('%Y-%m-%d %H:%M:%S')
    run['end'] = cfg.run_end.strftime('%Y-%m-%d %H:%M:%S')
    output['run'] = run

    version = OrderedDict()
    version['bakta'] = bakta.__version__
    version['db'] = {
        'version': f"{cfg.db_info['major']}.{cfg.db_info['minor']}",
        'type': cfg.db_info['type']
    }
    output['version'] = version

    with json_path.open('wt') as fh:
        json.dump(output, fh, indent=4)

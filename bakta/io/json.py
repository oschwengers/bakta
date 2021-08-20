import json
import logging

from collections import OrderedDict

import bakta
import bakta.constants as bc
import bakta.config as cfg


log = logging.getLogger('JSON')


def write_json(genome, features, json_path):
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
    ordered_genome = OrderedDict()
    ordered_genome['genus'] = genome['genus']
    ordered_genome['species'] = genome['species']
    ordered_genome['strain'] = genome['strain']
    if('plasmid' in genome):
        ordered_genome['plasmid'] = genome['plasmid']
    ordered_genome['complete'] = genome['complete']
    ordered_genome['gram'] = genome['gram']
    ordered_genome['translation_table'] = genome['translation_table']
    output['genome'] = ordered_genome

    stats = OrderedDict()
    stats['no_sequences'] = len(genome['contigs'])
    stats['size'] = genome['size']
    stats['gc'] = genome['gc']
    stats['n_ratio'] = genome['n_ratio']
    stats['n50'] = genome['n50']
    stats['coding_ratio'] = genome['coding_ratio']
    output['stats'] = stats

    output['features'] = features
    output['sequences'] = genome['contigs']

    version = OrderedDict()
    version['bakta'] = bakta.__version__
    version['db'] = f"{cfg.db_info['major']}.{cfg.db_info['minor']}"
    output['version'] = version

    with json_path.open('wt') as fh:
        json.dump(output, fh, indent=4)

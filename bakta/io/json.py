import json
from collections import OrderedDict

import bakta.constants as bc

def write_json(genome, features, json_path):
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
    ordered_genome = OrderedDict()
    ordered_genome['genus'] = genome['genus']
    ordered_genome['species'] = genome['species']
    ordered_genome['strain'] = genome['strain']
    if('plasmid' in genome):
        ordered_genome['plasmid'] = genome['plasmid']
    ordered_genome['gram'] = genome['gram']
    ordered_genome['translation_table'] = genome['translation_table']
    
    ordered_genome['no_sequences'] = len(genome['contigs'])
    ordered_genome['size'] = genome['size']
    ordered_genome['gc'] = genome['gc']
    ordered_genome['n_ratio'] = genome['n_ratio']
    ordered_genome['n50'] = genome['n50']
    ordered_genome['coding_ratio'] = genome['coding_ratio']

    ordered_genome['features'] = features
    ordered_genome['sequences'] = genome['contigs']

    with json_path.open('w') as fh:
        json.dump(ordered_genome, fh, indent=4)
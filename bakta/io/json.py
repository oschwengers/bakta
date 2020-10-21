import json

import bakta.constants as bc

def write_json(features, json_path):
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
    
    with json_path.open('w') as fh:
        json.dump(features, fh, indent=4)
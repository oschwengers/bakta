import binascii
import json

import bakta.constants as bc

def write_json(features, json_path):

    for feat in features:
        if(feat['type'] == bc.FEATURE_CDS):
            feat['aa_hash'] = binascii.b2a_hex(feat['aa_hash']).decode()
    
    with json_path.open('w') as fh:
        json.dump(features, fh, sort_keys=True, indent=4)
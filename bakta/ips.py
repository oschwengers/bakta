
import logging
import sqlite3

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu

############################################################################
# IPS DB columns
############################################################################
DB_IPS_COL_UNIREF90 = 'uniref90_id'
DB_IPS_COL_UNIREF100 = 'uniref100_id'
DB_IPS_COL_GENE = 'gene'
DB_IPS_COL_PRODUCT = 'product'
DB_IPS_COL_EC = 'ec_ids'
DB_IPS_COL_GO = 'go_ids'

log = logging.getLogger('IPS')


def lookup(features):
    """Lookup IPS by hash values."""
    try:
        features_found = []
        features_not_found = []
        with sqlite3.connect(f"file:{cfg.db_path.joinpath('bakta.db')}?mode=ro", uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                if('truncated' in feature):  # skip truncated CDS
                    continue
                uniref100_id = feature.get('ups', {}).get('uniref100_id', None)
                if(uniref100_id):
                    c.execute('select * from ips where uniref100_id=?', (feature['ups']['uniref100_id'][10:],))
                    rec = c.fetchone()
                    if(rec is not None):
                        ips = parse_annotation(rec)
                        feature['ips'] = ips
                        features_found.append(feature)
                        log.debug(
                            'lookup: contig=%s, start=%i, stop=%i, aa-length=%i, strand=%s, gene=%s, UniRef100=%s, UniRef90=%s',
                            feature['contig'], feature['start'], feature['stop'], len(feature['sequence']), feature['strand'], ips.get(DB_IPS_COL_GENE, ''), ips.get(DB_IPS_COL_UNIREF100, ''), ips.get(DB_IPS_COL_UNIREF90, '')
                        )
                    else:
                        features_not_found.append(feature)
                else:
                    features_not_found.append(feature)

        log.info('looked-up=%i', len(features_found))
        return features_found, features_not_found
    except Exception as ex:
        log.exception('Could not read IPSs from db!', ex)
        raise Exception('SQL error!', ex)


def parse_annotation(rec):
    ips = {
        DB_IPS_COL_UNIREF100: bc.DB_PREFIX_UNIREF_100 + rec[DB_IPS_COL_UNIREF100]  # must not be NULL/None
    }
    db_xrefs = [
        'SO:0001217',
        f'{bc.DB_XREF_UNIREF_100}:{ips[DB_IPS_COL_UNIREF100]}'
    ]

    # add non-empty PSC annotations and attach database prefixes to identifiers
    if(rec[DB_IPS_COL_GENE]):
        ips[DB_IPS_COL_GENE] = rec[DB_IPS_COL_GENE]
    if(rec[DB_IPS_COL_PRODUCT]):
        ips[DB_IPS_COL_PRODUCT] = rec[DB_IPS_COL_PRODUCT]
    if(rec[DB_IPS_COL_UNIREF90]):
        ips[DB_IPS_COL_UNIREF90] = bc.DB_PREFIX_UNIREF_90 + rec[DB_IPS_COL_UNIREF90]
        db_xrefs.append(f'{bc.DB_XREF_UNIREF_90}:{ips[DB_IPS_COL_UNIREF90]}')
    if(rec[DB_IPS_COL_EC]):
        ips[DB_IPS_COL_EC] = rec[DB_IPS_COL_EC]
        ecs = []
        for ec in rec[DB_IPS_COL_EC].split(','):
            if(ec != ''):
                ecs.append(ec)
                db_xrefs.append(f'{bc.DB_XREF_EC}:{ec}')
        if(len(ecs) != 0):
            ips[DB_IPS_COL_EC] = ecs
    if(rec[DB_IPS_COL_GO]):
        go_ids = []
        for go_id in rec[DB_IPS_COL_GO].split(','):
            if(go_id != ''):
                go_id = bc.DB_PREFIX_GO + go_id
                go_ids.append(go_id)
                db_xrefs.append(go_id)
        if(len(go_ids) != 0):
            ips[DB_IPS_COL_GO] = go_ids
    
    ips['db_xrefs'] = db_xrefs
    return ips

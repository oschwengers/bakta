
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

log = logging.getLogger('ips')


def lookup(features):
    """Lookup IPS by hash values."""
    try:
        features_found = []
        features_not_found = []
        with sqlite3.connect("file:%s?mode=ro" % str(cfg.db_path.joinpath('bakta.db')), uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                uniref100_id = feature.get('ups', {}).get('uniref100_id', None)
                if(uniref100_id):
                    c.execute("select * from ips where uniref100_id=?", (feature['ups']['uniref100_id'][10:],))
                    rec = c.fetchone()
                    if(rec is not None):
                        ips = parse_annotation(rec)
                        feature['ips'] = ips

                        if('db_xrefs' not in feature):
                            feature['db_xrefs'] = []
                        db_xrefs = feature['db_xrefs']
                        db_xrefs.append('SO:0001217')
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_100, ips[DB_IPS_COL_UNIREF100]))
                        if(bu.has_annotation(ips, DB_IPS_COL_UNIREF90)):
                            db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_90, ips[DB_IPS_COL_UNIREF90]))
                        if(bu.has_annotation(ips, DB_IPS_COL_GO)):
                            db_xrefs.append('%s:%s' % (bc.DB_XREF_GO, ips[DB_IPS_COL_GO]))
                        if(bu.has_annotation(ips, DB_IPS_COL_EC)):
                            db_xrefs.append('%s:%s' % (bc.DB_XREF_EC, ips[DB_IPS_COL_EC]))
                        features_found.append(feature)

                        log.debug(
                            'lookup: contig=%s, start=%i, stop=%i, aa-length=%i, strand=%s, gene=%s, UniRef100=%s, UniRef90=%s',
                            feature['contig'], feature['start'], feature['stop'], len(feature['sequence']), feature['strand'], ips.get(DB_IPS_COL_GENE, ''), ips.get(DB_IPS_COL_UNIREF100, ''), ips.get(DB_IPS_COL_UNIREF90, '')
                        )
                    else:
                        features_not_found.append(feature)
                else:
                    features_not_found.append(feature)

        log.info('# %i', len(features_found))
        return features_found, features_not_found
    except Exception as ex:
        log.exception('Could not read IPSs from db!', ex)
        raise Exception("SQL error!", ex)


def parse_annotation(rec):
    ips = {
        DB_IPS_COL_UNIREF100: bc.DB_PREFIX_UNIREF_100 + rec[DB_IPS_COL_UNIREF100]  # must not be NULL/None
    }

    # add non-empty PSC annotations and attach database prefixes to identifiers
    if(rec[DB_IPS_COL_GENE]):
        ips[DB_IPS_COL_GENE] = rec[DB_IPS_COL_GENE]
    if(rec[DB_IPS_COL_PRODUCT]):
        ips[DB_IPS_COL_PRODUCT] = rec[DB_IPS_COL_PRODUCT]
    if(rec[DB_IPS_COL_UNIREF90]):
        ips[DB_IPS_COL_UNIREF90] = bc.DB_PREFIX_UNIREF_90 + rec[DB_IPS_COL_UNIREF90]
    if(rec[DB_IPS_COL_EC]):
        ips[DB_IPS_COL_EC] = rec[DB_IPS_COL_EC]
    if(rec[DB_IPS_COL_GO]):
        go_ids = []
        for go_id in rec[DB_IPS_COL_GO].split(';'):
            if(go_id != ''):
                go_ids.append(bc.DB_PREFIX_GO + go_id)
        if(len(go_ids) != 0):
            ips[DB_PSC_COL_GO] = go_ids
    
    return ips

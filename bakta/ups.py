
import logging
import sqlite3

import bakta.constants as bc
import bakta.utils as bu

############################################################################
# UPS DB columns
############################################################################
DB_UPS_COL_HASH = 'hash'
DB_UPS_COL_LENGTH = 'length'
DB_UPS_COL_UNIREF90 = 'uniref90_id'
DB_UPS_COL_UNIREF100 = 'uniref100_id'
DB_UPS_COL_UNIPARC = 'uniparc_id'
DB_UPS_COL_REFSEQ_NRP = 'ncbi_nrp_id'
DB_UPS_COL_UNIPROTKB = 'uniprotkb_acc'
DB_UPS_COL_GENE = 'gene'
DB_UPS_COL_PRODUCT = 'product'

log = logging.getLogger('ups')


def lookup_upss(config, features):
    """Lookup UPS by hash values."""
    try:
        features_found = []
        features_not_found = []
        with sqlite3.connect("file:%s?mode=ro" % str(config['db'].joinpath('ups.db')), uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                rec = c.fetchone("select * from ups where hash=?", (feature['aa_hash'],))
                if(rec is not None):
                    ups = {
                        'uniref100_id': bc.DB_PREFIX_UNIREF_100 + rec[DB_UPS_COL_UNIREF100],  # must not be NULL/None
                        'uniprotkb_acc': rec[DB_UPS_COL_UNIPROTKB],
                        'gene': rec[DB_UPS_COL_GENE],
                        'product': rec[DB_UPS_COL_PRODUCT]
                    }

                    # reattach database prefixes to identifiers
                    if(rec[DB_UPS_COL_UNIREF90] is not None):
                        ups['uniref90_id'] = bc.DB_PREFIX_UNIREF_90 + rec[DB_UPS_COL_UNIREF90]
                    if(rec[DB_UPS_COL_UNIPARC] is not None):
                        ups['uniparc_id'] = bc.DB_PREFIX_UNIPARC + rec[DB_UPS_COL_UNIPARC]
                    if(rec[DB_UPS_COL_REFSEQ_NRP] is not None):
                        ups['refseq_nrp_id'] = bc.DB_PREFIX_REFSEQ_NRP + rec[DB_UPS_COL_REFSEQ_NRP]

                    feature['ups'] = ups

                    if('db_xrefs' not in feature):
                        feature['db_xrefs'] = []
                    db_xrefs = feature['db_xrefs']
                    db_xrefs.append('SO:0001217')
                    db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_100, ups['uniref100_id']))
                    if(bu.has_annotation(ups, 'cluster_id')):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_90, ups['uniref90_id']))
                    if(bu.has_annotation(ups, 'uniparc_id')):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIPARC, ups['uniparc_id']))
                    if(bu.has_annotation(ups, 'ncbi_nrp_id')):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_REFSEQ_NRP, ups['ncbi_nrp_id']))
                    if(bu.has_annotation(ups, 'uniprotkb_acc')):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIPROTKB, ups['uniprotkb_acc']))
                    features_found.append(feature)

                    log.debug(
                        'UPS: contig=%s, start=%i, stop=%i, strand=%s, UniRef100=%s, NCBI NRP=%s, UniRef90=%s',
                        feature['contig'], feature['gene'], feature['start'], feature['stop'], feature['strand'], ups['uniref100_id'], ups['ncbi_nrp_id'], ups['cluster_id']
                    )
                else:
                    features_not_found.append(feature)

        log.info('UPSs: # %i', len(features_found))
        return features_found, features_not_found
    except Exception:
        log.exception('Could not read UPSs from db!')

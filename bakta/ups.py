
import logging
import sqlite3

import bakta.config as cfg
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


def lookup_upss(features):
    """Lookup UPS by hash values."""
    try:
        features_found = []
        features_not_found = []
        with sqlite3.connect("file:%s?mode=ro" % str(cfg.db_path.joinpath('bakta.db')), uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                c.execute("select * from ups where hash=?", (feature['aa_hash'],))
                rec = c.fetchone()
                if(rec is not None and rec[DB_UPS_COL_LENGTH] == len(feature['sequence'])):
                    ups = {
                        DB_UPS_COL_UNIREF100: bc.DB_PREFIX_UNIREF_100 + rec[DB_UPS_COL_UNIREF100],  # must not be NULL/None
                        DB_UPS_COL_UNIPROTKB: rec[DB_UPS_COL_UNIPROTKB],
                        DB_UPS_COL_GENE: rec[DB_UPS_COL_GENE],
                        DB_UPS_COL_PRODUCT: rec[DB_UPS_COL_PRODUCT],
                        DB_UPS_COL_UNIREF90: '',
                        DB_UPS_COL_UNIPARC: '',
                        DB_UPS_COL_REFSEQ_NRP: '',
                    }

                    # reattach database prefixes to identifiers
                    if(rec[DB_UPS_COL_UNIREF90] is not None):
                        ups[DB_UPS_COL_UNIREF90] = bc.DB_PREFIX_UNIREF_90 + rec[DB_UPS_COL_UNIREF90]
                    if(rec[DB_UPS_COL_UNIPARC] is not None):
                        ups[DB_UPS_COL_UNIPARC] = bc.DB_PREFIX_UNIPARC + rec[DB_UPS_COL_UNIPARC]
                    if(rec[DB_UPS_COL_REFSEQ_NRP] is not None):
                        ups[DB_UPS_COL_REFSEQ_NRP] = bc.DB_PREFIX_REFSEQ_NRP + rec[DB_UPS_COL_REFSEQ_NRP]

                    feature['ups'] = ups

                    if('db_xrefs' not in feature):
                        feature['db_xrefs'] = []
                    db_xrefs = feature['db_xrefs']
                    db_xrefs.append('SO:0001217')
                    db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_100, ups[DB_UPS_COL_UNIREF100]))
                    if(bu.has_annotation(ups, DB_UPS_COL_UNIREF90)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_90, ups[DB_UPS_COL_UNIREF90]))
                    if(bu.has_annotation(ups, DB_UPS_COL_UNIPARC)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIPARC, ups[DB_UPS_COL_UNIPARC]))
                    if(bu.has_annotation(ups, DB_UPS_COL_REFSEQ_NRP)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_REFSEQ_NRP, ups[DB_UPS_COL_REFSEQ_NRP]))
                    if(bu.has_annotation(ups, DB_UPS_COL_UNIPROTKB)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIPROTKB, ups[DB_UPS_COL_UNIPROTKB]))
                    features_found.append(feature)

                    log.debug(
                        'UPS: contig=%s, start=%i, stop=%i, aa-length=%i, strand=%s, gene=%s, UniRef100=%s, NCBI NRP=%s, UniRef90=%s',
                        feature['contig'], feature['start'], feature['stop'], len(feature['sequence']), feature['strand'], ups[DB_UPS_COL_GENE], ups[DB_UPS_COL_UNIREF100], ups[DB_UPS_COL_REFSEQ_NRP], ups[DB_UPS_COL_UNIREF90]
                    )
                else:
                    features_not_found.append(feature)

        log.info('UPSs: # %i', len(features_found))
        return features_found, features_not_found
    except Exception as ex:
        log.exception('Could not read UPSs from db!', ex)
        raise Exception("SQL error!", ex)


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
DB_UPS_COL_UNIPARC = 'uniparc_id'
DB_UPS_COL_REFSEQ_NRP = 'ncbi_nrp_id'
DB_UPS_COL_UNIREF100 = 'uniref100_id'

log = logging.getLogger('UPS')


def lookup(features):
    """Lookup UPS by hash values."""
    try:
        features_found = []
        features_not_found = []
        with sqlite3.connect(f"file:{cfg.db_path.joinpath('bakta.db')}?mode=ro", uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                if('truncated' in feature):  # skip truncated CDS
                    continue
                c.execute('select * from ups where hash=?', (feature['aa_digest'],))
                rec = c.fetchone()
                if(rec is not None and rec[DB_UPS_COL_LENGTH] == len(feature['sequence'])):
                    ups = parse_annotation(rec)
                    feature['ups'] = ups
                    features_found.append(feature)
                    log.debug(
                        'lookup: contig=%s, start=%i, stop=%i, aa-length=%i, strand=%s, UniParc=%s, UniRef100=%s, NCBI NRP=%s',
                        feature['contig'], feature['start'], feature['stop'], len(feature['sequence']), feature['strand'], ups.get(DB_UPS_COL_UNIPARC, ''), ups.get(DB_UPS_COL_UNIREF100, ''), ups.get(DB_UPS_COL_REFSEQ_NRP, '')
                    )
                else:
                    features_not_found.append(feature)

        log.info('looked-up=%i', len(features_found))
        return features_found, features_not_found
    except Exception as ex:
        log.exception('Could not read UPSs from db!', ex)
        raise Exception("SQL error!", ex)


def parse_annotation(rec):
    ups = {}
    db_xrefs = ['SO:0001217']

    # add non-empty PSC annotations and attach database prefixes to identifiers
    if(rec[DB_UPS_COL_UNIPARC]):
        ups[DB_UPS_COL_UNIPARC] = bc.DB_PREFIX_UNIPARC + rec[DB_UPS_COL_UNIPARC]
        db_xrefs.append(f'{bc.DB_XREF_UNIPARC}:{ups[DB_UPS_COL_UNIPARC]}')
    if(rec[DB_UPS_COL_REFSEQ_NRP]):
        ups[DB_UPS_COL_REFSEQ_NRP] = bc.DB_PREFIX_REFSEQ_NRP + rec[DB_UPS_COL_REFSEQ_NRP]
        db_xrefs.append(f'{bc.DB_XREF_REFSEQ_NRP}:{ups[DB_UPS_COL_REFSEQ_NRP]}')
    if(rec[DB_UPS_COL_UNIREF100]):
        ups[DB_UPS_COL_UNIREF100] = bc.DB_PREFIX_UNIREF_100 + rec[DB_UPS_COL_UNIREF100]
        db_xrefs.append(f'{bc.DB_XREF_UNIREF_100}:{ups[DB_UPS_COL_UNIREF100]}')
    
    ups['db_xrefs'] = db_xrefs
    return ups

import logging
import sqlite3

from concurrent.futures import ThreadPoolExecutor

import bakta.config as cfg
import bakta.constants as bc


############################################################################
# PSCC DB columns
############################################################################
DB_PSCC_COL_UNIREF50 = 'uniref50_id'
DB_PSCC_COL_PRODUCT = 'product'


log = logging.getLogger('PSCC')


def lookup(features):
    """Lookup PSCC information"""
    no_pscc_lookups = 0
    try:
        rec_futures = []
        with sqlite3.connect(f"file:{cfg.db_path.joinpath('bakta.db')}?mode=ro", uri=True, check_same_thread=False) as conn:
            conn.execute('PRAGMA omit_readlock;')
            conn.row_factory = sqlite3.Row
            with ThreadPoolExecutor(max_workers=max(10, cfg.threads)) as tpe:  # use min 10 threads for IO bound non-CPU lookups
                for feature in features:
                    if('psc' in feature):
                        uniref50_id = feature['psc'].get(DB_PSCC_COL_UNIREF50, None)
                        if(uniref50_id is not None):
                            if(bc.DB_PREFIX_UNIREF_50 in uniref50_id):
                                uniref50_id = uniref50_id[9:]  # remove 'UniRef50_' prefix
                            future = tpe.submit(fetch_db_pscc_result, conn, uniref50_id)
                            rec_futures.append((feature, future))

        for (feature, future) in rec_futures:
            rec = future.result()
            if(rec is not None):
                pscc = parse_annotation(rec)
                feature['pscc'] = pscc
                no_pscc_lookups += 1
                log.debug(
                    'lookup: contig=%s, start=%i, stop=%i, strand=%s, UniRef50=%s, product=%s',
                    feature['contig'], feature['start'], feature['stop'], feature['strand'], pscc.get(DB_PSCC_COL_UNIREF50, ''), pscc.get(DB_PSCC_COL_PRODUCT, '')
                )
            else:
                log.debug('lookup failed! uniref50_id=%s', uniref50_id)
    except Exception as ex:
        log.exception('Could not read PSCCs from db!', ex)
        raise Exception('SQL error!', ex)
    log.info('looked-up=%i', no_pscc_lookups)


def fetch_db_pscc_result(conn, uniref50_id):
    c = conn.cursor()
    c.execute('select * from pscc where uniref50_id=?', (uniref50_id,))
    rec = c.fetchone()
    c.close()
    return rec


def parse_annotation(rec):
    uniref_full_id = bc.DB_PREFIX_UNIREF_50 + rec[DB_PSCC_COL_UNIREF50]
    pscc = {
        DB_PSCC_COL_UNIREF50: uniref_full_id,  # must not be NULL/None
        'db_xrefs': [
            'SO:0001217',
            f'{bc.DB_XREF_UNIPROTKB}:{uniref_full_id}'
        ]
    }
    # add non-empty PSCC annotations and attach database prefixes to identifiers
    if(rec[DB_PSCC_COL_PRODUCT]):
        pscc[DB_PSCC_COL_PRODUCT] = rec[DB_PSCC_COL_PRODUCT]
    return pscc

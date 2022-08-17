import logging
import subprocess as sp
import sqlite3

from concurrent.futures import ThreadPoolExecutor
from typing import Sequence, Tuple

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf


############################################################################
# PSC DB columns
############################################################################
DB_PSC_COL_UNIREF90 = 'uniref90_id'
DB_PSC_COL_GENE = 'gene'
DB_PSC_COL_PRODUCT = 'product'
DB_PSC_COL_UNIREF50 = 'uniref50_id'
DB_PSC_COL_COG_ID = 'cog_id'
DB_PSC_COL_COG_CAT = 'cog_category'
DB_PSC_COL_EC = 'ec_ids'
DB_PSC_COL_KOFAM_ID = 'kegg_orthology_id'
DB_PSC_COL_GO = 'go_ids'


log = logging.getLogger('PSC')


def search(cdss: Sequence[dict]) -> Tuple[Sequence[dict], Sequence[dict], Sequence[dict]]:
    """Conduct homology search of CDSs against PCS db."""
    cds_aa_path = cfg.tmp_path.joinpath('cds.psc.faa')
    orf.write_internal_faa(cdss, cds_aa_path)
    diamond_output_path = cfg.tmp_path.joinpath('diamond.cds.tsv')
    diamond_db_path = cfg.db_path.joinpath('psc.dmnd')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(diamond_db_path),
        '--query', str(cds_aa_path),
        '--out', str(diamond_output_path),
        '--id', str(int(bc.MIN_PSCC_IDENTITY * 100)),  # '50',
        '--query-cover', str(int(bc.MIN_PSC_COVERAGE * 100)),  # '80'
        '--subject-cover', str(int(bc.MIN_PSC_COVERAGE * 100)),  # '80'
        '--max-target-seqs', '1',  # single best output
        '--outfmt', '6',
        '--threads', str(cfg.threads),
        '--tmpdir', str(cfg.tmp_path),  # use tmp folder
        '--block-size', '3',  # slightly increase block size for faster executions
        '--fast'
    ]
    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('PSC failed! diamond-error-code=%d', proc.returncode)
        raise Exception(f'diamond error! error code: {proc.returncode}')

    cds_by_hexdigest = orf.get_orf_dictionary(cdss)
    with diamond_output_path.open() as fh:
        for line in fh:
            (aa_identifier, cluster_id, identity, alignment_length, align_mismatches,
                align_gaps, query_start, query_end, subject_start, subject_end,
                evalue, bitscore) = line.split('\t')
            cds = cds_by_hexdigest[aa_identifier]
            query_cov = int(alignment_length) / len(cds['aa'])
            identity = float(identity) / 100
            if(query_cov >= bc.MIN_PSC_COVERAGE and identity >= bc.MIN_PSCC_IDENTITY):
                cds['psc'] = {
                    DB_PSC_COL_UNIREF90: cluster_id,
                    'query_cov': query_cov,
                    'identity': identity,
                    'valid': identity >= bc.MIN_PSC_IDENTITY
                }
                log.debug(
                    'homology: contig=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, identity=%0.3f, UniRef90=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], len(cds['aa']), query_cov, identity, cluster_id
                )

    pscs_found = []
    psccs_found = []
    cds_not_found = []
    for cds in cdss:
        if('psc' in cds):
            if(cds['psc']['valid']):
                pscs_found.append(cds)
            else:
                psccs_found.append(cds)
        else:
            cds_not_found.append(cds)
    log.info('found: PSC=%i, PSCC=%i', len(pscs_found), len(psccs_found))
    return pscs_found, psccs_found, cds_not_found


def lookup(features: Sequence[dict], pseudo: bool = False):
    """Lookup PCS information"""
    no_psc_lookups = 0
    try:
        rec_futures = []
        with sqlite3.connect(f"file:{cfg.db_path.joinpath('bakta.db')}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
            conn.execute('PRAGMA omit_readlock;')
            conn.row_factory = sqlite3.Row
            with ThreadPoolExecutor(max_workers=max(10, cfg.threads)) as tpe:  # use min 10 threads for IO bound non-CPU lookups
                for feature in features:
                    uniref90_id = None
                    if(pseudo):  # if pseudogene use pseudogene info 
                        uniref90_id = feature[bc.PSEUDOGENE]['inference'][DB_PSC_COL_UNIREF90]
                    else:
                        if('psc' in feature):
                            uniref90_id = feature['psc'].get(DB_PSC_COL_UNIREF90, None)
                        elif('ips' in feature):
                            uniref90_id = feature['ips'].get(DB_PSC_COL_UNIREF90, None)

                    if(uniref90_id is not None):
                        if(bc.DB_PREFIX_UNIREF_90 in uniref90_id):
                            uniref90_id = uniref90_id[9:]  # remove 'UniRef90_' prefix
                        future = tpe.submit(fetch_db_psc_result, conn, uniref90_id)
                        rec_futures.append((feature, future))

        for (feature, future) in rec_futures:
            rec = future.result()
            if(rec is not None):
                psc = parse_annotation(rec)
                if(pseudo):
                    feature[bc.PSEUDOGENE]['psc'] = psc  # store PSC info for pseudogene annotation
                    feature[bc.PSEUDOGENE]['inference'][DB_PSC_COL_UNIREF90] = bc.DB_PREFIX_UNIREF_90 + feature[bc.PSEUDOGENE]['inference'][DB_PSC_COL_UNIREF90]  # restore full UniRef90 identifier
                else:
                    if('psc' in feature):
                        feature['psc'] = {**feature['psc'], **psc}  # merge dicts to store alignment info for later PSC/PSCC annotations
                    else:
                        feature['psc'] = psc
                no_psc_lookups += 1
                log.debug(
                    'lookup: contig=%s, start=%i, stop=%i, strand=%s, UniRef90=%s, EC=%s, gene=%s, product=%s',
                    feature['contig'], feature['start'], feature['stop'], feature['strand'], psc.get(DB_PSC_COL_UNIREF90, ''), psc.get(DB_PSC_COL_EC, ''), psc.get(DB_PSC_COL_GENE, ''), psc.get(DB_PSC_COL_PRODUCT, '')
                )
            else:
                log.debug('lookup: ID not found! uniref90_id=%s', uniref90_id)
    except Exception as ex:
        log.exception('Could not read PSCs from db!', ex)
        raise Exception('SQL error!', ex)
    log.info('looked-up=%i', no_psc_lookups)


def fetch_db_psc_result(conn: sqlite3.Connection, uniref90_id: str):
    c = conn.cursor()
    c.execute('select * from psc where uniref90_id=?', (uniref90_id,))
    rec = c.fetchone()
    c.close()
    return rec


def parse_annotation(rec) -> dict:
    psc = {
        DB_PSC_COL_UNIREF90: bc.DB_PREFIX_UNIREF_90 + rec[DB_PSC_COL_UNIREF90]  # must not be NULL/None
    }
    db_xrefs = [
        'SO:0001217',
        f'{bc.DB_XREF_UNIREF}:{psc[DB_PSC_COL_UNIREF90]}'
    ]

    # add non-empty PSC annotations and attach database prefixes to identifiers
    if(rec[DB_PSC_COL_GENE]):
        psc[DB_PSC_COL_GENE] = rec[DB_PSC_COL_GENE]
    if(rec[DB_PSC_COL_PRODUCT]):
        psc[DB_PSC_COL_PRODUCT] = rec[DB_PSC_COL_PRODUCT]
    if(rec[DB_PSC_COL_EC]):
        ecs = []
        for ec in rec[DB_PSC_COL_EC].split(','):
            if(ec != ''):
                ecs.append(ec)
                db_xrefs.append(f'{bc.DB_XREF_EC}:{ec}')
        if(len(ecs) != 0):
            psc[DB_PSC_COL_EC] = ecs
    if(rec[DB_PSC_COL_UNIREF50]):
        psc[DB_PSC_COL_UNIREF50] = bc.DB_PREFIX_UNIREF_50 + rec[DB_PSC_COL_UNIREF50]
        db_xrefs.append(f'{bc.DB_XREF_UNIREF}:{psc[DB_PSC_COL_UNIREF50]}')
    if(rec[DB_PSC_COL_COG_ID]):
        psc[DB_PSC_COL_COG_ID] = bc.DB_PREFIX_COG + rec[DB_PSC_COL_COG_ID]
        db_xrefs.append(f'{bc.DB_XREF_COG}:{psc[DB_PSC_COL_COG_ID]}')
    if(rec[DB_PSC_COL_COG_CAT]):
        psc[DB_PSC_COL_COG_CAT] = rec[DB_PSC_COL_COG_CAT]
        db_xrefs.append(f'{bc.DB_XREF_COG}:{psc[DB_PSC_COL_COG_CAT]}')
    if(rec[DB_PSC_COL_KOFAM_ID]):
        psc[DB_PSC_COL_KOFAM_ID] = bc.DB_PREFIX_KEGG_ORTHOLOGY + rec[DB_PSC_COL_KOFAM_ID]
        db_xrefs.append(f'{bc.DB_XREF_KOFAM}:{psc[DB_PSC_COL_KOFAM_ID]}')
    if(rec[DB_PSC_COL_GO]):
        go_ids = []
        for go_id in rec[DB_PSC_COL_GO].split(','):
            if(go_id != ''):
                go_id = bc.DB_PREFIX_GO + go_id
                go_ids.append(go_id)
                db_xrefs.append(go_id)
        if(len(go_ids) != 0):
            psc[DB_PSC_COL_GO] = go_ids

    psc['db_xrefs'] = db_xrefs
    return psc

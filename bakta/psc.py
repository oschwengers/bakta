
import logging
import subprocess as sp
import sqlite3

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu

############################################################################
# PSC DB columns
############################################################################
DB_PSC_COL_UNIREF90 = 'uniref90_id'
DB_PSC_COL_UNIREF50 = 'uniref50_id'
DB_PSC_COL_COG_ID = 'cog_id'
DB_PSC_COL_COG_CAT = 'cog_category'
DB_PSC_COL_EC = 'ec_id'
DB_PSC_COL_GO = 'go_ids'
DB_PSC_COL_GENE = 'gene'
DB_PSC_COL_PRODUCT = 'product'

log = logging.getLogger('psc')


def search_pscs(cdss):
    """Conduct homology search of CDSs against PCS db."""
    cds_fasta_path = cfg.tmp_path.joinpath('cds.faa')
    with cds_fasta_path.open(mode='w') as fh:
        for cds in cdss:
            fh.write(">%s\n%s\n" % (cds['tmp_id'], cds['sequence']))
    diamond_output_path = cfg.tmp_path.joinpath('diamond.tsv')
    diamond_db_path = cfg.db_path.joinpath('psc')
    cmd = [
        'diamond',
        'blastp',
        '--db', str(diamond_db_path),
        '--query', str(cds_fasta_path),
        '--out', str(diamond_output_path),
        '--id', str(int(bc.MIN_PROTEIN_IDENTITY * 100)), # '90',
        '--query-cover', str(int(bc.MIN_PROTEIN_COVERAGE * 100)), # '80'
        '--subject-cover', str(int(bc.MIN_PROTEIN_COVERAGE * 100)), # '80'
        '--max-target-seqs', '1',  # single best output
        '--outfmt', '6',
        '--threads', str(cfg.threads)
    ]
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning(
            'PSC failed! diamond-error-code=%d', proc.returncode
        )
        log.debug(
            'PSC: cmd=%s stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("diamond error! error code: %i" % proc.returncode)

    cds_by_id = {cds['tmp_id']: cds for cds in cdss}
    with diamond_output_path.open() as fh:
        for line in fh:
            (cds_id, cluster_id, identity, alignment_length, align_mismatches,
                align_gaps, query_start, query_end, subject_start, subject_end,
                evalue, bitscore) = line.split('\t')
            cds = cds_by_id[int(cds_id)]
            query_cov = int(alignment_length) / len(cds['sequence'])
            identity = float(identity) / 100
            if(query_cov >= bc.MIN_PROTEIN_COVERAGE and identity >= bc.MIN_PROTEIN_IDENTITY):
                cds['psc'] = {
                    DB_PSC_COL_UNIREF90: cluster_id,
                    'query-cov': query_cov,
                    'identity': identity
                }
                log.debug(
                    'homology: contig=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, identity=%0.3f, UniRef90=%s',
                     cds['contig'], cds['start'], cds['stop'], cds['strand'], len(cds['sequence']), query_cov, identity, cluster_id
                )

    pscs_found = []
    pscs_not_found = []
    for cds in cdss:
        if('psc' in cds):
            pscs_found.append(cds)
        else:
            pscs_not_found.append(cds)
    log.info('PSCs: # %i', len(pscs_found))
    return pscs_found, pscs_not_found


def lookup_pscs(features):
    """Lookup PCS information"""
    try:
        with sqlite3.connect("file:%s?mode=ro" % str(cfg.db_path.joinpath('bakta.db')), uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for feature in features:
                if('psc' in feature):
                    uniref90_id = feature['psc'].get(DB_PSC_COL_UNIREF90, None)
                elif('ups' in feature):
                    uniref90_id = feature['ups'].get(DB_PSC_COL_UNIREF90, None)
                else:
                    continue  # skip PSC lookup for this feature object

                if(not uniref90_id):
                    continue  # skip PSC lookup for this feature without a valid UniRef90 id
                if(bc.DB_PREFIX_UNIREF_90 in uniref90_id):
                    uniref90_id = uniref90_id[9:]  # remove 'UniRef90_' prefix
                
                c.execute("select * from psc where uniref90_id=?", (uniref90_id,))
                rec = c.fetchone()
                if(rec is not None):
                    psc = parse_psc_annotation(rec)
                    feature['psc'] = psc
                    
                    if('db_xrefs' not in feature):
                        feature['db_xrefs'] = []
                    db_xrefs = feature['db_xrefs']
                    db_xrefs.append('SO:0001217')
                    if(bu.has_annotation(psc, DB_PSC_COL_UNIREF90)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_90, psc[DB_PSC_COL_UNIREF90]))
                    if(bu.has_annotation(psc, DB_PSC_COL_UNIREF50)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_UNIREF_50, psc[DB_PSC_COL_UNIREF50]))
                    if(bu.has_annotation(psc, DB_PSC_COL_COG_ID)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_COG, psc[DB_PSC_COL_COG_ID]))
                    if(bu.has_annotation(psc, DB_PSC_COL_COG_CAT)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_COG, psc[DB_PSC_COL_COG_CAT]))
                    if(bu.has_annotation(psc, DB_PSC_COL_GO)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_GO, psc[DB_PSC_COL_GO]))
                    if(bu.has_annotation(psc, DB_PSC_COL_EC)):
                        db_xrefs.append('%s:%s' % (bc.DB_XREF_EC, psc[DB_PSC_COL_EC]))

                    log.debug(
                        'lookup: contig=%s, start=%i, stop=%i, strand=%s, UniRef90=%s, EC=%s, gene=%s, product=%s',
                        feature['contig'], feature['start'], feature['stop'], feature['strand'], psc.get(DB_PSC_COL_UNIREF90, ''), psc.get(DB_PSC_COL_EC, ''), psc.get(DB_PSC_COL_GENE, ''), psc.get(DB_PSC_COL_PRODUCT, '')
                    )
                else:
                    log.debug('lookup failed! uniref90_id=%s', uniref90_id)

    except Exception as ex:
        log.exception('Could not read PSCs from db!', ex)
        raise Exception("SQL error!", ex)


def parse_psc_annotation(rec):
    psc = {
        DB_PSC_COL_UNIREF90: bc.DB_PREFIX_UNIREF_90 + rec[DB_PSC_COL_UNIREF90],  # must not be NULL/None
    }

    # add non-empty PSC annotations and attach database prefixes to identifiers
    if(rec[DB_PSC_COL_GENE]):
        psc[DB_PSC_COL_GENE] = rec[DB_PSC_COL_GENE]
    if(rec[DB_PSC_COL_PRODUCT]):
        psc[DB_PSC_COL_PRODUCT] = rec[DB_PSC_COL_PRODUCT]
    if(rec[DB_PSC_COL_EC]):
        psc[DB_PSC_COL_EC] = rec[DB_PSC_COL_EC]
    # if(rec[DB_PSC_COL_UNIREF50]):
        # psc[DB_PSC_COL_UNIREF50] = bc.DB_PREFIX_UNIREF_50 + rec[DB_PSC_COL_UNIREF50]
    if(rec[DB_PSC_COL_COG_ID]):
        psc[DB_PSC_COL_COG_ID] = bc.DB_PREFIX_COG + rec[DB_PSC_COL_COG_ID]
    if(rec[DB_PSC_COL_COG_CAT]):
        psc[DB_PSC_COL_COG_CAT] = rec[DB_PSC_COL_COG_CAT]
    if(rec[DB_PSC_COL_GO]):
        go_ids = []
        for go_id in rec[DB_PSC_COL_GO].split(';'):
            if(go_id is not ''):
                go_ids.append(bc.DB_PREFIX_GO + go_id)
        if(len(go_ids) != 0):
            psc[DB_PSC_COL_GO] = go_ids
    
    return psc

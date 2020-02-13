
import logging
import subprocess as sp
import sqlite3

import bakta.config as cfg
import bakta.constants as bc

############################################################################
# PSC DB columns
############################################################################
DB_PSC_COL_UNIREF90 = 'uniref90_id'
DB_PSC_COL_COG = 'cog_id_function'
DB_PSC_COL_EC = 'ec'
DB_PSC_COL_GO = 'go'
DB_PSC_COL_IS = 'is'
DB_PSC_COL_GENE = 'gene'
DB_PSC_COL_PRODUCT = 'product'

log = logging.getLogger('psc')


def search_pscs(cdss):
    """Conduct homology search of CDSs against PCS db."""
    cds_fasta_path = cfg.tmp_path.joinpath('cds.faa')
    with cds_fasta_path.open(mode='w') as fh:
        for cds in cdss:
            fh.write(">%s\n%s\n" % (cds['tmp_id'], cds['sequence']))
    ghostz_output_path = cfg.tmp_path.joinpath('ghostz.tsv')
    ghostz_db_path = cfg.db_path.joinpath('psc')
    cmd = [
        'ghostz',
        'aln',
        '-d', str(ghostz_db_path),
        '-i', str(cds_fasta_path),
        '-o', str(ghostz_output_path),
        '-b', '1',  # single best output
        '-a', str(cfg.threads)
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
            'PSC failed! ghostz-error-code=%d', proc.returncode
        )
        log.debug(
            'PSC: cmd=%s stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
        raise Exception("GhostZ error! error code: %i" % proc.returncode)

    cds_by_id = {cds['tmp_id']: cds for cds in cdss}
    with ghostz_output_path.open() as fh:
        for line in fh:
            (cds_id, cluster_id, identity, alignment_length, align_mismatches,
                align_gaps, query_start, query_stop, subject_start, subject_stop,
                evalue, bitscore) = line.split('\t')
            cds = cds_by_id[cds_id]
            if(int(alignment_length) / len(cds['sequence']) >= bc.MIN_PROTEIN_COVERAGE
                    and float(identity) >= bc.MIN_PROTEIN_IDENTITY):
                cds['psc'] = {
                    'uniref90_id': cluster_id
                }

    pscs_found, pscs_not_found = []
    for cds in cdss:
        if('psc' in cds):
            pscs_found.append(cds)
        else:
            cds['hypothetical'] = True
            pscs_not_found.append(cds)
    log.info('PSCs: # %i', len(pscs_found))
    return pscs_found, pscs_not_found


def lookup_pscs(cdss):
    """Lookup PCS information"""
    try:
        with sqlite3.connect("file:%s?mode=ro" % str(cfg.db_path.joinpath('psc.db')), uri=True) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for cds in cdss:
                if('psc' in cds):
                    uniref90_id = cds['psc']['uniref90_id']
                elif('ups' in cds):
                    uniref90_id = cds['ups']['uniref90_id']
                else:
                    continue  # skip this cds object
                rec = c.fetchone("select * from psc where uniref90_id=?", (uniref90_id,))
                if(rec is not None):
                    psc = {
                        'uniref90_id': bc.DB_PREFIX_UNIREF_90 + rec[DB_PSC_COL_UNIREF90],  # must not be NULL/None
                        'ec_id': rec[DB_PSC_COL_EC],
                        'is_id': rec[DB_PSC_COL_IS],
                        'gene': rec[DB_PSC_COL_GENE],
                        'product': rec[DB_PSC_COL_PRODUCT]
                    }

                    # reattach database prefixes to identifiers
                    if(rec[DB_PSC_COL_COG] is not None):
                        cog_id, cog_category = rec[DB_PSC_COL_COG].split(';')
                        cog_id = bc.DB_PREFIX_COG + cog_id
                        psc['cog_id'] = cog_id
                        psc['cog_category'] = cog_category
                    if(rec[DB_PSC_COL_GO] is not None):
                        go_ids = []
                        for go_id in rec[DB_PSC_COL_GO].split(';'):
                            go_id.append("%s:%s" % (bc.DB_PREFIX_GO, go_id))
                        psc['go_ids'] = go_ids

                    cds['psc'] = psc
                    if('db_xrefs' not in cds):
                        cds['db_xrefs'] = []
                    db_xrefs = cds['db_xrefs']
                    db_xrefs.append('SO:0001217')

                    log.debug(
                        'PSC: contig=%s, start=%i, stop=%i, strand=%s, UniRef90=%s',
                        cds['contig'], cds['gene'], cds['start'], cds['stop'], cds['strand'], psc['uniref90_id']
                    )

    except Exception:
        log.exception('Could not read UPSs from db!')

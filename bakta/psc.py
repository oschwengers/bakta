
import logging
import subprocess as sp
import sqlite3

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


def search_pscs(config, cdss):
    """Conduct homology search of CDSs against PCS db."""
    cds_fasta_path = config['tmp'].joinpath('cds.faa')
    with cds_fasta_path.open(mode='w') as fh:
        for cds in cdss:
            fh.write(">%s\n%s\n" % (cds['tmp_id'], cds['sequence']))
    ghostz_output_path = config['tmp'].joinpath('ghostz.tsv')
    ghostz_db_path = config['db'].joinpath('psc')
    cmd = [
        'ghostz',
        'aln',
        '-d', str(ghostz_db_path),
        '-i', str(cds_fasta_path),
        '-o', str(ghostz_output_path),
        '-b', '1',  # single best output
        '-a', str(config['threads'])
    ]
    proc = sp.run(
        cmd,
        cwd=str(config['tmp']),
        env=config['env'],
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
    for cds in cdss:
        if('psc' not in cds):
            cds['hypothetical'] = True


def lookup_pscs(config, cdss):
    """Lookup PCS information"""
    try:
        seqs_found = []
        seqs_not_found = []
        with sqlite3.connect("file:%s?mode=ro" % str(config['db'].joinpath('psc.db')), uri=True) as conn:
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
                        'uniref90_id': bc.DB_PREFIX_UNIREF_90 + rec[DB_PSC_COL_UNIREF90],
                        'cog_id': bc.DB_PREFIX_COG + rec[DB_PSC_COL_COG].split(';')[0],
                        'cog_function': rec[DB_PSC_COL_COG],
                        'ec': rec[DB_PSC_COL_EC],
                        'is_id': rec[DB_PSC_COL_IS],
                        'gene': rec[DB_PSC_COL_GENE],
                        'product': rec[DB_PSC_COL_PRODUCT],
                        'go': rec[DB_PSC_COL_GO]
                    }
                    cds['psc'] = psc
                    if('db_xrefs' not in cds):
                        cds['db_xrefs'] = []
                    seqs_found.append(cds)

                    log.debug(
                        'PSC: contig=%s, start=%i, stop=%i, strand=%s, UniRef90=%s',
                        cds['contig'], cds['gene'], cds['start'], cds['stop'], cds['strand'], psc['uniref90_id']
                    )
                else:
                    seqs_not_found.append(cds)

        log.info('UPSs: # %i', len(seqs_found))
        return seqs_found, seqs_not_found
    except Exception:
        log.exception('Could not read UPSs from db!')

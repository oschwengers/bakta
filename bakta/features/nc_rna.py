
import logging
import subprocess as sp

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:nc_rna')


def predict_nc_rnas(data, contigs_path):
    """Search for non-coding RNA genes."""

    output_path = cfg.tmp_path.joinpath('ncrna-genes.tsv')
    cmd = [
        'cmsearch',
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--rfam',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path),
        str(cfg.db_path.joinpath('ncRNA-genes')),
        str(contigs_path)
    ]
    if(data['genome_size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(data['genome_size'] // 1000000))
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
        log.warning('ncRNAs failed! cmscan-error-code=%d', proc.returncode)
        raise Exception("cmsearch error! error code: %i" % proc.returncode)

    rfam2go = {}
    rfam2go_path = cfg.db_path.joinpath('rfam-go.tsv')
    with rfam2go_path.open() as fh:
        for line in fh:
            (rfam, go) = line.split('\t')
            if(rfam in rfam2go):
                rfam2go[rfam].append(go)
            else:
                rfam2go[rfam] = [go]

    ncrnas = []
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, accession, subject, subject_id, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = line.strip().split()
                
                if(strand == '-'):
                    (start, stop) = (stop, start)
                
                rfam_id = "RFAM:%s" % subject_id
                db_xrefs = [rfam_id, 'SO:0001263']
                if(rfam_id in rfam2go):
                    db_xrefs += rfam2go[rfam_id]
                ncrna = {
                    'type': bc.FEATURE_NC_RNA,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'gene': subject,
                    'product': description,
                    'score': float(score),
                    'evalue': float(evalue),
                    'db_xrefs': db_xrefs
                }
                ncrnas.append(ncrna)
                log.debug(
                    'contig=%s, start=%i, stop=%i, strand=%s, gene=%s',
                    ncrna['contig'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['gene']
                )
    log.info('# %i', len(ncrnas))
    return ncrnas

import logging
import subprocess as sp

import bakta.config as cfg
import bakta.constants as bc

log = logging.getLogger('features:nc_rna_regions')


def predict_nc_rna_regions(data, contigs_path):
    """Search for non-coding RNA regions."""

    output_path = cfg.tmp_path.joinpath('ncrna-regions.tsv')
    cmd = [
        'cmsearch',
        '--noali',
        '--cut_tc',
        '--notrunc',
        '--rfam',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path),
        str(cfg.db_path.joinpath('ncRNA-regions')),
        str(contigs_path)
    ]
    if(data['genome_size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(data['genome_size'] // 1000000))
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.warning('ncRNA regions failed! cmscan-error-code=%d', proc.returncode)
        log.debug(
            'ncRNA regions: cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            cmd, proc.stdout, proc.stderr
        )
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
                    'type': bc.FEATURE_NC_RNA_REGION,
                    'contig': contig,
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'product': description,
                    'score': float(score),
                    'evalue': float(evalue),
                    'db_xrefs': db_xrefs
                }
                ncrnas.append(ncrna)
                log.debug(
                    'ncRNA regions: contig=%s, start=%i, stop=%i, strand=%s, product=%s',
                    ncrna['contig'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['product']
                )
    log.info('ncRNA regions: # %i', len(ncrnas))
    return ncrnas
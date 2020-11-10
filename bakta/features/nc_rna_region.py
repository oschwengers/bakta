
import logging
import re
import subprocess as sp
from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so

log = logging.getLogger('NC_RNA_REGION')


def predict_nc_rna_regions(genome, contigs_path):
    """Search for non-coding RNA regions."""

    output_path = cfg.tmp_path.joinpath('ncrna-regions.tsv')
    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '-g',  # activate glocal mode
        '--nohmmonly',  # strictly use CM models
        '--rfam',
        '--cpu', str(cfg.threads),
        '--tblout', str(output_path)
    ]
    if(genome['size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(2 * genome['size'] // 1000000))
    cmd.append(str(cfg.db_path.joinpath('ncRNA-regions')))
    cmd.append(str(contigs_path))
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
        log.warning('ncRNA regions failed! cmscan-error-code=%d', proc.returncode)
        raise Exception(f'cmscan error! error code: {proc.returncode}')

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
                (subject, accession, contig_id, contig_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = re.split('\s+', line.strip(), maxsplit=17)
                
                if(strand == '-'):
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                length = stop - start + 1
                partial = trunc != 'no'
                
                if(evalue > 1E-6):
                    log.debug(
                        'discard low E value: contig=%s, start=%i, stop=%i, strand=%s, gene=%s, partial=%s, length=%i, evalue=%f',
                        contig_id, start, stop, strand, subject, partial, length, evalue
                    )
                else:
                    rfam_id = f'RFAM:{accession}'
                    db_xrefs = [rfam_id]
                    if(rfam_id in rfam2go):
                        db_xrefs += rfam2go[rfam_id]
                    
                    ncrna_region = OrderedDict()
                    ncrna_region['type'] = bc.FEATURE_NC_RNA_REGION
                    ncrna_region['contig'] = contig_id
                    ncrna_region['start'] = start
                    ncrna_region['stop'] = stop
                    ncrna_region['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    ncrna_region['gene'] = subject
                    ncrna_region['product'] = f'(partial) {description}' if partial else description
                    
                    if(partial):
                        ncrna_region['partial'] = partial
                    
                    ncrna_region['score'] = float(score)
                    ncrna_region['evalue'] = evalue
                    
                    if('5' in trunc):
                        ncrna_region['trunc_5'] = True
                    if('3' in trunc):
                        ncrna_region['trunc_3'] = True

                    ncrna_region['db_xrefs'] = db_xrefs

                    ncrnas.append(ncrna_region)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, product=%s, length=%i, evalue=%f',
                        ncrna_region['contig'], ncrna_region['start'], ncrna_region['stop'], ncrna_region['strand'], ncrna_region['product'], length, ncrna_region['evalue']
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas
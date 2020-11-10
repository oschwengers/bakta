
import logging
import re
import subprocess as sp
from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so

log = logging.getLogger('NC_RNA')


def predict_nc_rnas(genome, contigs_path):
    """Search for non-coding RNA genes."""

    output_path = cfg.tmp_path.joinpath('ncrna-genes.tsv')
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
    cmd.append(str(cfg.db_path.joinpath('ncRNA-genes')))
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
        log.warning('ncRNAs failed! cmscan-error-code=%d', proc.returncode)
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
                    db_xrefs = [rfam_id, so.SO_NCRNA_GENE.id]
                    if(rfam_id in rfam2go):
                        db_xrefs += rfam2go[rfam_id]
                    
                    ncrna = OrderedDict()
                    ncrna['type'] = bc.FEATURE_NC_RNA
                    ncrna['contig'] = contig_id
                    ncrna['start'] = start
                    ncrna['stop'] = stop
                    ncrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    ncrna['gene'] = subject
                    ncrna['product'] = f'(partial) {description}' if partial else description
                    
                    if(partial):
                        ncrna['partial'] = partial
                    
                    ncrna['score'] = float(score)
                    ncrna['evalue'] = evalue
                    
                    if('5' in trunc):
                        ncrna['trunc_5'] = True
                    if('3' in trunc):
                        ncrna['trunc_3'] = True

                    ncrna['db_xrefs'] = db_xrefs

                    ncrnas.append(ncrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, product=%s, length=%i, evalue=%f',
                        ncrna['contig'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['product'], length, ncrna['evalue']
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas
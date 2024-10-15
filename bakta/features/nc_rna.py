import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.features.annotation as ba
import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_EVALUE = 1E-4


log = logging.getLogger('NC_RNA')


def predict_nc_rnas(data: dict, sequences_path: Path):
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
    if(data['stats']['size'] >= 1000000):
        cmd.append('-Z')
        cmd.append(str(2 * data['stats']['size'] // 1000000))
    cmd.append(str(cfg.db_path.joinpath('ncRNA-genes')))
    cmd.append(str(sequences_path))
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
    sequences = {seq['id']: seq for seq in data['sequences']}
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (
                    subject, accession, sequence_id, sequence_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description
                ) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

                if(strand == '-'):
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                score = float(score)
                length = stop - start + 1
                if(trunc == "5'"):
                    truncated = bc.FEATURE_END_5_PRIME
                elif(trunc == "3'"):
                    truncated = bc.FEATURE_END_3_PRIME
                else:
                    truncated = None

                if(evalue > HIT_EVALUE):
                    log.debug(
                        'discard low E value: seq=%s, start=%i, stop=%i, strand=%s, gene=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        sequence_id, start, stop, strand, subject, length, truncated, score, evalue
                    )
                else:
                    rfam_id = f'{bc.DB_XREF_RFAM}:{accession}'
                    db_xrefs = [rfam_id]
                    if(rfam_id in rfam2go):
                        db_xrefs += rfam2go[rfam_id]

                    ncrna = OrderedDict()
                    ncrna['type'] = bc.FEATURE_NC_RNA
                    ncrna['class'] = determine_class(description)
                    ncrna['sequence'] = sequence_id
                    ncrna['start'] = start
                    ncrna['stop'] = stop
                    ncrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    
                    gene = subject
                    if(ba.RE_PROTEIN_SYMBOL.fullmatch(gene)):
                        gene = gene[0].lower() + gene[1:]
                        log.debug('fix gene: lowercase first char. new=%s, old=%s', gene, subject)
                    ncrna['gene'] = gene
                    ncrna['product'] = description

                    if(ncrna['class'] is not None):
                        db_xrefs.append(ncrna['class'].id)
                    else:
                        db_xrefs.append(so.SO_NCRNA_GENE.id)

                    if(truncated):
                        ncrna['truncated'] = truncated

                    ncrna['score'] = score
                    ncrna['evalue'] = evalue
                    ncrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(ncrna, sequences[sequence_id])  # extract nt sequences
                    ncrna['nt'] = nt

                    ncrnas.append(ncrna)
                    log.info(
                        'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        ncrna['sequence'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['gene'], ncrna['product'], length, truncated, ncrna['score'], ncrna['evalue'], nt[:10], nt[-10:]
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas


def determine_class(description: str) -> str:
    description = description.lower()
    if('ribozyme' in description):
        return so.SO_NCRNA_GENE_RIBOZYME
    elif('rnase p' in description):
        return so.SO_NCRNA_GENE_RNASEP
    elif('antisense' in description):
        return so.SO_NCRNA_GENE_ANTISENSE
    else:
        None

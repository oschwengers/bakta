import logging
import subprocess as sp

from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


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
    contigs = {c['id']: c for c in genome['contigs']}
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (
                    subject, accession, contig_id, contig_acc, mdl, mdl_from, mdl_to,
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

                if(evalue > 1E-4):
                    log.debug(
                        'discard low E value: contig=%s, start=%i, stop=%i, strand=%s, gene=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        contig_id, start, stop, strand, subject, length, truncated, score, evalue
                    )
                else:
                    rfam_id = f'{bc.DB_XREF_RFAM}:{accession}'
                    db_xrefs = [rfam_id]
                    if(rfam_id in rfam2go):
                        db_xrefs += rfam2go[rfam_id]

                    ncrna = OrderedDict()
                    ncrna['type'] = bc.FEATURE_NC_RNA
                    ncrna['class'] = determine_class(description)
                    ncrna['contig'] = contig_id
                    ncrna['start'] = start
                    ncrna['stop'] = stop
                    ncrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    ncrna['gene'] = subject

                    if(truncated is None):
                        ncrna['product'] = description
                    elif(truncated == bc.FEATURE_END_UNKNOWN):
                        ncrna['product'] = f'(partial) {description}'
                    elif(truncated == bc.FEATURE_END_5_PRIME):
                        ncrna['product'] = f"(5' truncated) {description}"
                    elif(truncated == bc.FEATURE_END_3_PRIME):
                        ncrna['product'] = f"(3' truncated) {description}"

                    if(ncrna['class'] is not None):
                        db_xrefs.append(ncrna['class'].id)
                    else:
                        db_xrefs.append(so.SO_NCRNA_GENE.id)

                    if(truncated):
                        ncrna['truncated'] = truncated

                    ncrna['score'] = score
                    ncrna['evalue'] = evalue
                    ncrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(ncrna, contigs[contig_id])  # extract nt sequences
                    ncrna['nt'] = nt

                    ncrnas.append(ncrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        ncrna['contig'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['gene'], ncrna['product'], length, truncated, ncrna['score'], ncrna['evalue'], nt[:10], nt[-10:]
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas


def determine_class(description):
    description = description.lower()
    if('ribozyme' in description):
        return so.SO_NCRNA_GENE_RIBOZYME
    elif('rnase p' in description):
        return so.SO_NCRNA_GENE_RNASEP
    elif('antisense' in description):
        return so.SO_NCRNA_GENE_ANTISENSE
    else:
        None

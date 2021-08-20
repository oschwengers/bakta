import logging
import subprocess as sp

from collections import OrderedDict

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('R_RNA')


def predict_r_rnas(genome, contigs_path):
    """Search for ribosomal RNA sequences."""

    output_path = cfg.tmp_path.joinpath('rrna.tsv')
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
    cmd.append(str(cfg.db_path.joinpath('rRNA')))
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
        log.warning('rRNAs failed! cmscan-error-code=%d', proc.returncode)
        raise Exception(f'cmscan error! error code: {proc.returncode}')

    rrnas = []
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

                db_xrefs = ['GO:0005840', 'GO:0003735']
                if(accession == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += ['RFAM:RF00001', so.SO_RRNA_5S.id]
                    consensus_length = 119
                elif(accession == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += ['RFAM:RF00177', so.SO_RRNA_16S.id]
                    consensus_length = 1533
                elif(accession == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += ['RFAM:RF02541', so.SO_RRNA_23S.id]
                    consensus_length = 2925
                else:
                    log.warning(
                        'unknown rRNA detected! accession=%s, contig=%s, start=%i, stop=%i, strand=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        accession, contig_id, start, stop, strand, length, truncated, score, evalue
                    )
                    continue

                coverage = length / consensus_length
                if(coverage < 0.8):
                    truncated = bc.FEATURE_END_UNKNOWN

                if(coverage < 0.3):
                    log.debug(
                        'discard low coverage: contig=%s, rRNA=%s, start=%i, stop=%i, strand=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e',
                        contig_id, rrna_tag, start, stop, strand, length, coverage, truncated, score, evalue
                    )
                else:
                    rrna = OrderedDict()
                    rrna['type'] = bc.FEATURE_R_RNA
                    rrna['contig'] = contig_id
                    rrna['start'] = start
                    rrna['stop'] = stop
                    rrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    rrna['gene'] = f'{rrna_tag}_rrna'

                    if(truncated is None):
                        rrna['product'] = f'{rrna_tag} ribosomal RNA'
                    elif(truncated == bc.FEATURE_END_UNKNOWN):
                        rrna['product'] = f'(partial) {rrna_tag} ribosomal RNA'
                    elif(truncated == bc.FEATURE_END_5_PRIME):
                        rrna['product'] = f"(5' truncated) {rrna_tag} ribosomal RNA"
                    elif(truncated == bc.FEATURE_END_3_PRIME):
                        rrna['product'] = f"(3' truncated) {rrna_tag} ribosomal RNA"

                    if(truncated):
                        rrna['truncated'] = truncated

                    rrna['coverage'] = coverage
                    rrna['score'] = score
                    rrna['evalue'] = evalue
                    rrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(rrna, contigs[contig_id])  # extract nt sequences
                    rrna['nt'] = nt

                    rrnas.append(rrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        rrna['contig'], rrna['start'], rrna['stop'], rrna['strand'], rrna['gene'], rrna['product'], length, coverage, truncated, score, evalue, nt[:10], nt[-10:]
                    )

    log.info('predicted=%i', len(rrnas))
    return rrnas


import logging
import subprocess as sp
from collections import OrderedDict

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu
import bakta.so as so

log = logging.getLogger('CDS')


def predict(genome, filtered_contigs_path):
    """Predict open reading frames with Prodigal."""

    proteins_path = cfg.tmp_path.joinpath('proteins.faa')
    gff_path = cfg.tmp_path.joinpath('prodigal.gff')
    cmd = [
        'prodigal',
        '-i', str(filtered_contigs_path),
        '-a', str(proteins_path),
        '-g', str(cfg.translation_table),  # set translation table
        '-f', 'gff',  # GFF output
        '-o', str(gff_path)  # prodigal output
    ]
    if(genome['complete'] == False):
        cmd.append('-c')  # closed ends
    if(cfg.prodigal_tf):
        cmd.append('-t')  # use supplied prodigal training file
        cmd.append(str(cfg.prodigal_tf))
    elif(genome['size'] < 20000):  # not enough sequence information and no trainings file provided
        cmd.append('-p')  # run prodigal in meta mode
        cmd.append('meta')
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
        log.warning('ORFs failed! prodigal-error-code=%d', proc.returncode)
        raise Exception(f'prodigal error! error code: {proc.returncode}')

    # parse orfs
    # TODO: replace code by BioPython GFF3 parser
    contigs = {k['id']: k for k in genome['contigs']}
    cdss = {}
    partial_cdss = {}
    partial_cdss_per_contig = {k['id']: [] for k in genome['contigs']}
    with gff_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (contig, inference, _, start, stop, score, strand, _, annotations_raw) = line.strip().split('\t')
                gff_annotations = split_gff_annotation(annotations_raw)
                contig_orf_id = gff_annotations['ID'].split('_')[1]

                cds = OrderedDict()
                cds['type'] = bc.FEATURE_CDS
                cds['contig'] = contig
                cds['start'] = int(start)
                cds['stop'] = int(stop)
                cds['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                cds['gene'] = None
                cds['product'] = None
                cds['start_type'] = gff_annotations['start_type']
                cds['rbs_motif'] = gff_annotations['rbs_motif']
                cds['db_xrefs'] = [so.SO_CDS.id]
                
                if(cds['strand'] == bc.STRAND_FORWARD):
                    cds['frame'] = (cds['start'] - 1) % 3 + 1
                else:
                    cds['frame'] = (contigs[cds['contig']]['length'] - cds['stop']) % 3 + 1
                
                if(gff_annotations['partial'] == '10'):
                    cds['truncated'] = '5-prime' if cds['strand'] == bc.STRAND_FORWARD else '3-prime'
                    partial_cdss[f"{cds['contig']}_{contig_orf_id}"] = cds
                    partial_cdss_per_contig[cds['contig']].append(cds)
                elif(gff_annotations['partial'] == '01'):
                    cds['truncated'] = '3-prime' if cds['strand'] == bc.STRAND_FORWARD else '5-prime'
                    partial_cdss[f"{cds['contig']}_{contig_orf_id}"] = cds
                    partial_cdss_per_contig[cds['contig']].append(cds)
                else:
                    cdss[f"{cds['contig']}_{contig_orf_id}"] = cds
                
                log.info(
                    'contig=%s, start=%i, stop=%i, strand=%s, frame=%s, truncated=%s, start-type=%s, RBS-motif=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds.get('truncated', '-'), cds['start_type'], cds['rbs_motif']
                )

    # extract translated orf sequences
    with proteins_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            cds = cdss.get(record.id, None)
            if(cds):
                seq = str(record.seq)[:-1]  # discard trailing asterisk
                cds['sequence'] = seq
                cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(seq)
            else:
                partial_cds = partial_cdss.get(record.id, None)
                if(partial_cds):
                    seq = str(record.seq)
                    partial_cds['sequence'] = seq
                    log.debug(
                        'store trunc CDS: contig=%s, start=%i, stop=%i, strand=%s, trunc=%s, seq=%s',
                        partial_cds['contig'], partial_cds['start'], partial_cds['stop'], partial_cds['strand'], partial_cds.get('truncated', '-'), seq
                    )
    cdss = list(cdss.values())

    # merge matching partial features on sequence edges
    for partial_cdss_pairs in partial_cdss_per_contig.values():
        if(len(partial_cdss_pairs) >= 2):  # skip unpaired edge features
            first_partial_cds = partial_cdss_pairs[0]  # first partial CDS per contig
            last_partial_cds = partial_cdss_pairs[-1]  # last partial CDS per contig
            # check if partial CDS are on same strand
            # and have opposite truncated edges
            # and firtst starts at 1
            # and last ends at contig end (length)
            if(first_partial_cds['strand'] == last_partial_cds['strand']
                and first_partial_cds['truncated'] != last_partial_cds['truncated']
                and first_partial_cds['start'] == 1
                and last_partial_cds['stop'] == contigs[last_partial_cds['contig']]['length']):
                cds = last_partial_cds
                cds['stop'] = first_partial_cds['stop']
                if(last_partial_cds['truncated'] == '3-prime'):
                    seq = last_partial_cds['sequence'] + first_partial_cds['sequence']  # merge sequence
                else:
                    seq = first_partial_cds['sequence'] + last_partial_cds['sequence']  # merge sequence
                    cds['start_type'] = first_partial_cds['start_type']
                    cds['rbs_motif'] = first_partial_cds['rbs_motif']
                log.debug(f'trunc seq: seq-start={seq[:10]}, seq-end={seq[-10:]}')

                cds['edge'] = True  # mark CDS as edge feature
                cds.pop('truncated')

                seq = seq[:-1]  # discard trailing asterisk
                cds['sequence'] = seq
                cds['aa_digest'], cds['aa_hexdigest'] = bu.calc_aa_hash(seq)
                cdss.append(cds)
                log.info(
                    'edge CDS: contig=%s, start=%i, stop=%i, strand=%s, frame=%s, start-type=%s, RBS-motif=%s, aa-hexdigest=%s, seq=[%s..%s]',
                    cds['contig'], cds['start'], cds['stop'], cds['strand'], cds['frame'], cds['start_type'], cds['rbs_motif'], cds['aa_hexdigest'], seq[:10], seq[-10:]
                )

    log.info('predicted=%i', len(cdss))
    return cdss


def split_gff_annotation(annotation_string):
    annotations = {}
    for expr in annotation_string.split(';'):
        if(expr != ''):
            try:
                key, value = expr.split('=')
                annotations[key] = value
            except:
                log.error('expr=%s', expr)
    return annotations


def mark_hypotheticals(cdss):
    for cds in cdss:
        if('ips' not in cds and 'psc' not in cds):
            cds['hypothetical'] = True

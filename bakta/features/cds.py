
import logging
import subprocess as sp
from collections import OrderedDict

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc
import bakta.utils as bu
import bakta.so as so

log = logging.getLogger('features:cds')


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
    if(cfg.complete == False):
        cmd.append('-c')  # closed ends
    if(cfg.prodigal_tf):
        cmd.append('-t')  # use supplied prodigal training file
        cmd.append(str(cfg.prodigal_tf))
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
        raise Exception("prodigal error! error code: %i" % proc.returncode)

    # parse orfs
    # TODO: replace code by BioPython GFF3 parser
    contigs = {k['id']: k for k in genome['contigs']}
    cdss = {}
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
                cdss["%s_%s" % (cds['contig'], contig_orf_id)] = cds
                log.debug(
                    'contig=%s, start=%i, stop=%i, strand=%s',
                    cds['contig'], cds['start'], cds['stop'], cds['strand']
                )

    # extract translated orf sequences
    with proteins_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            cds = cdss[record.id]
            seq = str(record.seq)[:-1]  # discard trailing asterisk
            cds['sequence'] = seq
            (aa_digest, aa_hexdigest) = bu.calc_aa_hash(seq)
            cds['aa_digest'] = aa_digest
            cds['aa_hexdigest'] = aa_hexdigest

    gff_path.unlink()
    proteins_path.unlink()
    log.info('# %i', len(cdss))
    return list(cdss.values())


def split_gff_annotation(annotation_string):
    annotations = {}
    for expr in annotation_string.split(';'):
        if(expr != ''):
            try:
                key, value = expr.split('=')
                annotations[key] = value
            except:
                log.error('expr=%s' % expr)
    return annotations


def mark_hypotheticals(cdss):
    for cds in cdss:
        if('ips' not in cds and 'psc' not in cds):
            cds['hypothetical'] = True


import logging

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('io:fasta')

FASTA_LINE_WRAPPING = 60

def import_contigs(contigs_path, min_length):
    """Apply min-length filters and rename contig headers."""

    contig_counter = 1
    contigs = []
    discarded = []
    contig_prefix = cfg.locus if cfg.locus else 'contig'
    with contigs_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq)
            contig_name = "%s_%i" % (contig_prefix, contig_counter)
            contig_counter += 1
            contig = {
                'id': contig_name,
                'org_id': record.id,
                'org_desc': record.description,
                'sequence': seq,
                'length': len(seq)
            }
            if len(seq) >= min_length:
                contigs.append(contig)
            else:
                discarded.append(contig)

    return contigs, discarded


def export_contigs(contigs, fasta_path):
    """Write valid contigs to Fasta file."""
    with fasta_path.open('w') as fh:
        for c in contigs:
            fh.write(">%s\n%s\n" % (c['id'], c['sequence']))


def format_fasta(contig, line_wrapping=False):
    lines = '>%s\n' % contig['id']
    if line_wrapping:
        seq = contig['sequence']
        while len(seq) > FASTA_LINE_WRAPPING:
            lines += seq[:FASTA_LINE_WRAPPING] + '\n'
            seq = seq[FASTA_LINE_WRAPPING:]
        lines += seq + '\n'
    else:
        lines += contig['sequence'] + '\n'
    return lines


def write_faa(features, faa_path):
    """Write translated CDS sequences to Fasta file."""
    with faa_path.open('w') as fh:
        for feat in features:
            if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):
                fh.write(">%s %s\n%s\n" % (feat['locus'], feat['product'], feat['sequence']))

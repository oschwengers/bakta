
import logging

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('io:fasta')

FASTA_LINE_WRAPPING = 60

def import_contigs(contigs_path):
    """Import raw contigs."""
    contigs = []
    with contigs_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            contig = {
                'id': record.id,
                'desc': record.description,
                'sequence': seq,
                'length': len(seq),
                'complete': False,
                'type': bc.REPLICON_CONTIG
            }
            contigs.append(contig)
    return contigs


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

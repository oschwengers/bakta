
import logging

from Bio import SeqIO


log = logging.getLogger('io:fasta')


def import_contigs(contigs_path, min_length):
    """Apply min-length filters and rename contig headers."""

    contig_counter = 1
    contigs = []
    discarded = []
    with contigs_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq)
            contig_name = "contig_%d" % contig_counter
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



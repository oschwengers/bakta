
import logging

from Bio import SeqIO
from xopen import xopen

import bakta.config as cfg
import bakta.constants as bc


log = logging.getLogger('FASTA')

FASTA_LINE_WRAPPING = 60

def import_contigs(contigs_path):
    """Import raw contigs."""
    contigs = []
    # with contigs_path.open() as fh:
    with xopen(str(contigs_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            contig = {
                'id': record.id,
                'description': record.description,
                'sequence': seq,
                'length': len(seq),
                'complete': False,
                'type': bc.REPLICON_CONTIG,
                'topology': bc.TOPOLOGY_LINEAR
            }
            log.info(
                'imported: id=%s, length=%i, complete=%s, topology=%s, description=%s', 
                contig['id'], contig['length'], contig['complete'], contig['topology'], contig['description']
            )
            contigs.append(contig)
    return contigs


def export_contigs(contigs, fasta_path, description=False, wrap=False):
    """Write contigs to Fasta file."""
    with fasta_path.open('w') as fh:
        for contig in contigs:
            if(description):
                fh.write(f">{contig['id']} {contig['description']}\n")
            else:
                fh.write(f">{contig['id']}\n")
            if(wrap):
                fh.write(wrap_sequence(contig['sequence']))
            else:
                fh.write(contig['sequence'])
                fh.write('\n')


def wrap_sequence(sequence):
    lines = ''
    while len(sequence) > FASTA_LINE_WRAPPING:
        lines += sequence[:FASTA_LINE_WRAPPING] + '\n'
        sequence = sequence[FASTA_LINE_WRAPPING:]
    lines += sequence + '\n'
    return lines


def write_faa(features, faa_path):
    """Write translated CDS sequences to Fasta file."""
    with faa_path.open('w') as fh:
        for feat in features:
            if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):
                fh.write(f">{feat['locus']} {feat['product']}\n{feat['sequence']}\n")

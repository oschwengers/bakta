import logging
import re

from Bio import SeqIO
from xopen import xopen

import bakta.constants as bc


log = logging.getLogger('FASTA')


FASTA_DNA_SEQUENCE_PATTERN = re.compile(r'[ATGCNMRWSYKVHDBN]+', re.IGNORECASE)
FASTA_LINE_WRAPPING = 60


def import_contigs(contigs_path):
    """Import raw contigs."""
    contigs = []
    # with contigs_path.open() as fh:
    with xopen(str(contigs_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            if(FASTA_DNA_SEQUENCE_PATTERN.fullmatch(seq) is None):
                log.error('import: Fasta sequence contains invalid DNA characters! id=%s')
                raise ValueError(f'Fasta sequence contains invalid DNA characters! id={record.id}')
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
    log.info('write genome sequences: path=%s, description=%s, wrap=%s', fasta_path, description, wrap)

    with fasta_path.open('wt') as fh:
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
    lines = []
    for i in range(0, len(sequence), FASTA_LINE_WRAPPING):
        lines.append(sequence[i:i + FASTA_LINE_WRAPPING])
    return '\n'.join(lines) + '\n'


def write_faa(features, faa_path):
    """Write translated CDS sequences to Fasta file."""
    log.info('write translated CDS/sORF: path=%s', faa_path)

    with faa_path.open('wt') as fh:
        for feat in features:
            if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):
                fh.write(f">{feat['locus']} {feat['product']}\n{feat['aa']}\n")


def write_ffn(features, ffn_path):
    """Write translated CDS sequences to Fasta file."""
    log.info('write feature nucleotide sequences: path=%s', ffn_path)

    with ffn_path.open('wt') as fh:
        for feat in features:
            if('locus' in feat):
                if(feat.get('product', '') != ''):
                    fh.write(f">{feat['locus']} {feat['product']}\n{feat['nt']}\n")
                else:
                    fh.write(f">{feat['locus']}\n{feat['nt']}\n")

import logging
import re

from pathlib import Path
from typing import Sequence

from Bio import SeqIO
from xopen import xopen

import bakta.constants as bc


log = logging.getLogger('FASTA')


FASTA_AA_SEQUENCE_PATTERN = re.compile(r'[ARNDCQEGHILKMFPOSUTWYVBZXJ]+', re.IGNORECASE)
FASTA_DNA_SEQUENCE_PATTERN = re.compile(r'[ATGCNMRWSYKVHDBN]+', re.IGNORECASE)
FASTA_LINE_WRAPPING = 60


def import_contigs(contigs_path: Path, is_genomic: bool=True, is_dna: bool=True) -> Sequence[dict]:
    """Import raw contigs."""
    contigs = []
    with xopen(str(contigs_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            if('-' in seq):
                dash_count = seq.count('-')
                seq = seq.replace('-', '')
                log.info('import: Discarded alignment gaps (dashes): id=%s, occurences=%i', record.id, dash_count)
            if(is_dna):
                if(FASTA_DNA_SEQUENCE_PATTERN.fullmatch(seq) is None):
                    log.error('import: Fasta sequence contains invalid DNA characters! id=%s', record.id)
                    raise ValueError(f'Fasta sequence contains invalid DNA characters! id={record.id}')
            else:
                if(seq[-1] == '*'):  # remove trailing stop asterik
                    seq = seq[:-1]
                    log.debug('import: Removed trailing asterik! id=%s, seq=%s', record.id, seq)
                if(FASTA_AA_SEQUENCE_PATTERN.fullmatch(seq) is None):
                    log.error('import: Fasta sequence contains invalid AA characters! id=%s, seq=%s', record.id, seq)
                    raise ValueError(f'Fasta sequence contains invalid AA characters! id={record.id}')

            contig = {
                'id': record.id,
                'description': record.description.split(' ', maxsplit=1)[1] if ' ' in record.description else '',
                'sequence': seq,
                'length': len(seq)
            }
            if(is_genomic):
                contig['complete'] = False
                contig['type'] = bc.REPLICON_CONTIG
                contig['topology'] = bc.TOPOLOGY_LINEAR
            log.info(
                'imported: id=%s, length=%i, description=%s, genomic=%s, dna=%s',
                contig['id'], contig['length'], contig['description'], is_genomic, is_dna
            )
            contigs.append(contig)
    return contigs


def export_contigs(contigs: Sequence[dict], fasta_path: Path, description: bool=False, wrap: bool=False):
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


def wrap_sequence(sequence: str):
    lines = []
    for i in range(0, len(sequence), FASTA_LINE_WRAPPING):
        lines.append(sequence[i:i + FASTA_LINE_WRAPPING])
    return '\n'.join(lines) + '\n'


def write_faa(features: Sequence[dict], faa_path: Path):
    """Write translated CDS sequences to Fasta file."""
    log.info('write translated CDS/sORF: path=%s', faa_path)
    with faa_path.open('wt') as fh:
        for feat in features:
            if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):
                fh.write(f">{feat['locus']} {feat['product']}\n{feat['aa']}\n")


def write_ffn(features: Sequence[dict], ffn_path: Path):
    """Write translated CDS sequences to Fasta file."""
    log.info('write feature nucleotide sequences: path=%s', ffn_path)
    with ffn_path.open('wt') as fh:
        for feat in features:
            if('locus' in feat):
                if(feat.get('product', '') != ''):
                    fh.write(f">{feat['locus']} {feat['product']}\n{feat['nt']}\n")
                else:
                    fh.write(f">{feat['locus']}\n{feat['nt']}\n")

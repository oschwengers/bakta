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


def import_sequences(sequences_path: Path, is_genomic: bool=True, is_dna: bool=True) -> Sequence[dict]:
    """Import raw sequences from Fasta file."""
    sequences = []
    with xopen(str(sequences_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            sequence = {
                'id': record.id,
                'description': record.description.split(' ', maxsplit=1)[1] if ' ' in record.description else ''
            }
            
            raw_sequence = str(record.seq).upper()
            if('-' in raw_sequence):
                dash_count = raw_sequence.count('-')
                raw_sequence = raw_sequence.replace('-', '')
                log.info('import: Discarded alignment gaps (dashes): id=%s, occurences=%i', record.id, dash_count)
            if(is_dna):
                if(FASTA_DNA_SEQUENCE_PATTERN.fullmatch(raw_sequence) is None):
                    log.error('import: Fasta sequence contains invalid DNA characters! id=%s', record.id)
                    raise ValueError(f'Fasta sequence contains invalid DNA characters! id={record.id}')
                sequence['nt'] = raw_sequence
            else:
                if(raw_sequence[-1] == '*'):  # remove trailing stop asterik
                    raw_sequence = raw_sequence[:-1]
                    log.debug('import: Removed trailing asterik! id=%s, seq=%s', record.id, raw_sequence)
                if(FASTA_AA_SEQUENCE_PATTERN.fullmatch(raw_sequence) is None):
                    log.error('import: Fasta sequence contains invalid AA characters! id=%s, seq=%s', record.id, raw_sequence)
                    raise ValueError(f'Fasta sequence contains invalid AA characters! id={record.id}')
                sequence['aa'] = raw_sequence
            sequence['length'] = len(raw_sequence)
            if(is_genomic):
                sequence['complete'] = False
                sequence['type'] = bc.REPLICON_CONTIG
                sequence['topology'] = bc.TOPOLOGY_LINEAR
            log.info(
                'imported: id=%s, length=%i, description=%s, genomic=%s, dna=%s',
                sequence['id'], sequence['length'], sequence['description'], is_genomic, is_dna
            )
            sequences.append(sequence)
    return sequences


def export_sequences(sequences: Sequence[dict], fasta_path: Path, description: bool=False, wrap: bool=False):
    """Write sequences to Fasta file."""
    log.info('write genome sequences: path=%s, description=%s, wrap=%s', fasta_path, description, wrap)

    with fasta_path.open('wt') as fh:
        for seq in sequences:
            if(description):
                fh.write(f">{seq['id']} {seq['description']}\n")
            else:
                fh.write(f">{seq['id']}\n")
            if(wrap):
                fh.write(wrap_sequence(seq['nt'] if 'nt' in seq else seq['sequence']))  # <1.10.0 compatibility
            else:
                fh.write(seq['nt'])
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
    """Write annotated nucleotide sequences to Fasta file."""
    log.info('write feature nucleotide sequences: path=%s', ffn_path)
    with ffn_path.open('wt') as fh:
        for feat in features:
            if(feat['type'] in [bc.FEATURE_T_RNA, bc.FEATURE_TM_RNA, bc.FEATURE_R_RNA, bc.FEATURE_NC_RNA, bc.FEATURE_NC_RNA_REGION, bc.FEATURE_CRISPR, bc.FEATURE_CDS, bc.FEATURE_SORF, bc.FEATURE_ORIC, bc.FEATURE_ORIV, bc.FEATURE_ORIT]):
                identifier = feat['locus'] if 'locus' in feat else feat['id']
                if(feat.get('product', '') != ''):
                    fh.write(f">{identifier} {feat['product']}\n{feat['nt']}\n")
                else:
                    fh.write(f">{identifier}\n{feat['nt']}\n")

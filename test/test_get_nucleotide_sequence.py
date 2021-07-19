import pytest

import atexit
import shutil
import tempfile
import subprocess as sp
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

import bakta.config as cfg
import bakta.utils as bu
from bakta.features import cds as feat_cds
from bakta.io import fasta as fasta
import bakta.features.s_orf as s_orf


# Todo test for rRNA, ori, ncRNA, gaps, crispr
# Todo write cds and non cds sequence to fna


def prodigal_extract_nucletide_sequence(tmp, original_genome):
    """
    Generate test cases for all cds predicted by prodigal.
    """
    tmp_genome_path = tmp.joinpath('genome.fasta')
    fasta.export_contigs(original_genome, tmp_genome_path)
    nucleotide_path = tmp.joinpath('nucleotides.fna')

    cmd = [
        'prodigal',
        '-i', str(tmp_genome_path),
        '-d', str(nucleotide_path)
    ]

    proc = sp.run(
        cmd,
        cwd=str(tmp),
        env=cfg.env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        raise Exception(f'Prodigal error! error code: {proc.returncode}')

    sequences = []
    with open(nucleotide_path) as fh:
        for seqs in FastaIO.SimpleFastaParser(fh):
            header = [x.strip() for x in seqs[0].split('#')]
            sequences.append(
                ({'type': 'cds',
                  'contig': '_'.join(header[0].split('_')[:2]),
                  'start': int(header[1]),
                  'stop': int(header[2]),
                  'strand': '+' if header[3] == '1' else '-',
                  }, genome, ''.join([x.rstrip('\n') for x in seqs[1]]))
            )
    return sequences


def sorf_sequence_nt_aa(sorfs_entries):
    """
    Generate test cases for all sorfs.
    """
    ls = []
    for sorf in sorfs_entries:
        ls.append((sorf['nucleotide_sequence'], sorf['sequence']))
    return ls


tmp_path = Path(tempfile.mkdtemp())
atexit.register(shutil.rmtree, tmp_path)

cfg.min_contig_length = 1
cfg.translation_table = 11

genome, unused = bu.qc_contigs(fasta.import_contigs('test/data/GCF_000008865.2.fna.gz'), cfg.replicons)
cdss = prodigal_extract_nucletide_sequence(tmp_path, genome)  # (json entry, genome, nucleotide sequence)
sorfs = sorf_sequence_nt_aa(s_orf.extract({'contigs': genome}))


@pytest.mark.parametrize(
    'feature, contigs, expected',
    [
        # all cds nucleotide sequences output by prodigal -d
        *cdss,
        # cds on 5’ & 3’ contig edge; + strand
        ({'type': 'cds',
          'contig': 'contig_2',
          'start': 92563,
          'stop': 2502,
          'edge': True,
          'strand': '+'}, genome, 'ATGAAATTAAAGTATCTGTCATGTACGATCCTTGCCCCTCTGGCGATTGGGGTATTTTCTGCAACAGCTGCTGATAATAATTCA'
                                  'GCCATTTATTTCAATACCTCCCAGCCTATAAATGATCTGCAGGGTTCGTTGGCCGCAGAGGTGAAATTTGCACAAAGCCAGATT'
                                  'TTACCCGCCCATCCTAAAGAAGGGGATAGTCAACCACATCTGACCAGCCTGCGGAAAAGTCTGCTGCTTGTCCGTCCGGTGAAA'
                                  'GCTGATGATAAAACACCTGTTCAGGTGGAAGCCCGCGATGATAATAATAAAATTCTCGGTACGTTAACCCTTTATCCTCCTTCA'
                                  'TCACTACCGGATACAATCTACCATCTGGATGGTGTTCCGGAAGGTGGTATCGATTTCACACCTCATAATGGAACGAAAAAGATC'
                                  'ATTAATACGGTGGCTGAAGTAAACAAACTCAGTGATGCCAGCGGGAGTTCTATTCATAGCCATCTAACAAATAATGCACTGGTG'
                                  'GAGATCCATACTGCAAATGGTCGTTGGGTAAGAGACATTTATCTGCCGCAGGGACCCGACCTTGAAGGTAAGATGGTTCGCTTT'
                                  'GTTTCGTCTGCAGGCTATAGTTCAACGGTTTTTTATGGTGATCGAAAAGTCACACTCTCGGTGGGTAACACTCTTCTGTTCAAA'
                                  'TATGTAAATGGTCAGTGGTTCCGCTCCGGTGAACTGGAGAATAATCGAATCACTTATGCTCAGCATATTTGGAGTGCTGAACTG'
                                  'CCTGCGCACTGGATCGTGCCTGGTTTAAACTTGGTGATTAAACAGGGCAATCTGAGCGGTCGCCTAAATGATATCAAGATTGGA'
                                  'GCACCGGGTGAGCTGTTGTTGCATACAATTGATATCGGGATGTTGACCACTCCCCGGGATCGCTTTGATTTTGCCAAAGACAAA'
                                  'GAAGCACATAGGGAATATTTCCAGACCATTCCTGTAAGTCGTATGATTGTTAATAATTATGCGCCTCTACACCTAAAGGAAGTT'
                                  'ATGTTACCAACCGGAGAGTTATTGACAGATATGGATCCAGGAAATGGTGGGTGGCATAGTGGTACAATGCGTCAAAGAATAGGT'
                                  'AAAGAATTGGTTTCGCATGGCATTGATAATGCTAACTATGGTTTAAATAGTACCGCAGGCTTAGGGGAGAATAGTCATCCATAT'
                                  'GTAGTTGCGCAATTAGCGGCACATAATAGCCGCGGTAATTATGCTAATGGCATCCAGGTTCATGGTGGCTCCGGAGGTGGGGGA'
                                  'ATTGTTACTTTAGATTCCACATTGGGGAATGAGTTCAGTCATGAAGTTGGTCATAATTATGGTCTTGGTCATTATGTAGATGGT'
                                  'TTCAAGGGTTCTGTACATCGTAGTGCAGAAAATAACAACTCAACTTGGGGATGGGATGGTGATAAAAAACGGTTTATTCCTAAC'
                                  'TTTTATCCGTCTCAAACAAATGAAAAGAGTTGTCTGAATAATCAGTGTCAAGAACCGTTTGATGGACACAAATTTGGTTTTGAC'
                                  'GCCATGGCGGGAGGCAGCCCTTTCTCTGCTGCAAACCGTTTCACAATGTATACTCCGAATTCATCGGCTATCATCCAGCGTTTT'
                                  'TTTGAAAATAAAGCTGTGTTCGATAGCCGTTCCTCCACCGGCTTCAGCAAGTGGAATGCAGATACGCAGGAAATGGAACCGTAT'
                                  'GAACACACCATTGACCGTGCGGAGCAGATTACGGCTTCAGTCAATGAGCTAAGTGAAAGCAAAATGGCTGAGCTGATGGCAGAG'
                                  'TACGCTGTCGTCAAAGTGCATATGTGGAACGGTAACTGGACAAGAAACATCTATATCCCTACAGCCTCCGCAGATAATAGAGGC'
                                  'AGTATCCTGACCATCAACCATGAGGCCGGTTATAATAGTTATCTGTTTATAAATGGTGACGAAAAGGTCGTTTCCCAGGGGTAT'
                                  'AAAAAGAGCTTTGTTTCCGATGGTCAGTTCTGGAAAGAACGTGATGTGGTTGATACTCGTGAAGCGCGTAAGCCAGAGCAGTTT'
                                  'GGTGTTCCTGTGACGACTCTGGTGGGGTATTACGATCCGGAAGGCACGCTGTCAAGCTACATCTATCCTGCGATGTATGGTGCC'
                                  'TATGGCTTCACTTATTCCGATGATAGTCAGAATCTATCCGATAACGACTGCCAGCTGCAGGTGGATACGAAAGAAGGGCAGTTG'
                                  'CGATTCAGACTGGCTAATCACCGGGCTAACAACACTGTAATGAATAAGTTCCATATTAACGTGCCAACAGAAAGTCAGCCCACA'
                                  'CAGGCCACATTGGTTTGCAATAACAAGATACTGGATACCAAATCGCTCACACCTGCGCCAGAAGGACTTACCTATACTGTAAAT'
                                  'GGGCAGGCACTTCCAGCAAAAGAAAACGAGGGATGCATCGTGTCCGTGAATTCAGGTAAACGTTACTGTTTGCCGGTTGGTCAA'
                                  'CGGTCAGGATATAGCCTTCCTGACTGGATTGTTGGGCAGGAAGTCTATGTCGACAGCGGGGCTAAAGCGAAAGTGCTGCTTTCT'
                                  'GACTGGGATAACCTGTCCTATAACAGGATTGGTGAGTTTGTAGGTAATGTGAACCCAGCTGATATGAAAAAAGTTAAAGCCTGG'
                                  'AACGGACAGTATTTGGACTTCAGTAAACCTAGGTCAATGAGGGTTGTATATAAATAA')
    ]
)
def test_get_nucleotide_seqeunce(feature, contigs, expected):
    assert feat_cds.get_nucleotide_sequence(feature['start'], feature['stop'], feature['strand'], feature['contig'],
                                            contigs) == expected


@pytest.mark.parametrize(
    'sorf_dna, expected_sorf_aa',
    sorfs
)
def test_sorf_nucleotide_seqeunces(sorf_dna, expected_sorf_aa):
    assert Seq.translate(sorf_dna, table=cfg.translation_table,
                         stop_symbol='', to_stop=False, cds=False) == expected_sorf_aa

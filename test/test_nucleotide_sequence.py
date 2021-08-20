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
from bakta.io import fasta as fasta
import bakta.features.s_orf as s_orf
from bakta.features import cds as feat_cds


# Todo test for rRNA, ori, ncRNA, gaps, crispr

# Functions to create test cases

def sorf_sequence_nt_aa(sorfs_entries):
    """
    Generate test cases for all sorfs.
    """
    ls = []
    for sorf in sorfs_entries:
        ls.append((sorf['nt'], sorf['aa']))
    return ls


def prodigal_extract_nucletide_sequence(original_genome):
    """
    Generate test cases for all cds predicted by prodigal.
    """
    # create prodigal training file and get standard bakta cds predictions
    contigs_path = cfg.tmp_path.joinpath('contigs.fna')
    fasta.export_contigs(original_genome, contigs_path)
    all_cds = feat_cds.predict(genome, contigs_path)

    prodigal_ouputs = []
    cmds = []

    # execute prodigal for non-complete sequences (contigs)
    non_complete_contigs = [c for c in genome['contigs'] if not c['complete']]
    if len(non_complete_contigs) > 0:
        incomplete_contigs_path = cfg.tmp_path.joinpath('incomplete_contigs.fna')
        fasta.export_contigs(non_complete_contigs, incomplete_contigs_path)
        nucleotide_path = cfg.tmp_path.joinpath('incomplete_cds.fna')
        prodigal_ouputs.append(nucleotide_path)

        non_complete_cmd = [
            'prodigal',
            '-i', str(incomplete_contigs_path),
            '-g', str(cfg.translation_table),  # set translation table
            '-t', str(cfg.tmp_path.joinpath('prodigal.tf')),
            '-d', str(nucleotide_path),
            '-c'  # closed ends
        ]
        cmds.append(non_complete_cmd)

    # execute prodigal for complete replicons (chromosomes/plasmids)
    replicons = [c for c in genome['contigs'] if c['complete']]
    if len(replicons) > 0:
        replicons_path = cfg.tmp_path.joinpath('replicons.fasta')
        fasta.export_contigs(replicons, replicons_path)
        replicon_nucleotide_path = cfg.tmp_path.joinpath('replicon_cds.fna')
        prodigal_ouputs.append(replicon_nucleotide_path)

        cmd_replicon = [
            'prodigal',
            '-i', str(replicons_path),
            '-g', str(cfg.translation_table),  # set translation table
            '-t', str(cfg.tmp_path.joinpath('prodigal.tf')),
            '-d', str(replicon_nucleotide_path)
            ]
        cmds.append(cmd_replicon)

    # execute prodigal -d commands
    for c in cmds:
        proc = sp.run(
            c,
            cwd=str(cfg.tmp_path),
            env=cfg.env,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            raise Exception(f'Prodigal error! error code: {proc.returncode}')

    # extract nucleotide sequences from prodigal -d output
    nt = {c: {'+': {}, '-': {}} for c in clean_contigs.keys()}
    for p_out in prodigal_ouputs:
        with open(p_out) as fh:
            for seqs in FastaIO.SimpleFastaParser(fh):
                header = [x.strip() for x in seqs[0].split('#')]
                contig = '_'.join(header[0].split('_')[:2])
                strand = '+' if header[3] == '1' else '-'
                start = int(header[1])
                stop = int(header[2])
                nt[contig][strand][f'{start}-{stop}'] = ''.join([x.rstrip('\n') for x in seqs[1]])

    # create tests for all cds
    sequences = []
    for cds in all_cds:
        seq = ''
        if f"{cds['start']}-{cds['stop']}" in nt[cds['contig']][cds['strand']]:
            seq = nt[cds['contig']][cds['strand']][f"{cds['start']}-{cds['stop']}"]
        elif f"{cds['stop']}-{cds['start']}" in nt[cds['contig']][cds['strand']]:
            seq = nt[cds['contig']][cds['strand']][f"{cds['stop']}-{cds['start']}"]
        elif 'edge' in cds:
            print(cds['edge'])  # Todo test case is not predicted by prodigal -d; is added manually as test below
            print(f"{cds['start']}-{cds['stop']}" in nt[cds['contig']][cds['strand']])
            print(f"{cds['stop']}-{cds['start']}" in nt[cds['contig']][cds['strand']])
        else:
            print('error: cds was not predicted by prodigal -d')
            continue

        sequences.append(
            ({'type': 'cds',
              'contig': cds['contig'],
              'start': cds['start'],
              'stop': cds['stop'],
              'strand': cds['strand'],
              }, seq)
        )

    return sequences


# Create test cases.
# set necessary config options
cfg.tmp_path = Path(tempfile.mkdtemp())
atexit.register(shutil.rmtree, cfg.tmp_path)
cfg.min_contig_length = 1
cfg.translation_table = 11

# extract the test genome
contigs, complete_genome = bu.qc_contigs(fasta.import_contigs('test/data/GCF_000008865.2.fna.gz'), cfg.replicons)
clean_contigs = {c['id']: {'sequence': c['sequence']} for c in contigs}

# create dummmy genome
genome = {
    'genus': cfg.genus,
    'species': cfg.species,
    'strain': cfg.strain,
    'taxon': cfg.taxon,
    'gram': cfg.gram,
    'translation_table': cfg.translation_table,
    'size': sum([c['length'] for c in contigs]),
    'complete': cfg.complete or complete_genome,
    'features': {},
    'contigs': contigs
}

# get sORF nt and aa sequences (sorf['nt'], sorf['aa'])
sorfs = sorf_sequence_nt_aa(s_orf.extract(genome))

# generate cds tests
cdss = prodigal_extract_nucletide_sequence(contigs)  # (json entry, genome, nucleotide sequence)


@pytest.mark.parametrize(
    'sorf_dna, expected_sorf_aa',
    sorfs
)
def test_sorf_nucleotide_seqeunces(sorf_dna, expected_sorf_aa):
    assert Seq.translate(sorf_dna, table=cfg.translation_table,
                         stop_symbol='', to_stop=False, cds=False) == expected_sorf_aa


@pytest.mark.parametrize(
    'feature, expected',
    [
        # all cds nucleotide sequences output by prodigal -d
        *cdss,
        # CDS on 5’ & 3’ contig edge; forward strand
        ({'type': 'cds',
          'contig': 'contig_2',
          'start': 92563,
          'stop': 2502,
          'edge': True,
          'strand': '+'},
         'ATGAAATTAAAGTATCTGTCATGTACGATCCTTGCCCCTCTGGCGATTGGGGTATTTTCTGCAACAGCTGCTGATAATAATTCA'
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
def test_get_nucleotide_seqeunce(feature, expected):
    assert bu.extract_feature_sequence(feature, clean_contigs[feature['contig']]) == expected

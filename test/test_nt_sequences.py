import json

from pathlib import Path
from subprocess import run

CDS = 'TTGACTACGCCATTGAAAAAGATTGTGATTGTCGGCGGCGGTGCTGGTGGGCTGGAAATGGCAACACAGCTGGGGCATAAGCTGGGACGCAAGAAAAAAGCCAAAATTACGCTGGTCGATCGTAACCACAGCCATCTGTGGAAACCGCTGCTGCACGAAGTGGCGACTGGCTCGCTTGATGAAGGCGTCGATGCGTTGAGCTATCTGGCCCATGCGCGCAATCATGGTTTCCAGTTCCAGCTGGGTTCCGTCATTGATATTGATCGTGAAGCGAAAACAATCACTATTGCAGAACTGCGCGATGAGAAAGGTGAACTGCTGGTTCCGGAACGTAAAATCGCCTATGACACCCTGGTAATGGCGCTGGGTAGCACCTCTAACGATTTCAATACGCCAGGTGTCAAAGAGAACTGCATTTTCCTCGATAACCCGCACCAGGCGCGTCGCTTTCACCAGGAGATGCTGAATCTCTTCCTGAAATACTCCGCCAACCTGGGCGCAAATGGCAAAGTGAACATTGCGATTGTCGGCGGCGGCGCGACGGGTGTAGAACTCTCCGCTGAATTGCACAACGCGGTCAAGCAACTGCACAGCTACGGTTACAAAGGCCTGACCAACGAAGCCCTGAACGTAACGCTGGTAGAAGCGGGAGAACGTATTTTGCCTGCATTACCGCCACGTATCTCTGCTGCGGCCCACAACGAGCTAACGAAACTTGGCGTTCGCGTGCTGACGCAAACCATGGTCACCAGTGCTGATGAAGGCGGCCTGCACACTAAAGATGGCGAATATATTGAGGCTGATCTGATGGTGTGGGCAGCCGGGATCAAAGCGCCAGACTTCCTGAAAGATATCGGTGGTCTTGAAACTAACCGTATCAACCAGCTGGTGGTGGAACCGACGCTGCAAACCACCCGCGATCCAGACATTTACGCTATTGGCGACTGCGCGTCATGCCCGCGTCCGGAAGGGGGCTTTGTTCCGCCGCGTGCTCAGGCTGCACACCAGATGGCGACTTGCGCAATGAACAACATTCTGGCGCAGATGAATGGTAAGCCGCTGAAAAATTATCAGTATAAAGATCATGGTTCGCTGGTATCGCTGTCGAACTTCTCCACCGTTGGTAGCCTGATGGGTAACCTGACGCGCGGCTCAATGATGATTGAAGGACGAATTGCGCGCTTTGTATATATCTCGCTATACCGAATGCATCAGATTGCGCTGCATGGTTACTTTAAAACCGGATTAATGATGCTGGTGGGGAGTATTAACCGCGTTATCCGTCCGCGTTTGAAGTTGCATTAA'
SORF = 'ATGGTGAATACCGGCGGCAATAAACGTCAGGTGCCGGCGAAACGTCAGAATCGTGGCTCCCGTAATTCCAAAGATGATGGCGGCTAA'


def test_bakta_cds_nt_sequence(tmpdir):
    # test extracted nucleotide sequences of cds
    proc = run(
        [
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test',
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot',
            '--keep-contig-headers', '--complete', 'test/data/cds.fna'
        ]
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    output_path = tmpdir_path.joinpath('test.json')
    assert Path.exists(output_path)
    assert output_path.stat().st_size > 0

    with output_path.open() as fh:
        results = json.load(fh)

    for feat in results['features']:
        if(feat['sequence'] != 'dummy'):
            assert feat['nt'] == CDS


def test_bakta_sorf_nt_sequence(tmpdir):
    # test extracted nucleotide sequences of sorfs
    proc = run(
        [
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test',
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-ori', '--skip-gap', '--skip-plot',
            '--keep-contig-headers', '--complete', 'test/data/sorf.fna'
        ]
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    output_path = tmpdir_path.joinpath('test.json')
    assert Path.exists(output_path)
    assert output_path.stat().st_size > 0

    with output_path.open() as fh:
        results = json.load(fh)

    for feat in results['features']:
        assert feat['nt'] == SORF

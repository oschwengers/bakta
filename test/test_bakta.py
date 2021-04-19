import pytest

from pathlib import Path
from subprocess import run


FILES = [
    'test.log',
    'test.json',
    'test.tsv',
    'test.gff3',
    'test.gbff',
    'test.embl',
    'test.fna',
    'test.faa'
]


@pytest.mark.slow
def test_bakta_mock_skipped_features(tmpdir):
    # fast test skipping all feature detections
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.slow
def test_bakta_plasmid(tmpdir):
    # full test on plasmid
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.slow
def test_bakta_genome(tmpdir):
    # full test on plasmid
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/GCF_000008865.2.fna.gz'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))
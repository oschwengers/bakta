import pytest

from pathlib import Path
from subprocess import run

from .conftest import FILES, SKIP_PARAMETERS


@pytest.mark.slow
def test_bakta_mock_skipped_features(tmpdir):
    # fast test skipping all feature detections
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test'] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.slow
def test_bakta_plasmid(tmpdir):
    # full test on plasmid
    proc = run(['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--prefix', 'test', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0
        if 'tsv' in str(output_path):
            feat_count = count_features(output_path)
            assert feat_count == 3


@pytest.mark.slow
def test_bakta_genome(tmpdir):
    # full test on plasmid
    proc = run(['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--prefix', 'test', 'test/data/GCF_000008865.2.fna.gz'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0
        if 'tsv' in str(output_path):
            feat_count = count_features(output_path)
            assert feat_count == 5552

def count_features(file_path):
    with open (file_path, 'r') as fh:
        feat_count = 0
        for line in (line for line in fh if not line.startswith('#')):
            feat_count += 1
        return feat_count

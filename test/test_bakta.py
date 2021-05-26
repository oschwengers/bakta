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
            feat_test = {
                'tRNA': 0,
                'tmRNA': 0,
                'rRNA': 0,
                'ncRNA': 0,
                'ncRNA-region': 0,
                'crispr': 0,
                'sorf': 0,
                'oriV': 0,
                'oriC': 0,
                'oriT': 0,
                'cds': 3
            }
            feat_file = features_individual(output_path)
            for type in feat_file:
                assert feat_file[type] == feat_test[type]


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
            feat_test = {
                'tRNA': 107,
                'tmRNA': 1,
                'rRNA': 7,
                'ncRNA': 57,
                'ncRNA-region': 1,
                'crispr': 1,
                'sorf': 2,
                'oriV': 0,
                'oriC': 0,
                'oriT': 0,
                'cds': 5376
            }
            feat_file = features_individual(output_path)
            for type in feat_file:
                assert feat_file[type] == feat_test[type]

def count_features(file_path):
    with open (file_path, 'r') as fh:
        feat_count = 0
        for line in (line for line in fh if not line.startswith('#')):
            feat_count += 1
        return feat_count

def features_individual(file_path):
    with open (file_path, 'r') as fh:
        feat_file = {
            'tRNA': 0,
            'tmRNA': 0,
            'rRNA': 0,
            'ncRNA': 0,
            'ncRNA-region': 0,
            'crispr': 0,
            'sorf': 0,
            'oriV': 0,
            'oriC': 0,
            'oriT': 0,
            'cds': 0
        }
        for line in (line for line in fh if not line.startswith('#')):
            feature = line.split('\t')[1]
            if feature in feat_file:
                feat_file[feature] += 1
    return feat_file

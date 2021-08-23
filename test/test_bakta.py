from pathlib import Path
from subprocess import run

import pytest

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
    proc = run(['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--prefix', 'test', '--complete', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0

    output_path = tmpdir_path.joinpath('test.tsv')
    feature_count, feature_counts = count_features(output_path)
    assert feature_count == 3
    feature_counts_expected = {
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
    for type in feature_counts:
        assert feature_counts[type] == feature_counts_expected[type]


@pytest.mark.slow
def test_bakta_genome(tmpdir):
    # full test on complete genome in compliant mode
    proc = run(['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--prefix', 'test', '--complete', '--compliant', 'test/data/GCF_000008865.2.fna.gz'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0

    output_path = tmpdir_path.joinpath('test.tsv')
    feature_count, feature_counts = count_features(output_path)
    assert feature_count == 5551
    feature_counts_expected = {
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
        'cds': 5375
    }
    for type in feature_counts:
        assert feature_counts[type] == feature_counts_expected[type]


def count_features(file_path):
    with open(file_path, 'r') as fh:
        feature_count = 0
        feature_counts = {
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
        for line in fh:
            if not line.startswith('#'):
                feature_count += 1
                feature = line.split('\t')[1]
                if feature in feature_counts:
                    feature_counts[feature] += 1
    return feature_count, feature_counts

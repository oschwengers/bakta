import subprocess as sp
import sys

from pathlib import Path
from subprocess import run

import pytest

from .conftest import FILES, SKIP_PARAMETERS


def check_deepsig():
    command = ('deepsig', '--version')
    is_available = False
    try:
        # tool_output = str(sp.check_output(command, stderr=sp.STDOUT))  # stderr must be added in case the tool output is not piped into stdout
        sp.check_output(command, stderr=sp.STDOUT)  # stderr must be added in case the tool output is not piped into stdout
        is_available = True
    except:
        print(f"WARNING: {command[0]} not found or not executable! Skip dependen test.")
    return is_available


@pytest.mark.parametrize(
    'db',
    [
        ('db'),  # full DB
        ('db-light')  # light DB
    ]
)
def test_bakta_mock_skipped_features(db, tmpdir):
    # fast test full features
    proc = run(
        ['bin/bakta', '--db', f'test/{db}', '--output', tmpdir, '--force', '--prefix', 'test', '--min-contig-length', '200', '--proteins', 'test/data/user-proteins.faa'] +
        ['--genus', 'Foo gen. nov.', '--species', 'bar sp. nov.', '--strain', 'test 1'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.skipif(check_deepsig() == False, reason=f'Skip on unavailable DeepSig')
def test_bakta_plasmid(tmpdir):
    # full test on plasmid
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--force', '--prefix', 'test', '--min-contig-length', '200', '--complete', '--gram', '-', '--proteins', 'test/data/user-proteins.faa'] +
        ['--genus', 'Foo gen. nov.', '--species', 'bar sp. nov.', '--strain', 'test 1', 'test/data/NC_002127.1.fna.gz']
    )
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


@pytest.mark.parametrize(
    'db',
    [
        ('db'),  # full DB
        ('db-light')  # light DB
    ]
)
def test_bakta_genome(db, tmpdir):
    # full test on complete genome in compliant mode
    proc = run(['bin/bakta', '--db', f'test/{db}', '--verbose', '--output', tmpdir, '--force', '--force', '--prefix', 'test', '--min-contig-length', '200', '--complete', '--compliant', '--proteins', 'test/data/user-proteins.faa', '--genus', 'Foo gen. nov.', '--species', 'bar sp. nov.', '--strain', 'test 1', 'test/data/GCF_000008865.2.fna.gz'])
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

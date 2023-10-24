import json
import subprocess as sp

from pathlib import Path
from subprocess import run

import pytest

import bakta.constants as bc

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

    results_path = tmpdir_path.joinpath('test.json')
    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    features = results['features']
    assert len(features) == 3
    feature_counts_expected = {
        bc.FEATURE_T_RNA: 0,
        bc.FEATURE_TM_RNA: 0,
        bc.FEATURE_R_RNA: 0,
        bc.FEATURE_NC_RNA: 0,
        bc.FEATURE_NC_RNA_REGION: 0,
        bc.FEATURE_CRISPR: 0,
        bc.FEATURE_CDS: 3,
        bc.FEATURE_SORF: 0,
        bc.FEATURE_ORIC: 0,
        bc.FEATURE_ORIV: 0,
        bc.FEATURE_ORIT: 0
    }
    for type, count in feature_counts_expected.items():
        assert len([f for f in features if f['type'] == type]) == count


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

    results_path = tmpdir_path.joinpath('test.json')
    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    features = results['features']
    assert len(features) == 5551
    feature_counts_expected = {
        bc.FEATURE_T_RNA: 107,
        bc.FEATURE_TM_RNA: 1,
        bc.FEATURE_R_RNA: 7,
        bc.FEATURE_NC_RNA: 57,
        bc.FEATURE_NC_RNA_REGION: 1,
        bc.FEATURE_CRISPR: 1,
        bc.FEATURE_CDS: 5375,
        bc.FEATURE_SORF: 2,
        bc.FEATURE_ORIC: 0,
        bc.FEATURE_ORIV: 0,
        bc.FEATURE_ORIT: 0
    }
    for type, count in feature_counts_expected.items():
        assert len([f for f in features if f['type'] == type]) == count


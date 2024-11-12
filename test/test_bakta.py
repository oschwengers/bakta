import json
import subprocess as sp

from pathlib import Path
from subprocess import run

import pytest

import bakta.constants as bc

from .conftest import FILES, FILES_IO, SKIP_PARAMETERS


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
        assert len([feat for feat in features if feat['type'] == type]) == count


@pytest.mark.parametrize(
    'db',
    [
        pytest.param('db', marks=pytest.mark.dependency(name="test_bakta_genome")),  # full DB
        ('db-light')  # light DB
    ]
)
def test_bakta_genome(db, tmpdir, request):
    # full test on complete genome in compliant mode
    output_path = Path(tmpdir).joinpath(db)
    proc = run(
        [
            'bin/bakta',
            '--db', f'test/{db}',
            '--verbose',
            '--output', str(output_path),
            '--force',
            '--prefix', 'test',
            '--min-contig-length', '200',
            '--complete',
            '--genus', 'Foo gen. nov.',
            '--species', 'bar sp. nov.',
            '--strain', 'test 1',
            '--replicons', 'test/data/replicons.tsv',
            '--regions', 'test/data/NC_002127.1-region.gff3',
            '--proteins', 'test/data/user-proteins.faa',
            '--hmms', 'test/data/NF000185.2.HMM',
            '--compliant',
            'test/data/GCF_000008865.2.fna.gz']
    )
    assert proc.returncode == 0

    for file in FILES:
        file_path = output_path.joinpath(file)
        assert Path.exists(file_path)
        assert file_path.stat().st_size > 0

    results_path = output_path.joinpath('test.json')
    if(db == 'db'):
        request.config.cache.set('test_bakta_json', str(results_path))  # cache for later re-usage

    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    features = results['features']
    assert len(features) == 5550
    feature_counts_expected = {
        bc.FEATURE_T_RNA: 107,
        bc.FEATURE_TM_RNA: 1,
        bc.FEATURE_R_RNA: 7,
        bc.FEATURE_NC_RNA: 57,
        bc.FEATURE_NC_RNA_REGION: 1,
        bc.FEATURE_CRISPR: 1,
        bc.FEATURE_CDS: 5374,
        bc.FEATURE_SORF: 2,
        bc.FEATURE_ORIC: 0,
        bc.FEATURE_ORIV: 0,
        bc.FEATURE_ORIT: 0
    }
    for type, count in feature_counts_expected.items():
        assert len([feat for feat in features if feat['type'] == type]) == count


@pytest.mark.dependency(depends=['test_bakta_genome'])
def test_bakta_io(tmpdir, request):
    test_bakta_json = request.config.cache.get('test_bakta_json', None)
    assert test_bakta_json is not None
    result_path = Path(test_bakta_json)
    assert Path.exists(result_path)

    output_path = Path(tmpdir)
    proc = run(
        [
            'bin/bakta_io',
            '--verbose',
            '--output', str(output_path),
            '--force',
            '--prefix', 'test',
            test_bakta_json
        ]
    )
    assert proc.returncode == 0

    for file in FILES_IO:
        file_path = output_path.joinpath(file)
        assert Path.exists(file_path)
        assert file_path.stat().st_size > 0

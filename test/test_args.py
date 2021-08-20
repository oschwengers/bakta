import os
import pytest

from pathlib import Path
from subprocess import run

from .conftest import FILES, SKIP_PARAMETERS


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo.fasta']),  # not existing
        (['test/data/invalid.fasta'])  # invalid fasta DNA alphabet
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/bakta', '--db', 'test/db', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['--db']),  # missing path
        (['--db', '', ]),  # empty
        (['--db', 'test/foo']),  # not existing
    ]
)
def test_database_failing_parameter(parameters, tmpdir):
    # test database arguments

    cmd_line = ['bin/bakta', '--output', tmpdir] + parameters + ['test/data/NC_002127.1.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'env_key,env_value',
    [
        ('foo', ''),  # not provided
        ('BAKTA_DB', ''),  # missing path
        ('BAKTA_DB', 'test/foo')  # not existing path
    ]
)
def test_database_failing_environment(env_key, env_value, tmpdir):
    # test database arguments

    env = os.environ
    env[env_key] = env_value
    cmd_line = ['bin/bakta', '--output', tmpdir, 'test/data/NC_002127.1.fna']
    proc = run(cmd_line, env=env)
    assert proc.returncode != 0


def test_database_ok(tmpdir):
    # test database arguments

    # parameter OK
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    # environment OK
    env = os.environ
    env['BAKTA_DB'] = 'test/db'
    proc = run(['bin/bakta', '--output', tmpdir] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'], env=env)
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--prodigal-tf']),  # not provided
        (['--prodigal-tf', '']),  # empty
        (['--prodigal-tf', 'foo'])  # not existing
    ]
)
def test_prodigal_tf_failiing(parameters, tmpdir):
    # test prodigal training file arguments

    # missing path
    cmd_line = ['bin/bakta', '--db', 'test/db', '--output', tmpdir] + parameters + ['test/data/NC_002127.1.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.slow
def test_prodigal_tf_ok(tmpdir):
    # test prodigal training file arguments

    # OK
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/prodigal.tf'] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.parametrize(
    'parameters',
    [
        (['--replicons']),  # not provided
        (['--replicons', '']),  # empty
        (['--replicons', 'foo'])  # not existing
    ]
)
def test_replicons_failiing(parameters, tmpdir):
    # test replicons file arguments

    # missing path
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir] + parameters + ['test/data/NC_002127.1.fna'])
    assert proc.returncode != 0


@pytest.mark.slow
def test_replicons_ok(tmpdir):
    # test replicons file arguments

    # OK: replicons as CSV
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'test/data/replicons.csv'] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))

    # OK: replicons as TSV
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'test/data/replicons.tsv'] + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


def test_output_failing():
    # test database arguments
    cmd_line = ['bin/bakta', '--output', '/', 'test/data/draft-w-plasmids.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--threads']),  # not provided
        (['--threads', '']),  # empty
        (['--threads', 'foo']),  # string
        (['--threads', '-1']),  # smaller than zero
        (['--threads', '0']),  # zero
        (['--threads', '1.1']),  # float
        (['--threads', '1000'])  # larger than available threads
    ]
)
def test_threads_failing(parameters, tmpdir):
    # test threads arguments
    cmd_line = ['bin/bakta', '--db', 'test/db', '--output', tmpdir] + parameters + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--min-contig-length']),  # not provided
        (['--min-contig-length', '']),  # empty
        (['--min-contig-length', 'foo']),  # string
        (['--min-contig-length', '-1']),  # smaller than zero
        (['--min-contig-length', '0']),  # zero
        (['--min-contig-length', '1.1']),  # float
    ]
)
def test_min_contig_length_failing(parameters, tmpdir):
    # test min-contig-length arguments
    cmd_line = ['bin/bakta', '--db', 'test/db', '--output', tmpdir] + parameters + SKIP_PARAMETERS + ['test/data/NC_002127.1.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0

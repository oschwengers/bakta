import os
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


@pytest.mark.parametrize(
    'parameters',
    [
        (['--db', 'test/db']),  # not provided
        (['--db', 'test/db', '']),  # empty
        (['--db', 'test/db', 'foo.fasta'])  # not existing
    ]
)
def test_genome_failing(parameters):
    # test genome arguments
    cmd_line = ['bin/bakta'] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters,env_key,env_value',
    [
        (['test/data/NC_002127.1.fna'], 'foo', ''),  # not provided
        (['--db', 'test/data/NC_002127.1.fna'], 'foo', ''),  # parameter missing path
        (['--db', '', 'test/data/NC_002127.1.fna'], 'foo', ''),  # parameter empty
        (['--db', 'test/foo', 'test/data/NC_002127.1.fna'], 'foo', ''),  # parameter not existing
        (['test/data/NC_002127.1.fna'], 'BAKTA_DB', ''),  # environment missing path
        (['test/data/NC_002127.1.fna'], 'BAKTA_DB', 'test/foo')  # environment not existing
    ]
)
def test_database_failing(parameters, env_key, env_value):
    # test database arguments

    env = os.environ
    env[env_key] = env_value
    cmd_line = ['bin/bakta'] + parameters
    proc = run(cmd_line, env=env)
    assert proc.returncode != 0


def test_database_ok(tmpdir):
    # test database arguments

    # parameter OK
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    # environment OK
    env = os.environ
    env['BAKTA_DB'] = 'test/db'
    proc = run(['bin/bakta', '--output', tmpdir, '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'], env=env)
    assert proc.returncode == 0


def test_prodigal_tf_failiing(tmpdir):
    # test prodigal training file arguments

    # missing path
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # empty
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', '', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # not existing
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'foo', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0


@pytest.mark.slow
def test_prodigal_tf_ok(tmpdir):
    # test prodigal training file arguments

    # OK
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/prodigal.tf', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


def test_replicons_failiing(tmpdir):
    # test replicons file arguments

    # missing path
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # empty
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', '', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # not existing
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'foo', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0


@pytest.mark.slow
def test_replicons_ok(tmpdir):
    # test replicons file arguments

    # OK: replicons as CSV
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'test/data/replicons.csv', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))

    # OK: replicons as TSV
    proc = run(['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--replicons', 'test/data/replicons.tsv', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


def test_output_failing():
    # test database arguments
    cmd_line = ['bin/bakta', '--output', '/', 'test/data/draft-w-plasmids.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0
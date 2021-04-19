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
    "cmd_line",
    [
        (["bin/bakta", '--db', 'test/db']),  # no provided
        (["bin/bakta", '--db', 'test/db', 'foo.fasta'])  # non-existing
    ]
)
def test_genome_failing(cmd_line):
    # test genome arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "cmd_line,env_key,env_value",
    [
        (["bin/bakta", 'test/data/NC_002127.1.fna'], '', ''),  # no provided
        (["bin/bakta", '--db', 'test/db', 'foo.fasta'], '', ''),  # non-existing
        (["bin/bakta", '--db', 'test/data/NC_002127.1.fna'], '', ''),  # parameter missing path
        (["bin/bakta", '--db', 'test/', 'test/data/NC_002127.1.fna'], '', ''),  # parameter wrong path
        (["bin/bakta", '--db', 'test/', 'test/data/NC_002127.1.fna'], 'BAKTA_DB', ''),  # environment missing path
        (["bin/bakta", '--db', 'test/', 'test/data/NC_002127.1.fna'], 'BAKTA_DB', 'test/'),  # parameter wrong path
    ]
)
def test_database_failing(cmd_line, env_key, env_value):
    # test database arguments

    env = os.environ
    env[env_key] = env_value
    proc = run(cmd_line, env=env)
    assert proc.returncode != 0


def test_database_ok(tmpdir):
    # test database arguments

    # parameter OK
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    # parameter OK
    env = os.environ
    env['BAKTA_DB'] = 'test/db'
    proc = run(["bin/bakta", '--output', tmpdir, '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'], env=env)
    assert proc.returncode == 0


def test_prodigal_tf(tmpdir):
    # test prodigal training file arguments

    # missing
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # non-existing
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/foo', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # OK
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/prodigal.tf', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))
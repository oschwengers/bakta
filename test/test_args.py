import os
from pathlib import Path
from subprocess import run


def test_genome():
    # test genome arguments

    # no provided
    proc = run(["bin/bakta", '--db', 'test/db'])
    assert proc.returncode != 0

    # non-existing
    proc = run(["bin/bakta", '--db', 'test/db', 'foo.fasta'])
    assert proc.returncode != 0


def test_database(tmpdir):
    # test database arguments

    # no provided
    proc = run(["bin/bakta", 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # parameter missing path
    proc = run(["bin/bakta", '--db', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # parameter wrong path
    proc = run(["bin/bakta", '--db', 'test/', 'test/data/NC_002127.1.fna'])
    assert proc.returncode != 0

    # parameter OK
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0

    # environment missing path
    env = os.environ
    env['BAKTA_DB'] = ''
    proc = run(["bin/bakta", 'test/data/NC_002127.1.fna'], env=env)
    assert proc.returncode != 0

    # parameter wrong path
    env['BAKTA_DB'] = 'test/'
    proc = run(["bin/bakta", 'test/data/NC_002127.1.fna'], env=env)
    assert proc.returncode != 0

    # parameter OK
    env['BAKTA_DB'] = 'test/db'
    proc = run(["bin/bakta", '--output', tmpdir, 'test/data/NC_002127.1.fna'], env=env)
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
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--prodigal-tf', 'test/data/prodigal.tf', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0
    assert Path.exists(Path(tmpdir).joinpath('test.log'))
    assert Path.exists(Path(tmpdir).joinpath('test.json'))
    assert Path.exists(Path(tmpdir).joinpath('test.tsv'))
    assert Path.exists(Path(tmpdir).joinpath('test.gff3'))
    assert Path.exists(Path(tmpdir).joinpath('test.gbff'))
    assert Path.exists(Path(tmpdir).joinpath('test.fna'))
    assert Path.exists(Path(tmpdir).joinpath('test.faa'))
from pathlib import Path
from subprocess import run


def test_bakta_mock_skipped_features(tmpdir):
    # fast test skipping all feature detections
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0
    assert Path.exists(Path(tmpdir).joinpath('test.log'))
    assert Path.exists(Path(tmpdir).joinpath('test.json'))
    assert Path.exists(Path(tmpdir).joinpath('test.tsv'))
    assert Path.exists(Path(tmpdir).joinpath('test.gff3'))
    assert Path.exists(Path(tmpdir).joinpath('test.gbff'))
    assert Path.exists(Path(tmpdir).joinpath('test.fna'))
    assert Path.exists(Path(tmpdir).joinpath('test.faa'))


def test_bakta_plasmid(tmpdir):
    # full test on plasmid
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/NC_002127.1.fna'])
    assert proc.returncode == 0
    assert Path.exists(Path(tmpdir).joinpath('test.log'))
    assert Path.exists(Path(tmpdir).joinpath('test.json'))
    assert Path.exists(Path(tmpdir).joinpath('test.tsv'))
    assert Path.exists(Path(tmpdir).joinpath('test.gff3'))
    assert Path.exists(Path(tmpdir).joinpath('test.gbff'))
    assert Path.exists(Path(tmpdir).joinpath('test.fna'))
    assert Path.exists(Path(tmpdir).joinpath('test.faa'))


def test_bakta_genome(tmpdir):
    # full test on plasmid
    proc = run(["bin/bakta", '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/GCF_000008865.2.fna.gz'])
    assert proc.returncode == 0
    assert Path.exists(Path(tmpdir).joinpath('test.log'))
    assert Path.exists(Path(tmpdir).joinpath('test.json'))
    assert Path.exists(Path(tmpdir).joinpath('test.tsv'))
    assert Path.exists(Path(tmpdir).joinpath('test.gff3'))
    assert Path.exists(Path(tmpdir).joinpath('test.gbff'))
    assert Path.exists(Path(tmpdir).joinpath('test.fna'))
    assert Path.exists(Path(tmpdir).joinpath('test.faa'))
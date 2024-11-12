import os

from pathlib import Path
from subprocess import run

import pytest

from .conftest import FILES, SKIP_PARAMETERS


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo.fasta']),  # not existing
        (['fo o.fasta']),  # not existing (whitespace)
        (['test/data/empty']),  # empty file
        (['test/data/nonsense.txt']),  # nonsense file
        (['test/data/invalid.fasta']),  # invalid fasta DNA alphabet
        (['test/data/NC_002127.1.fna', 'foo']),  # additional argument
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        parameters
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'genome',
    [
        (['test/data/NC_002127.1.fna']),  # Linux/Unix format
        (['test/data/NC_002127.1.fna.gz']),  # Linux/Unix format
        (['test/data/NC_002127.1-win.fna']),  # Windows format (\r\n)
        (['test/data/NC_002127.1-win.fna.gz']),  # Windows format (\r\n)
        (['test/data/valid.fasta'])  # valid minimal DNA sequences
    ]
)
def test_genome_ok(genome, tmpdir):
    # test genome arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        genome
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['--db']),  # missing path
        (['--db', '', ]),  # empty
        (['--db', 'test/foo']),  # not existing
        (['--db', 'test/data/empty']),  # empty file
        (['--db', 'test/data/nonsense.txt'])  # nonsense file
    ]
)
def test_database_failing_parameter(parameters, tmpdir):
    # test database arguments
    proc = run(
        ['bin/bakta', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'env_key,env_value',
    [
        ('foo', ''),  # not provided
        ('BAKTA_DB', ''),  # missing path
        ('BAKTA_DB', 'test/foo'),  # not existing path
        ('BAKTA_DB', 'test/data/empty'),  # empty file
        ('BAKTA_DB', 'test/data/nonsense.txt')  # nonsense file
    ]
)
def test_database_failing_environment(env_key, env_value, tmpdir):
    # test database arguments

    env = os.environ
    env[env_key] = env_value
    proc = run(
        ['bin/bakta', '--output', tmpdir, '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna'],
        env=env
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'db',
    [
        ('db'),  # full DB
        ('db-light')  # light DB
    ]
)
def test_database_ok(db, tmpdir):
    # test database arguments

    # parameter OK
    proc = run(
        ['bin/bakta', '--db', f'test/{db}', '--output', tmpdir, '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0

    # environment OK
    env = os.environ
    env['BAKTA_DB'] = f'test/{db}'
    proc = run(
        ['bin/bakta', '--output', tmpdir, '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
        , env=env
    )
    assert proc.returncode == 0


def test_output_failing():
    # test database arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', '/', '--force', '--skip-plot'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--tmp-dir'])  # not provided
    ]
)
def test_tmp_dir_failiing(parameters, tmpdir):
    # test tmp dir arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


def test_tmp_dir_ok(tmpdir):
    # test tmp dir arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--tmp-dir', '', '--skip-plot'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--prodigal-tf']),  # not provided
        (['--prodigal-tf', '']),  # empty
        (['--prodigal-tf', 'foo']),  # not existing
        (['--prodigal-tf', 'test/data/empty'])  # empty file
    ]
)
def test_prodigal_tf_failiing(parameters, tmpdir):
    # test prodigal training file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force'] +
        parameters +
        ['--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot'] +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


def test_prodigal_tf_ok(tmpdir):
    # test prodigal training file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test', '--prodigal-tf', 'test/data/prodigal.tf'] +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.parametrize(
    'parameters',
    [
        (['--replicons']),  # not provided
        (['--replicons', '']),  # empty
        (['--replicons', 'foo']),  # not existing
        (['--replicons', 'test/data/empty']),  # empty file
        (['--replicons', 'test/data/nonsense.txt'])  # nonsense file
    ]
)
def test_replicons_failiing(parameters, tmpdir):
    # test replicons file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--replicons', 'test/data/replicons.csv']),  # CSV
        (['--replicons', 'test/data/replicons.tsv'])  # TSV
    ]
)
def test_replicons_ok(parameters, tmpdir):
    # test replicons file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        assert Path.exists(tmpdir_path.joinpath(file))


@pytest.mark.parametrize(
    'parameters',
    [
        (['--regions']),  # not provided
        (['--regions', '']),  # empty
        (['--regions', 'foo']),  # not existing
        (['--regions', 'test/data/empty'])  # empty file
    ]
)
def test_regions_failiing(parameters, tmpdir):
    # test regions file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--proteins']),  # not provided
        (['--proteins', '']),  # empty
        (['--proteins', 'foo']),  # not existing
        (['--proteins', 'test/data/empty']),  # empty file
        (['--proteins', 'test/data/nonsense.txt'])  # nonsense file
    ]
)
def test_proteins_failiing(parameters, tmpdir):
    # test proteins file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force'] + 
        parameters + 
        ['--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot'] +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--hmms']),  # not provided
        (['--hmms', '']),  # empty
        (['--hmms', 'foo']),  # not existing
        (['--hmms', 'test/data/empty']),  # empty file
        (['--hmms', 'test/data/nonsense.txt'])  # nonsense file
    ]
)
def test_hmms_failiing(parameters, tmpdir):
    # test HMM file arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force'] + 
        parameters + 
        ['--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot'] +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus']),  # not provided
        (['--locus', '']),  # empty
        (['--locus', ' ']),  # whitespace only
        (['--locus', '  ']),  # whitespaces only
        (['--locus', 'fo o']),  # containing whitespace
        (['--locus', 'ABCDEFGHIJKLMNOPQRSTU'])  # more than 20 characters
    ]
)
def test_locus_failiing(parameters, tmpdir):
    # test locus prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus', 'ABC']),
        (['--locus', 'ABCDEFGHIJKLMNOPQRST']),
        (['--locus', 'A123_.*-:#'])
    ]
)
def test_locus_ok(parameters, tmpdir):
    # test locus prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag']),  # not provided
        (['--locus-tag', '']),  # empty
        (['--locus-tag', ' ']),  # whitespace only
        (['--locus-tag', '  ']),  # whitespaces only
        (['--locus-tag', 'ABCDEFGHIJKLMNOPQRSTUVWXZ']),  # more than 24 characters
        (['--locus-tag', 'ABC!']),  # wrong characters
        (['--locus-tag', 'ABC?']),  # wrong characters
        (['--locus-tag', 'ABC*']),  # wrong characters
        (['--locus-tag', 'ABC,']),  # wrong characters
        (['--locus-tag', 'ABC;']),  # wrong characters
        (['--locus-tag', 'ABC:']),  # wrong characters
        (['--locus-tag', 'ABC§']),  # wrong characters
        (['--locus-tag', 'ABC$']),  # wrong characters
        (['--locus-tag', 'ABC%']),  # wrong characters
        (['--locus-tag', 'ABC&']),  # wrong characters
        (['--locus-tag', 'ABC/']),  # wrong characters
        (['--locus-tag', 'ABC=']),  # wrong characters
        (['--locus-tag', 'ABC#'])  # wrong characters
    ]
)
def test_locustag_failiing(parameters, tmpdir):
    # test locus-tag prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag', 'A']),
        (['--locus-tag', '1']),
        (['--locus-tag', 'ABCDEFGHIJKLMNOPQRSTUVWX']),
        (['--locus-tag', 'A12']),
        (['--locus-tag', 'ABC.']),
        (['--locus-tag', 'ABC-']),
        (['--locus-tag', 'ABC_']),
        (['--locus-tag', 'GCF_014267685.1']),
        (['--locus-tag', 'ASM25969v1']),
        (['--locus-tag', 'DAESDI010000001.1'])
    ]
)
def test_locustag_ok(parameters, tmpdir):
    # test locus-tag prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag']),  # not provided
        (['--locus-tag', '']),  # empty
        (['--locus-tag', ' ']),  # whitespace only
        (['--locus-tag', '  ']),  # whitespaces only
        (['--locus-tag', 'fo o']),  # containing whitespace
        (['--locus-tag', '123ABC']),  # first character is not a letter
        (['--locus-tag', 'abc']),  # lower case letters
        (['--locus-tag', 'AB']),  # less than 3 characters
        (['--locus-tag', 'ABCDEFGHIJKLM']),  # more than 12 characters
        (['--locus-tag', 'ABC_']),  # wrong characters
        (['--locus-tag', 'ABC-']),  # wrong characters
        (['--locus-tag', 'ABC.']),  # wrong characters
        (['--locus-tag', 'ABC!']),  # wrong characters
        (['--locus-tag', 'ABC?']),  # wrong characters
        (['--locus-tag', 'ABC*']),  # wrong characters
        (['--locus-tag', 'ABC,']),  # wrong characters
        (['--locus-tag', 'ABC;']),  # wrong characters
        (['--locus-tag', 'ABC:']),  # wrong characters
        (['--locus-tag', 'ABC§']),  # wrong characters
        (['--locus-tag', 'ABC$']),  # wrong characters
        (['--locus-tag', 'ABC%']),  # wrong characters
        (['--locus-tag', 'ABC&']),  # wrong characters
        (['--locus-tag', 'ABC/']),  # wrong characters
        (['--locus-tag', 'ABC=']),  # wrong characters
        (['--locus-tag', 'ABC#'])  # wrong characters
    ]
)
def test_locustag_compliant_failiing(parameters, tmpdir):
    # test locus-tag prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--compliant', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag', 'ABC']),
        (['--locus-tag', 'ABCDEFGHIJKL']),
        (['--locus-tag', 'A12']),
        (['--locus-tag', 'A23456789012'])
    ]
)
def test_locustag_compliant_ok(parameters, tmpdir):
    # test locus-tag prefix arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--compliant', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag-increment']),  # not provided
        (['--locus-tag-increment', '']),  # empty
        (['--locus-tag-increment', ' ']),  # whitespace only
        (['--locus-tag-increment', 'A']),  # wrong characters
        (['--locus-tag-increment', 'a']),  # wrong characters
        (['--locus-tag-increment', '0']),  # wrong number
        (['--locus-tag-increment', '11']),  # wrong number
    ]
)
def test_locustag_increment_failiing(parameters, tmpdir):
    # test locus-tag increment arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--locus-tag-increment', '1']),
        (['--locus-tag-increment', '5']),
        (['--locus-tag-increment', '10'])
    ]
)
def test_locustag_increment_ok(parameters, tmpdir):
    # test locus-tag increment arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--genus']),  # not provided
        (['--genus', '']),  # empty
        (['--genus', ' ']),  # whitespace only
        (['--genus', '  '])  # whitespaces only
    ]
)
def test_genus_failiing(parameters, tmpdir):
    # test genus arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--species', 'sp.']),
        (['--species', 'sp. nov.']),
        (['--species', 'a']),
        (['--species', '1']),
        (['--species', 'az_123'])
    ]
)
def test_species_ok(parameters, tmpdir):
    # test species arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--species']),  # not provided
        (['--species', '']),  # empty
        (['--species', ' ']),  # whitespace only
        (['--species', '  '])  # whitespaces only
    ]
)
def test_species_failiing(parameters, tmpdir):
    # test species arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--strain', 'a']),
        (['--strain', '1']),
        (['--strain', 'az. 123-123_234/123 = 123']),  # full blown crazy strain designation
        (['--strain', "0123456789'azAZ/.,#+* -_:(1)[]"])  # all RefSeq bacterial strain characters
    ]
)
def test_strain_ok(parameters, tmpdir):
    # test strain arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--strain']),  # not provided
        (['--strain', '']),  # empty
        (['--strain', ' ']),  # whitespace only
        (['--strain', '  '])  # whitespaces only
    ]
)
def test_strain_failiing(parameters, tmpdir):
    # test strain arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--plasmid', 'unnamed']),  # unnamed example
        (['--plasmid', 'unnamed1']),  # unnamed example
        (['--plasmid', 'unnamed12']),  # unnamed example
        (['--plasmid', 'unnamed123']),  # unnamed example
        (['--plasmid', 'pA']),  # minimal set
        (['--plasmid', 'pZ']),  # minimal set
        (['--plasmid', 'p1']),  # minimal set
        (['--plasmid', 'p2']),  # minimal set
        (['--plasmid', 'p.']),  # minimal set
        (['--plasmid', 'p_']),  # minimal set
        (['--plasmid', 'pAZ.az_09'])  # full blown plasmid designation
    ]
)
def test_plasmid_ok(parameters, tmpdir):
    # test plasmid arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode == 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--plasmid', '']),  # empty
        (['--plasmid', 'unnamed1234']),  # wrong unnamed example
        (['--plasmid', 'unname']),  # wrong unnamed example
        (['--plasmid', 'p']),  # p only
        (['--plasmid', 'pABCDEFGHIJKLMNOPQRST']),  # longer than 20 chars
        (['--plasmid', 'pAZ 123']),  # contains whitespace
        (['--plasmid', 'pAZ-123']),  # contains dash
        (['--plasmid', 'pAZ/123']),  # contains slash
        (['--plasmid', 'pAZ\\123']),  # contains backslash
        (['--plasmid', 'pAZ=123']),  # contains equals
        (['--plasmid', 'plasmid123']),  # contains word plasmid
        (['--plasmid', 'pPlasmidA1'])  # contains word plasmid
    ]
)
def test_plasmid_failiing(parameters, tmpdir):
    # test plasmid arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--gram']),  # not provided
        (['--gram', '']),  # empty
        (['--gram', 'foo']),  # wrong string
        (['--gram', '-1']),  # number
        (['--gram', '0']),  # number
        (['--gram', '1.1']),  # number
    ]
)
def test_threads_failing(parameters, tmpdir):
    # test gram arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--threads']),  # not provided
        (['--threads', '']),  # empty
        (['--threads', 'foo']),  # string
        (['--threads', '-1']),  # smaller than zero
        (['--threads', '1.1']),  # float
        (['--threads', '1000'])  # larger than available threads
    ]
)
def test_threads_failing(parameters, tmpdir):
    # test threads arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
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
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['--translation-table']),  # not provided
        (['--translation-table', '']),  # empty
        (['--translation-table', 'foo']),  # string
        (['--translation-table', '-1']),  # smaller than zero
        (['--translation-table', '0']),  # zero
        (['--translation-table', '1.1']),  # float
    ]
)
def test_min_contig_length_failing(parameters, tmpdir):
    # test min-contig-length arguments
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--skip-plot'] +
        parameters +
        SKIP_PARAMETERS +
        ['test/data/NC_002127.1.fna']
    )
    assert proc.returncode != 0

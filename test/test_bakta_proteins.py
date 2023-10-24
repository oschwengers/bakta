from pathlib import Path
from subprocess import run

import pytest


FILES = [
    'test.tsv',
    'test.hypotheticals.tsv',
    'test.faa'
]


@pytest.mark.parametrize(
    "user_proteins",
    [
        'test/data/user-proteins.faa',
        'test/data/user-proteins.gbff'
    ]
)
@pytest.mark.parametrize(
    "input",
    [
        'test/data/GCF_000008865.2.faa',
        'test/data/GCF_000008865.2.faa.gz'
    ]
)
def test_bakta_genome(user_proteins, input, tmpdir):
    # full test on complete genome in compliant mode
    proc = run(['bin/bakta_proteins', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test', '--proteins', user_proteins, input])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0

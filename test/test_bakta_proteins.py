from pathlib import Path
from subprocess import run

import pytest

FILES = [
    'test.tsv',
    'test.hypotheticals.tsv',
    'test.faa'
]


@pytest.mark.slow
@pytest.mark.parametrize(
    "parameters",
    [
        'test/data/user-proteins.faa',
        'test/data/user-proteins.gbff'
    ]
)
def test_bakta_genome(parameters, tmpdir):
    # full test on complete genome in compliant mode
    proc = run(['bin/bakta_proteins', '--db', 'test/db', '--output', tmpdir, '--tmp-dir', tmpdir, '--prefix', 'test', '--proteins', parameters, 'test/data/GCF_000008865.2.faa'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0

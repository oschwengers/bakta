import json
import subprocess as sp

from pathlib import Path
from subprocess import run

import pytest

from .conftest import FILES, SKIP_PARAMETERS


@pytest.mark.parametrize(
    'regions',
    [
        ('NC_002127.1-region.gff3'),  # GFF3
        ('NC_002127.1-region.gbff')  # Genbank
    ]
)
def test_bakta_plasmid(regions, tmpdir):
    # full test on plasmid
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--force', '--prefix', 'test', '--regions', f'test/data/{regions}'] +
        ['test/data/NC_002127.1.fna.gz']
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
    cdss = [feat for feat in results['features'] if feat['type'] == 'cds']
    assert len(cdss) == 3
    for cds in cdss:
        if(cds['strand'] == '+'  and  cds['stop'] == 736):  # test case 1: alternative downstream start codon +[32,736] superseding de novo-predicted CDS +[2,736]
            assert cds['start'] != 2  # de novo-predicted CDS start
            assert cds['start'] == 32  # user-provided CDS start
        elif(cds['strand'] == '-'  and  cds['start'] == 1348):  # test case 2: alternative downstream start codon -[1348,2229] superseding de novo-predicted CDS -[1348,2388]
            assert cds['stop'] != 2388  # de novo-predicted CDS start
            assert cds['stop'] == 2229  # user-provided CDS start
        elif(cds['strand'] == '+'  and  cds['start'] == 981):  # test case 3: user-provided CDS +[981,1121] overlapping de novo-predicted CDS -[971,1351]
            assert cds['stop'] == 1121  # user-provided overlapping CDS -> de novo-predicted was filtered-out

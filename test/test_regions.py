import json
import subprocess as sp

from pathlib import Path
from subprocess import run

import pytest

from .conftest import FILES


def test_wrong_seq_id_failiing(tmpdir):
    # CDS test on plasmid
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--force', '--prefix', 'test', '--regions', f'test/data/NC_002127.1-region.wrong-ids.gff3'] +
        ['--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot'] +
        ['test/data/NC_002127.1.fna.gz']
    )
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'regions',
    [
        ('NC_002127.1-region.gff3'),  # GFF3
        ('NC_002127.1-region.gbff')  # Genbank
    ]
)
@pytest.mark.parametrize(
    'keep_contig_headers',
    [
        ([]),  # autogenerate sequence ids
        (['--keep-contig-headers'])  # keep sequence headers
    ]
)
def test_regions_plasmid(regions, keep_contig_headers, tmpdir):
    # CDS test on plasmid
    proc = run(
        ['bin/bakta', '--db', 'test/db', '--verbose', '--output', tmpdir, '--force', '--prefix', 'test', '--regions', f'test/data/{regions}'] + keep_contig_headers +
        ['--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot'] +
        ['test/data/NC_002127.1.fna.gz']
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    results_path = tmpdir_path.joinpath('test.json')
    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    cdss = [feat for feat in results['features'] if feat['type'] == 'cds']
    assert len(cdss) == 3

    # test case 1: alternative downstream start codon +[32,736] superseding de novo-predicted CDS +[2,736]
    provided_cds_1 = [cds for cds in cdss if cds['start'] == 32 and cds['stop'] == 736 and cds['strand'] == '+']
    assert len(provided_cds_1) == 1
    denovo_cds_1 = [cds for cds in cdss if cds['start'] == 2 and cds['stop'] == 736 and cds['strand'] == '+']
    assert len(denovo_cds_1) == 0  # de novo-predicted was filtered-out

    # test case 2: alternative downstream start codon -[1348,2229] superseding de novo-predicted CDS -[1348,2388]
    provided_cds_2 = [cds for cds in cdss if cds['start'] == 1348 and cds['stop'] == 2229 and cds['strand'] == '-']
    assert len(provided_cds_2) == 1
    denovo_cds_2 = [cds for cds in cdss if cds['start'] == 1348 and cds['stop'] == 2388 and cds['strand'] == '-']
    assert len(denovo_cds_2) == 0  # de novo-predicted was filtered-out

    # test case 3: user-provided CDS +[981,1121] overlapping de novo-predicted CDS -[971,1351]
    provided_cds_3 = [cds for cds in cdss if cds['start'] == 981 and cds['stop'] == 1121 and cds['strand'] == '+']
    assert len(provided_cds_3) == 1
    denovo_cds_3 = [cds for cds in cdss if cds['start'] == 971 and cds['stop'] == 1351 and cds['strand'] == '-']
    assert len(denovo_cds_3) == 0  # de novo-predicted was filtered-out

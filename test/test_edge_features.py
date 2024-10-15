import json

from pathlib import Path
from subprocess import run

from bakta import constants as bc


def test_bakta_edge_features(tmpdir):
    # test edge lable on mock CDS sequence
    proc = run(
        [
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test',
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot',
            '--keep-contig-headers', '--complete', 'test/data/cds.fna'
        ]
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    output_path = tmpdir_path.joinpath('test.json')
    assert Path.exists(output_path)
    assert output_path.stat().st_size > 0

    with output_path.open() as fh:
        results = json.load(fh)

    for feat in results['features']:
        if(feat['sequence'] != 'dummy'):
            if('forward' in feat['sequence']):
                assert feat['strand'] == bc.STRAND_FORWARD
            elif('reverse' in feat['sequence']):
                assert feat['strand'] == bc.STRAND_REVERSE
            assert feat.get('edge', False) == ('edge' in feat['sequence'])

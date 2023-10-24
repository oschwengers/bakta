import json

from pathlib import Path
from subprocess import run


CRISPR_ARRAYS = [
    {
        'repeat_consensus': 'CGGTTTATCCCCGCTGGCGCGGGGAACACA',
        'spacers': [
            'AACCGAAACACACGATCAATCCGAATATGAG',
            'TTGGTGACAGTTTTTGTCACTGTTTTGGTGA',
            'CTAAGCATACATATCTGTTTTTAAACA'
        ],
        'repeats': 3
    }
]


def test_crispr_arrays(tmpdir):
    proc = run(
        [
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test',
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-cds', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot',
            'test/data/GCF_000008865.2.fna.gz'
        ]
    )
    assert proc.returncode == 0

    results_path = Path(tmpdir).joinpath('test.json')
    assert Path.exists(results_path)
    
    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    
    crispr_arrays = [feat for feat in results['features'] if feat['type'] == 'crispr']
    assert len(crispr_arrays) == 1

    for idx, crispr_array in enumerate(crispr_arrays):
        assert crispr_array['repeat_consensus'] == CRISPR_ARRAYS[idx]['repeat_consensus']
        assert len(crispr_array['repeats']) == CRISPR_ARRAYS[idx]['repeats']

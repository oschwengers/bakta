import pytest

import bakta.config as cfg

from bakta.features import s_orf as bu


CONTIG_1 = {
    'id': 1,
    'description': 'no sORFs',
    'nt': 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
}
CONTIG_2 = {
    'id': 2,
    'description': 'out of limits',
    'nt': 'ATGAAAAAATAGGGGATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAG'
}
CONTIG_3 = {
    'id': 3,
    'description': 'two sORFs',
    'nt': 'ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAG'
}

GENOME_1 = {
    'sequences': [CONTIG_1]
}
GENOME_2 = {
    'sequences': [CONTIG_2]
}
GENOME_3 = {
    'sequences': [CONTIG_3]
}


@pytest.mark.parametrize(
    "genome, expected",
    [
        (GENOME_1, 0),
        (GENOME_2, 0),
        (GENOME_3, 2)
    ]
)
def test_sORF(genome, expected):
    cfg.translation_table = 11
    assert len(bu.extract(genome)) == expected

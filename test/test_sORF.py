import pytest

import bakta.config as cfg

from bakta.features import s_orf as bu


contig_1 = {
    'id': 1,
    'description': 'no sORFs',
    'sequence': 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
}
contig_2 = {
    'id': 2,
    'description': 'out of limits',
    'sequence': 'ATGAAAAAATAGGGGATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAG'
}
contig_3 = {
    'id': 3,
    'description': 'two sORFs',
    'sequence': 'ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAG'
}

genome_1 = {
    'contigs': [contig_1]
}
genome_2 = {
    'contigs': [contig_2]
}
genome_3 = {
    'contigs': [contig_3]
}


@pytest.mark.parametrize(
    "genome, expected",
    [
        (genome_1, 0),
        (genome_2, 0),
        (genome_3, 2)
    ]
)
def test_sORF(genome, expected):
    cfg.translation_table = 11
    assert len(bu.extract(genome)) == expected

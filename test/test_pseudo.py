from typing import Dict, Set, Union

import bakta.constants as bc
import bakta.features.cds as feat_cds

import pytest


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
        (  # deletion
            'MKEGQFVGY/FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIH',
            'MKEGQFVGY-FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIN',
            {
                bc.PSEUDOGENE_INSERTION: set(),
                bc.PSEUDOGENE_DELETION: {28},
                bc.PSEUDOGENE_START: set(),
                bc.PSEUDOGENE_STOP: set(),
                bc.PSEUDOGENE_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_PYROLYSINE: set(),
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: True
            }
        ),
        (  # insertion
            'MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKNVIWDYPKDFIAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL\\KIIA\\AVANNYMPGVASYAG',
            'MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKDVIWDYPKDFMAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL-KIIA-PVANDYMPGVDSYAG',
            {
                bc.PSEUDOGENE_INSERTION: {430, 443},
                bc.PSEUDOGENE_DELETION: set(),
                bc.PSEUDOGENE_START: set(),
                bc.PSEUDOGENE_STOP: set(),
                bc.PSEUDOGENE_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_PYROLYSINE: set(),
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: True
            }
        ),
        (  # internal stop
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASQGSWYNF',
            {
                bc.PSEUDOGENE_INSERTION: set(),
                bc.PSEUDOGENE_DELETION: set(),
                bc.PSEUDOGENE_START: set(),
                bc.PSEUDOGENE_STOP: {277},
                bc.PSEUDOGENE_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_PYROLYSINE: set(),
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: True
            }
        ),
        (  # selenocysteine
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASUGSWYNF',
            {
                bc.PSEUDOGENE_INSERTION: set(),
                bc.PSEUDOGENE_DELETION: set(),
                bc.PSEUDOGENE_START: set(),
                bc.PSEUDOGENE_STOP: set(),
                bc.PSEUDOGENE_SELENOCYSTEINE: {277},
                bc.PSEUDOGENE_PYROLYSINE: set(),
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: False
            }
        ),
        (  # pyrolysine
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASOGSWYNF',
            {
                bc.PSEUDOGENE_INSERTION: set(),
                bc.PSEUDOGENE_DELETION: set(),
                bc.PSEUDOGENE_START: set(),
                bc.PSEUDOGENE_STOP: set(),
                bc.PSEUDOGENE_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_PYROLYSINE: {277},
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: False
            }
        )
    ]
)
def test_compare_alignments(alignment, ref_alignment, expected_result):
    # Includes test_downstream_elongation
    causes = {
        bc.PSEUDOGENE_INSERTION: set(),
        bc.PSEUDOGENE_DELETION: set(),
        bc.PSEUDOGENE_START: set(),
        bc.PSEUDOGENE_STOP: set(),
        bc.PSEUDOGENE_SELENOCYSTEINE: set(),
        bc.PSEUDOGENE_PYROLYSINE: set(),
        bc.FEATURE_END_3_PRIME: False,
        bc.FEATURE_END_5_PRIME: False
    }
    cds = {
        'start': 1,
        'contig': 'foo',
        'stop': 100,
        'strand': '+',
        'edge': False
    }
    direction = bc.FEATURE_END_5_PRIME
    feat_cds.compare_alignments(causes, alignment, ref_alignment, cds, direction)
    assert causes == expected_result


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
        (  # point mutation -> internal start codon
            'MINWRKVGMTSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
            'MINWRKVGATSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
            {
                bc.PSEUDOGENE_INSERTION: set(),
                bc.PSEUDOGENE_DELETION: set(),
                bc.PSEUDOGENE_START: {30},
                bc.PSEUDOGENE_STOP: set(),
                bc.PSEUDOGENE_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_PYROLYSINE: set(),
                bc.FEATURE_END_3_PRIME: False,
                bc.FEATURE_END_5_PRIME: True
            }
        )
    ]
)
def test_upstream_elongation(alignment, ref_alignment, expected_result):
    causes = {
        bc.PSEUDOGENE_INSERTION: set(),
        bc.PSEUDOGENE_DELETION: set(),
        bc.PSEUDOGENE_START: set(),
        bc.PSEUDOGENE_STOP: set(),
        bc.PSEUDOGENE_SELENOCYSTEINE: set(),
        bc.PSEUDOGENE_PYROLYSINE: set(),
        bc.FEATURE_END_3_PRIME: False,
        bc.FEATURE_END_5_PRIME: False
    }
    cds = {
        'start': 30,
        'stop': '',
        'strand': '',
        'contig': '',
        'rbs_motif': None
    }
    extended_positions = {'start': 1}
    elongated_edge = False
    feat_cds.upstream_elongation(causes, alignment, ref_alignment, 8, extended_positions, cds, elongated_edge)
    assert causes == expected_result


@pytest.mark.parametrize('cds, contig, expected_result', [
        (
            {
              'start': 310,  # linear fits cutoff
              'stop': 370,
              'strand': '+',
              'edge': False
            },
            {
              'sequence': 'ACGT' * 200,
              'topology': 'linear'
            },
            {
              'start': 10,
              'stop': 670,
              'strand': '+',
              'edge': False
            }
        ),
        (
            {
              'start': 100,  # linear does not fit cutoff
              'stop': 190,
              'strand': '+',
              'edge': False
            },
            {
              'sequence': 'ACGT' * 50,  # 200nt
              'topology': 'linear'
            },
            {
              'start': 0,
              'stop': 200,
              'strand': '+',
              'edge': False
            }
        ),
        (
            {
              'start': 100,  # circular does not fit cutoff
              'stop': 190,
              'strand': '+',
              'edge': True
            },
            {
              'sequence': 'ACGT' * 100,  # 400nt
              'topology': 'circular'
            },
            {
              'start': 200,
              'stop': 90,
              'strand': '+',
              'edge': True
            }
        )
    ]
)
def test_get_elongated_cds(cds, contig, expected_result):
    assert feat_cds.get_elongated_cds(cds, contig) == expected_result

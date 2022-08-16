from typing import Dict, Set, Union

import bakta.constants as bc
import bakta.features.cds as feat_cds

import pytest


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
        (  # deletion
            'MKEGQFVGY/FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIH',
            'MKEGQFVGY-FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIN',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: {28},
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_5_PRIME}
            }
        ),
        (  # insertion
            'MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKNVIWDYPKDFIAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL\\KIIA\\AVANNYMPGVASYAG',
            'MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKDVIWDYPKDFMAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL-KIIA-PVANDYMPGVDSYAG',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: {430, 443},
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_5_PRIME}
            }
        ),
        (  # internal stop
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASQGSWYNF',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),  #{277},
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: {277},
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_5_PRIME}
            }
        ),
        (  # selenocysteine
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASUGSWYNF',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: {277},
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': set()
            }
        ),
        (  # pyrolysine
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASOGSWYNF',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: {277},
                'directions': set()
            }
        )
    ]
)
def test_compare_alignments(alignment, ref_alignment, expected_result):
    # Includes test_downstream_elongation
    observations = {
        bc.PSEUDOGENE_CAUSE_INSERTION: set(),
        bc.PSEUDOGENE_CAUSE_DELETION: set(),
        bc.PSEUDOGENE_CAUSE_MUTATION: set(),
        bc.PSEUDOGENE_EFFECT_START: set(),
        bc.PSEUDOGENE_EFFECT_STOP: set(),
        bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
        bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
        'directions': set()
    }
    cds = {
        'start': 1,
        'contig': 'foo',
        'stop': 100,
        'strand': '+',
        'edge': False
    }
    feat_cds.compare_alignments(observations, alignment, ref_alignment, cds, bc.FEATURE_END_5_PRIME)
    assert observations == expected_result


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
        (  # point mutation -> internal start codon
            'MINWRKVGMTSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
            'MINWRKVGATSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
            {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: {30},
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_5_PRIME}
            }
        )
    ]
)
def test_upstream_elongation(alignment, ref_alignment, expected_result):
    observations = {
        bc.PSEUDOGENE_CAUSE_INSERTION: set(),
        bc.PSEUDOGENE_CAUSE_DELETION: set(),
        bc.PSEUDOGENE_CAUSE_MUTATION: set(),
        bc.PSEUDOGENE_EFFECT_START: set(),
        bc.PSEUDOGENE_EFFECT_STOP: set(),
        bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
        bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
        'directions': set()
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
    feat_cds.detect_pseudogenization_observations_upstream(observations, alignment, ref_alignment, 8, extended_positions, cds, elongated_edge)
    assert observations == expected_result


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

from typing import Dict, Set, Union

import bakta.constants as bc
import bakta.features.cds as feat_cds

import pytest


@pytest.mark.parametrize('alignment, ref_alignment, cds, coordinates, expected_result', [
        (  # deletion
            'MKEGQFVGY/FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIH',
            'MKEGQFVGY-FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIN',
            #            *
            {
                'sequence': 'foo',
                'start': 37,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'rbs_motif': 'AGGA',
                'edge': False
            }, {
                'upstream': -36,
                'downstream': 0
            }, {
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
                'sequence': 'foo',
                'start': 1,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: {430, 443},
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        ),
        (  # internal stop
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASQGSWYNF',
{
                'sequence': 'foo',
                'start': 1,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: {277},
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: {277},
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        ),
        (  # selenocysteine
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASUGSWYNF',
            {
                'sequence': 'foo',
                'start': 1,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
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
                'sequence': 'foo',
                'start': 1,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: {277},
                'directions': set()
            }
        ),
        (  # point mutation -> internal start codon
            'MLSIQSNRDWLSMSIFSDYSSSSEMHNNLTIDYYLALSSTKGSGITNIISIILQQAQDYDVAKIT',
            'MLSIQSNRDWLSASIFSDYSSSSEMHNNLTIDYYLALSSTKGSGITNIISIILQQAQDYDVAKIT',
            {
                'sequence': 'foo',
                'start': 40,
                'stop': 100,
                'strand': bc.STRAND_FORWARD,
                'rbs_motif': None,
                'edge': False
            }, {
                'upstream': -39,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: {40},
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        ),
        (  # deletion
            'MKEGQFVGY/FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIH',
            'MKEGQFVGY-FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIN',
            #            *
            {
                'sequence': 'foo',
                'start': 10,
                'stop': 200,
                'strand': bc.STRAND_REVERSE,
                'rbs_motif': 'AGGA',
                'edge': False
            }, {
                'upstream': -36,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: {209},
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
                'sequence': 'foo',
                'start': 1,
                'stop': 500,
                'strand': bc.STRAND_REVERSE,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: {71, 58},
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        ),
        (  # internal stop
            'MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
            'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASQGSWYNF',
{
                'sequence': 'foo',
                'start': 100,
                'stop': 500,
                'strand': bc.STRAND_REVERSE,
                'edge': False
            }, {
                'upstream': 0,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: {224},
                bc.PSEUDOGENE_EFFECT_START: set(),
                bc.PSEUDOGENE_EFFECT_STOP: {224},
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        ),
        (  # point mutation -> internal start codon
            'MLSIQSNRDWLSMSIFSDYSSSSEMHNNLTIDYYLALSSTKGSGITNIISIILQQAQDYDVAKIT',
            'MLSIQSNRDWLSASIFSDYSSSSEMHNNLTIDYYLALSSTKGSGITNIISIILQQAQDYDVAKIT',
            {
                'sequence': 'foo',
                'start': 40,
                'stop': 100,
                'strand': bc.STRAND_REVERSE,
                'rbs_motif': None,
                'edge': False
            }, {
                'upstream': -39,
                'downstream': 0
            }, {
                bc.PSEUDOGENE_CAUSE_INSERTION: set(),
                bc.PSEUDOGENE_CAUSE_DELETION: set(),
                bc.PSEUDOGENE_CAUSE_MUTATION: set(),
                bc.PSEUDOGENE_EFFECT_START: {100},
                bc.PSEUDOGENE_EFFECT_STOP: set(),
                bc.PSEUDOGENE_EXCEPTION_SELENOCYSTEINE: set(),
                bc.PSEUDOGENE_EXCEPTION_PYROLYSINE: set(),
                'directions': {bc.FEATURE_END_3_PRIME}
            }
        )
    ]
)
def test_compare_alignments(alignment, ref_alignment, cds, coordinates, expected_result):
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

    feat_cds.compare_alignments(observations, alignment, ref_alignment, cds, coordinates, edge=False)
    assert observations == expected_result


@pytest.mark.parametrize('cds, sequence, expected_result', [
        (
            {
              'start': 310,  # linear fits cutoff
              'stop': 370,
              'strand': bc.STRAND_FORWARD,
              'edge': False
            },
            {
              'nt': 'ACGT' * 200,
              'topology': 'linear'
            },
            {
              'start': 10,
              'stop': 670,
              'strand': bc.STRAND_FORWARD,
              'edge': False,
              'elongation_upstream': 300,
              'elongation_downstream': 300
            }
        ),
        (
            {
              'start': 100,  # linear does not fit cutoff
              'stop': 190,
              'strand': bc.STRAND_FORWARD,
              'edge': False
            },
            {
              'nt': 'ACGT' * 50,  # 200nt
              'topology': 'linear'
            },
            {
              'start': 1,
              'stop': 200,
              'strand': bc.STRAND_FORWARD,
              'edge': False,
              'elongation_upstream': 100,
              'elongation_downstream': 10
            }
        ),
        (
            {
              'start': 100,  # circular does not fit cutoff
              'stop': 190,
              'strand': bc.STRAND_FORWARD,
              'edge': True,
              'elongation_upstream': 300,
              'elongation_downstream': 300
            },
            {
              'nt': 'ACGT' * 100,  # 400nt
              'topology': 'circular'
            },
            {
              'start': 200,
              'stop': 90,
              'strand': bc.STRAND_FORWARD,
              'edge': True,
              'elongation_upstream': 300,
              'elongation_downstream': 300
            }
        )
    ]
)
def test_get_elongated_cds(cds, sequence, expected_result):
    assert feat_cds.get_elongated_cds(cds, sequence, offset=300) == expected_result

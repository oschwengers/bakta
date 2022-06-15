from typing import Dict, Set, Union

import bakta.features.cds as feat_cds

import pytest


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
    ('MKEGQFVGY/FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIH',
     'MKEGQFVGY-FKMKEQRKIPLTHIMIIGAFIFAFLQVVLLASLVHAVNVNNEIQEGLFQSGRIMVESLQHILSVQTGIN',
     {'insertion': set(),
      'deletion': {28},
      'start': set(),
      'stop': set(),
      'selenocysteine': set(),
      'pyrolysine': set(),
      3: False,
      5: True}
     ),  # deletion
    ('MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKNVIWDYPKDFIAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL\\KIIA\\AVANNYMPGVASYAG',
     'MTQRPWSKLQREIYDLLTPTINLQIHCTRYPMRSQNGGSTDLPRYWITLDKDVIWDYPKDFMAGNGGVRNFHGETCWYPYLTDICSISDLLREYIDTPKAELLTKQFTSDKWGLVNILRAADRRIGMRRLDQLRRKTHNIAAL-KIIA-PVANDYMPGVDSYAG',
     {'insertion': {430, 443},
      'deletion': set(),
      'start': set(),
      'stop': set(),
      'selenocysteine': set(),
      'pyrolysine': set(),
      3: False,
      5: True}
     ),  # insertion
    ('MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
     'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASQGSWYNF',
     {'insertion': set(),
      'deletion': set(),
      'start': set(),
      'stop': {277},
      'selenocysteine': set(),
      'pyrolysine': set(),
      3: False,
      5: True}
     ),  # internal stop
    ('MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
     'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASUGSWYNF',
     {'insertion': set(),
      'deletion': set(),
      'start': set(),
      'stop': set(),
      'selenocysteine': {277},
      'pyrolysine': set(),
      3: False,
      5: False}
     ),  # selenocysteine
    ('MSLYIKLILSIVREISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNDAS*GSWYNF',
     'MPLYIKLILSIVRRISVNTICSLIVVVALSLLSFSSVAKTITAVGSTINSTEKEISLQAEKQGKSYKILGAFFKNRVYMIAKLTPVSKNNASOGSWYNF',
     {'insertion': set(),
      'deletion': set(),
      'start': set(),
      'stop': set(),
      'selenocysteine': set(),
      'pyrolysine': {277},
      3: False,
      5: False}
     ),  # pyrolysine
])
def test_compare_alignments(alignment, ref_alignment, expected_result):
    # Includes test_downstream_elongation
    cause: Dict[Union[str, int], Union[Set[int], bool]] = {'insertion': set(),
                                                           'deletion': set(),
                                                           'start': set(),
                                                           'stop': set(),
                                                           'selenocysteine': set(),
                                                           'pyrolysine': set(),
                                                           3: False,
                                                           5: False}
    cds: Dict = {'start': 1,
                 'contig': '',
                 'stop': '',
                 'strand': ''}
    direction: int = 5

    assert feat_cds.compare_alignments(cause, alignment, ref_alignment, cds, direction) == expected_result


@pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
    ('MINWRKVGMTSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
     'MINWRKVGATSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
     {'insertion': set(),
      'deletion': set(),
      'start': {30},
      'stop': set(),
      'selenocysteine': set(),
      'pyrolysine': set(),
      3: False,
      5: True}
     ),  # point mutation -> internal start codon
    # ('VINWRKVGMTSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
    #  'MINWRKVGMTSSHHGPYDQGYTRATMAHTKRSDLARASGPHKVRRSPDWSLQLDSMKSESLVIVDQNATVNTFPGLVHTARHTMGVGCKRSR',
    #  {'insertion': set(),
    #   'deletion': set(),
    #   'start': {5},  # TODO check
    #   'stop': set(),
    #   'selenocysteine': set(),
    #   'pyrolysine': set(),
    #   3: False,
    #   5: True}
    #  )  # point mutation -> loss of original start codon
])
def test_upstream_elongation(alignment, ref_alignment, expected_result):
    cause: Dict[Union[str, int], Union[Set[int], bool]] = {'insertion': set(),
                                                           'deletion': set(),
                                                           'start': set(),
                                                           'stop': set(),
                                                           'selenocysteine': set(),
                                                           'pyrolysine': set(),
                                                           3: False,
                                                           5: False}
    cds: Dict = {'start': 30,
                 'stop': '',
                 'strand': '',
                 'contig': '',
                 'rbs_motif': None}
    extended_positions: Dict = {'start': 1}
    elongated_edge: bool = False

    assert feat_cds.upstream_elongation(cause, alignment, ref_alignment, 8, extended_positions, cds, elongated_edge) == expected_result


# @pytest.mark.parametrize('alignment, ref_alignment, expected_result', [
#     ('MAVKRDMPEESKNSKVVKKEHFSIVFPDDIKVPKSEKELEAEKAENKSEHD',
#      'MAVKRDMPEESKNSKVVKKEHFSIVFPDDIKEPSDKDEQKKKTIDTKKDND',
#      {'insertion': set(),
#       'deletion': set(),
#       'start': set(),
#       'stop': {156},
#       'selenocysteine': set(),
#       'pyrolysine': set(),
#       3: True,
#       5: False}
#      ),  # loss of stop codon
# ])
# def test_loss_of_stop_codon(alignment, ref_alignment, expected_result):
#     cds: Dict = {'start': 1,
#                  'stop': 165,
#                  'strand': '',
#                  'contig': '',
#                  'aa': 'MAVKRDMPEESKNSKVVKKEHFSIVFPDDIKVPKSEKELEAEKAENKSEHDKTN',
#                  'rbs_motif': None,
#                  'pseudo-candidate': {'cluster_sequence': 'MAVKRDMPEESKNSKVVKKEHFSIVFPDDIKEPSDKDEQKKKTIDTKKDND'}}
#     extended_positions: Dict = {'start': 1,
#                                 'len_psc': 51}
#
#     assert feat_cds.pseudogenization_cause(alignment, ref_alignment, 1, 153, extended_positions, cds) == expected_result

# EOF

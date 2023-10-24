import pytest

from bakta.features import signal_peptides as bsp


ONE_TWENTYONE_FWD = {'start': 1, 'stop': 21, 'strand': '+'}
ONE_TWENTYONE_REV = {'start': 1, 'stop': 21, 'strand': '-'}
FOUR_FIFTEEN_FWD = {'start': 4, 'stop': 15, 'strand': '+'}
FOUR_FIFTEEN_REV = {'start': 4, 'stop': 15, 'strand': '-'}


@pytest.mark.parametrize(
    "orf, start_aa, stop_aa, expected",
    [
        # Forward strand
        (ONE_TWENTYONE_FWD, 1, 3, (1,9)),   # ORF spans whole sequence length, signal peptide starts at AA 1
        (ONE_TWENTYONE_FWD, 2, 7, (4,21)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (FOUR_FIFTEEN_FWD, 1, 3, (4,12)),   # ORF does not span whole sequence length, signal peptide starts at AA 1
        (FOUR_FIFTEEN_FWD, 2, 4, (7,15)),   # ORF does not span whole sequence length, signal peptide stops at last AA
        # Reverse strand
        (ONE_TWENTYONE_REV, 1, 3, (13,21)), # ORF spans whole sequence length, signal peptide starts at AA 1
        (ONE_TWENTYONE_REV, 3, 7, (1,15)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (FOUR_FIFTEEN_REV, 1, 2, (10,15)),  # ORF does not span whole sequence length, signal peptide starts at AA 1
        (FOUR_FIFTEEN_REV, 2, 3, (7,12))    # ORF does not span whole sequence length, signal peptide in the middle of ORF
    ]
)
def test_start_stop(orf, start_aa, stop_aa, expected):
    assert bsp.orf_nt_start_stop(orf, start_aa, stop_aa) == expected
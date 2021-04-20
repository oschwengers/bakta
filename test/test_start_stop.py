import pytest

from bakta.features import signal_peptides as bu

oneTwentyone_fwd = {'start': 1, 'stop': 21, 'strand': '+'}
oneTwentyone_rev = {'start': 1, 'stop': 21, 'strand': '-'}
fourFifteen_fwd = {'start': 4, 'stop': 15, 'strand': '+'}
fourFifteen_rev = {'start': 4, 'stop': 15, 'strand': '-'}

@pytest.mark.parametrize(
    "orf, start, stop, expected",
    [
        # Forward strand
        (oneTwentyone_fwd, 1, 3, (1,9)),   # ORF spans whole sequence length, signal peptide starts at AA 1
        (oneTwentyone_fwd, 2, 7, (4,21)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (fourFifteen_fwd, 1, 3, (4,12)),   # ORF does not span whole sequence length, signal peptide starts at AA 1
        (fourFifteen_fwd, 2, 4, (7,15)),   # ORF does not span whole sequence length, signal peptide stops at last AA
        # Reverse strand
        (oneTwentyone_rev, 1, 3, (13,21)), # ORF spans whole sequence length, signal peptide starts at AA 1
        (oneTwentyone_rev, 3, 7, (1,15)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (fourFifteen_rev, 1, 2, (10,15)),  # ORF does not span whole sequence length, signal peptide starts at AA 1
        (fourFifteen_rev, 2, 3, (7,12))    # ORF does not span whole sequence length, signal peptide in the middle of ORF
    ]
)
def test_start_stop(orf, start, stop, expected):
    assert bu.start_stop_orf(orf, start, stop) == expected

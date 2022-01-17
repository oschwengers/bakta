import pytest

from bakta.features import signal_peptides as bsp


one_twentyone_fwd = {'start': 1, 'stop': 21, 'strand': '+'}
one_twentyone_rev = {'start': 1, 'stop': 21, 'strand': '-'}
four_fifteen_fwd = {'start': 4, 'stop': 15, 'strand': '+'}
four_fifteen_rev = {'start': 4, 'stop': 15, 'strand': '-'}


@pytest.mark.parametrize(
    "orf, start_aa, stop_aa, expected",
    [
        # Forward strand
        (one_twentyone_fwd, 1, 3, (1,9)),   # ORF spans whole sequence length, signal peptide starts at AA 1
        (one_twentyone_fwd, 2, 7, (4,21)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (four_fifteen_fwd, 1, 3, (4,12)),   # ORF does not span whole sequence length, signal peptide starts at AA 1
        (four_fifteen_fwd, 2, 4, (7,15)),   # ORF does not span whole sequence length, signal peptide stops at last AA
        # Reverse strand
        (one_twentyone_rev, 1, 3, (13,21)), # ORF spans whole sequence length, signal peptide starts at AA 1
        (one_twentyone_rev, 3, 7, (1,15)),  # ORF spans whole sequence length, signal peptide stops at last AA
        (four_fifteen_rev, 1, 2, (10,15)),  # ORF does not span whole sequence length, signal peptide starts at AA 1
        (four_fifteen_rev, 2, 3, (7,12))    # ORF does not span whole sequence length, signal peptide in the middle of ORF
    ]
)
def test_start_stop(orf, start_aa, stop_aa, expected):
    assert bsp.orf_nt_start_stop(orf, start_aa, stop_aa) == expected
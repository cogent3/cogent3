"""Tests for Clustal sequence format writer."""

from unittest import TestCase

import cogent3
from cogent3.format.clustal import clustal_from_alignment


class ClustalTests(TestCase):
    """Tests for Clustal writer."""

    def setUp(self):
        """Setup for Clustal tests."""
        self.unaligned_dict = {
            "1st": "AAA",
            "2nd": "CCCC",
            "3rd": "GGGG",
            "4th": "UUUU",
        }
        self.alignment_dict = {
            "1st": "AAAA",
            "2nd": "CCCC",
            "3rd": "GGGG",
            "4th": "UUUU",
        }
        # create alignment change order.
        self.alignment_object = cogent3.make_aligned_seqs(
            self.alignment_dict,
            moltype="text",
        )
        self.alignment_order = ["2nd", "4th", "3rd", "1st"]

        self.clustal_with_label = """CLUSTAL

1st    AAAA
2nd    CCCC
3rd    GGGG
4th    UUUU
"""
        self.clustal_with_label_lw2 = """CLUSTAL

1st    AA
2nd    CC
3rd    GG
4th    UU

1st    AA
2nd    CC
3rd    GG
4th    UU
"""

        self.clustal_with_label_reordered = """CLUSTAL

2nd    CCCC
4th    UUUU
3rd    GGGG
1st    AAAA
"""

        self.clustal_with_label_lw2_reordered = """CLUSTAL

2nd    CC
4th    UU
3rd    GG
1st    AA

2nd    CC
4th    UU
3rd    GG
1st    AA
"""

    def test_clustal_from_alignment_unaligned(self):
        """should raise error with unaligned seqs."""
        self.assertRaises(ValueError, clustal_from_alignment, self.unaligned_dict)

    def test_clustal_from_alignment(self):
        """should return correct clustal string."""
        assert clustal_from_alignment({}) == ""
        assert clustal_from_alignment(self.alignment_dict) == self.clustal_with_label
        assert (
            clustal_from_alignment(self.alignment_dict, wrap=2)
            == self.clustal_with_label_lw2
        )

    def test_clustal_from_alignment_reordered(self):
        """should return correct clustal string."""
        assert (
            clustal_from_alignment(self.alignment_dict, order=self.alignment_order)
            == self.clustal_with_label_reordered
        )
        assert (
            clustal_from_alignment(
                self.alignment_dict,
                order=self.alignment_order,
                wrap=2,
            )
            == self.clustal_with_label_lw2_reordered
        )

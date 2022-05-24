#!/usr/bin/env python
"""Tests for Clustal sequence format writer.
"""
from unittest import TestCase, main

from cogent3.core.alignment import Alignment
from cogent3.format.clustal import clustal_from_alignment


__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"


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
        self.alignment_object = Alignment(self.alignment_dict)
        self.alignment_order = ["2nd", "4th", "3rd", "1st"]
        self.alignment_object.RowOrder = self.alignment_order

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
        self.assertEqual(clustal_from_alignment({}), "")
        self.assertEqual(
            clustal_from_alignment(self.alignment_dict), self.clustal_with_label
        )
        self.assertEqual(
            clustal_from_alignment(self.alignment_dict, wrap=2),
            self.clustal_with_label_lw2,
        )

    def test_clustal_from_alignment_reordered(self):
        """should return correct clustal string."""
        self.assertEqual(
            clustal_from_alignment(self.alignment_object),
            self.clustal_with_label_reordered,
        )
        self.assertEqual(
            clustal_from_alignment(self.alignment_object, wrap=2),
            self.clustal_with_label_lw2_reordered,
        )


if __name__ == "__main__":
    main()

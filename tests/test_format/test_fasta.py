#!/usr/bin/env python
"""Tests for FASTA sequence format writer.
"""
from unittest import TestCase, main

from cogent3.core.alignment import Alignment
from cogent3.core.info import Info
from cogent3.core.sequence import Sequence
from cogent3.format.fasta import alignment_to_fasta


__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Gavin Huttley", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"


class FastaTests(TestCase):
    """Tests for Fasta writer."""

    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ["AAAA", "CCCC", "gggg", "uuuu"]
        self.labels = ["1st", "2nd", "3rd", "4th"]
        self.infos = ["Dog", "Cat", "Mouse", "Rat"]
        self.sequences_with_labels = list(map(Sequence, self.strings))
        self.sequences_with_names = list(map(Sequence, self.strings))
        for l, sl, sn in zip(
            self.labels, self.sequences_with_labels, self.sequences_with_names
        ):
            sl.label = l
            sn.name = l
        self.fasta_no_label = ">0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu\n"
        self.fasta_with_label = ">1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU\n"
        self.fasta_with_label_lw2 = (
            ">1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU\n"
        )
        self.alignment_dict = {
            "1st": "AAAA",
            "2nd": "CCCC",
            "3rd": "GGGG",
            "4th": "UUUU",
        }
        self.alignment_object = Alignment(self.alignment_dict)
        for label, info in zip(self.labels, self.infos):
            self.alignment_object.named_seqs[label].info = Info(species=info)
        self.fasta_with_label_species = (
            ">1st:Dog\nAAAA\n>2nd:Cat\nCCCC\n>3rd:Mouse\nGGGG\n>4th:Rat\nUUUU\n"
        )
        self.alignment_object.RowOrder = ["1st", "2nd", "3rd", "4th"]

    def test_alignment_to_fasta(self):
        """should return correct fasta string."""
        self.assertEqual(alignment_to_fasta({}), "")
        self.assertEqual(alignment_to_fasta(self.alignment_dict), self.fasta_with_label)
        self.assertEqual(
            alignment_to_fasta(self.alignment_dict, block_size=2),
            self.fasta_with_label_lw2,
        )


if __name__ == "__main__":
    main()

"""Tests for FASTA sequence format writer."""

from unittest import TestCase

import cogent3
from cogent3.core.info import Info
from cogent3.format.fasta import seqs_to_fasta


class FastaTests(TestCase):
    """Tests for Fasta writer."""

    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ["AAAA", "CCCC", "gggg", "uuuu"]
        self.labels = ["1st", "2nd", "3rd", "4th"]
        self.infos = ["Dog", "Cat", "Mouse", "Rat"]
        ASCII = cogent3.get_moltype("text")
        self.sequences_with_labels = [ASCII.make_seq(seq=seq) for seq in self.strings]
        self.sequences_with_names = [ASCII.make_seq(seq=seq) for seq in self.strings]
        for l, sl, sn in zip(
            self.labels,
            self.sequences_with_labels,
            self.sequences_with_names,
            strict=False,
        ):
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
        self.alignment_object = cogent3.make_aligned_seqs(
            self.alignment_dict,
            moltype="text",
        )
        for label, info in zip(self.labels, self.infos, strict=False):
            self.alignment_object.named_seqs[label].info = Info(species=info)
        self.fasta_with_label_species = (
            ">1st:Dog\nAAAA\n>2nd:Cat\nCCCC\n>3rd:Mouse\nGGGG\n>4th:Rat\nUUUU\n"
        )
        self.alignment_object.RowOrder = ["1st", "2nd", "3rd", "4th"]

    def test_alignment_to_fasta(self):
        """should return correct fasta string."""
        self.assertEqual(seqs_to_fasta({}), "")
        self.assertEqual(seqs_to_fasta(self.alignment_dict), self.fasta_with_label)
        self.assertEqual(
            seqs_to_fasta(self.alignment_dict, block_size=2),
            self.fasta_with_label_lw2,
        )

#!/usr/bin/env python
"""Tests for FASTA sequence format writer.
"""
from unittest import TestCase, main
from cogent3.format.fasta import fasta_from_sequences, fasta_from_alignment
from cogent3.core.alignment import Alignment
from cogent3.core.sequence import Sequence
from cogent3.core.info import Info

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"


class FastaTests(TestCase):
    """Tests for Fasta writer.
    """

    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ['AAAA', 'CCCC', 'gggg', 'uuuu']
        self.labels = ['1st', '2nd', '3rd', '4th']
        self.infos = ["Dog", "Cat", "Mouse", "Rat"]
        self.sequences_with_labels = list(map(Sequence, self.strings))
        self.sequences_with_names = list(map(Sequence, self.strings))
        for l, sl, sn in zip(self.labels, self.sequences_with_labels,
                           self.sequences_with_names):
            sl.label = l
            sn.name = l
        self.fasta_no_label = '>0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu'
        self.fasta_with_label =\
        '>1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU'
        self.fasta_with_label_lw2 =\
            '>1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU'
        self.alignment_dict = {'1st': 'AAAA', '2nd': 'CCCC', '3rd': 'GGGG',
                               '4th': 'UUUU'}
        self.alignment_object = Alignment(self.alignment_dict)
        for label, info in zip(self.labels, self.infos):
            self.alignment_object.named_seqs[label].info = Info(species=info)
        self.fasta_with_label_species =\
        '>1st:Dog\nAAAA\n>2nd:Cat\nCCCC\n>3rd:Mouse\nGGGG\n>4th:Rat\nUUUU'
        self.alignment_object.RowOrder = ['1st', '2nd', '3rd', '4th']

    def test_fastaFromSequence(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_sequences(''), '')
        self.assertEqual(fasta_from_sequences(self.strings),
                         self.fasta_no_label)
        self.assertEqual(fasta_from_sequences(self.sequences_with_labels),
                         self.fasta_with_label)
        self.assertEqual(fasta_from_sequences(self.sequences_with_names),
                         self.fasta_with_label)
        make_seqlabel = lambda seq: "%s:%s" % (seq.name, seq.info.species)
        seqs = [self.alignment_object.named_seqs[label]
            for label in self.labels]
        self.assertEqual(fasta_from_sequences(seqs,
                                              make_seqlabel=make_seqlabel), self.fasta_with_label_species)

    def test_fasta_from_alignment(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_alignment({}), '')
        self.assertEqual(fasta_from_alignment(self.alignment_dict),
                         self.fasta_with_label)
        self.assertEqual(fasta_from_alignment(self.alignment_dict,
                                              line_wrap=2), self.fasta_with_label_lw2)
        self.assertEqual(fasta_from_alignment(self.alignment_object),
                         self.fasta_with_label)

if __name__ == "__main__":
    main()

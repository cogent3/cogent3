#!/usr/bin/env python
"""Tests for FASTA sequence format writer.
"""
from cogent.util.unit_test import TestCase, main
from cogent.format.fasta import fasta_from_sequences, fasta_from_alignment
from cogent.core.alignment import Alignment
from cogent.core.sequence import Sequence
from cogent.core.info import Info

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class FastaTests(TestCase):
    """Tests for Fasta writer.
    """
    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ['AAAA','CCCC','gggg','uuuu']
        self.labels = ['1st','2nd','3rd','4th']
        self.infos = ["Dog", "Cat", "Mouse", "Rat"]
        self.sequences_with_labels = map(Sequence, self.strings)
        self.sequences_with_names = map(Sequence, self.strings)
        for l,sl,sn in zip(self.labels,self.sequences_with_labels,\
            self.sequences_with_names):
            sl.Label = l
            sn.Name = l
        self.fasta_no_label='>0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu'
        self.fasta_with_label=\
                '>1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU'
        self.fasta_with_label_lw2=\
            '>1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU'
        self.alignment_dict = {'1st':'AAAA','2nd':'CCCC','3rd':'GGGG',
            '4th':'UUUU'}
        self.alignment_object = Alignment(self.alignment_dict)
        for label, info in zip(self.labels, self.infos):
            self.alignment_object.NamedSeqs[label].Info = Info(species=info)
        self.fasta_with_label_species=\
              '>1st:Dog\nAAAA\n>2nd:Cat\nCCCC\n>3rd:Mouse\nGGGG\n>4th:Rat\nUUUU'
        self.alignment_object.RowOrder = ['1st','2nd','3rd','4th']
    
    def test_fastaFromSequence(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_sequences(''),'')
        self.assertEqual(fasta_from_sequences(self.strings),\
            self.fasta_no_label)
        self.assertEqual(fasta_from_sequences(self.sequences_with_labels),\
            self.fasta_with_label)
        self.assertEqual(fasta_from_sequences(self.sequences_with_names),\
            self.fasta_with_label)
        make_seqlabel = lambda seq: "%s:%s" % (seq.Name, seq.Info.species)
        seqs = [self.alignment_object.NamedSeqs[label] for label in self.labels]
        self.assertEqual(fasta_from_sequences(seqs, 
                        make_seqlabel=make_seqlabel), self.fasta_with_label_species)
    
    def test_fasta_from_alignment(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_alignment({}),'')
        self.assertEqual(fasta_from_alignment(self.alignment_dict),\
            self.fasta_with_label)
        self.assertEqual(fasta_from_alignment(self.alignment_dict,
                line_wrap=2),self.fasta_with_label_lw2)
        self.assertEqual(fasta_from_alignment(self.alignment_object),\
            self.fasta_with_label)

if __name__ == "__main__":
    main()

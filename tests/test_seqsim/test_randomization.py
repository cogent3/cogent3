#!/usr/bin/env python
"""Unit tests for the microarray module, dealing with fake expression data."""
from cogent.util.unit_test import TestCase, main
from cogent.seqsim.randomization import shuffle_range, shuffle_between, \
    shuffle_except_indices, shuffle_except

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class randomization_tests(TestCase):
    """Tests of the top-level functionality"""
    def setUp(self):
        """Make some standard objects to randomize"""
        self.numbers = list('123')
        self.letters = list('abcdef')
        self.to_test = self.numbers + 2*self.letters + self.numbers
        
    def test_shuffle_range(self):
        """shuffle_range should shuffle only inside range"""
        shuffle_range(self.to_test, 3, -3)
        self.assertEqual(self.to_test[:3],self.numbers)
        self.assertEqual(self.to_test[-3:], self.numbers)
        self.assertNotEqual(self.to_test[3:-3], 2*self.letters)
        self.assertEqualItems(self.to_test[3:-3], 2*self.letters)
        #this time, start is negative and end is positive
        shuffle_range(self.to_test, -15, 15)
        self.assertEqual(self.to_test[:3],self.numbers)
        self.assertEqual(self.to_test[-3:], self.numbers)
        self.assertNotEqual(self.to_test[3:-3], 2*self.letters)
        self.assertEqualItems(self.to_test[3:-3], 2*self.letters)

    def test_shuffle_between(self):
        """shuffle_between should shuffle between specified chars"""
        shuffle_peptides = shuffle_between('KR')
        seq1 = 'AGHCDSGAHF' #each 10 chars long 
        seq2 = 'PLMIDNYHGT'
        protein = seq1 + 'K' + seq2
        result = shuffle_peptides(protein)
        self.assertEqual(result[10], 'K')
        self.assertNotEqual(result[:10], seq1)
        self.assertEqualItems(result[:10], seq1)
        self.assertNotEqual(result[11:], seq2)
        self.assertEqualItems(result[11:], seq2)

    def test_shuffle_except_indices(self):
        """shuffle_except_indices should shuffle all except specified indices"""
        seq1 = 'AGHCDSGAHF' #each 10 chars long 
        seq2 = 'PLMIDNYHGT'
        protein = seq1 + 'K' + seq2
        result = list(protein)
        shuffle_except_indices(result, [10])
        self.assertEqual(result[10], 'K')
        self.assertNotEqual(''.join(result), protein)
        self.assertEqualItems(''.join(result), protein)
        self.assertNotEqualItems(''.join(result[:10]), seq1)
        
    def test_shuffle_except(self):
        """shuffle_except_indices should shuffle all except specified indices"""
        seq1 = 'AGHCDSGAHF' #each 10 chars long 
        seq2 = 'PLMIDNYHGT'
        protein = seq1 + 'K' + seq2
        prot = protein
        se = shuffle_except('K')
        result = se(prot)
        self.assertEqual(result[10], 'K')
        self.assertNotEqual(''.join(result), protein)
        self.assertEqualItems(''.join(result), protein)
        self.assertNotEqualItems(''.join(result[:10]), seq1)
 

if __name__ == '__main__':
    main()

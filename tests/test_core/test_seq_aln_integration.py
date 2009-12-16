#!/usr/bin/env python

from __future__ import division
from numpy import array, transpose, alltrue
from cogent.util.unit_test import TestCase, main
from cogent.core.moltype import RNA
from cogent.core.sequence import RnaSequence, Sequence, ModelSequence
from cogent.core.alignment import Alignment, DenseAlignment

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class AllTests(TestCase):

    def setUp(self):
        """setUp method for all tests"""
        # named sequences
        self.rna1 = RnaSequence('UCAGGG', Name='rna1')
        self.rna2 = RnaSequence('YCU-RG', Name='rna2')
        self.rna3 = RnaSequence('CAA-NR', Name='rna3')
        self.model1 = ModelSequence('UCAGGG', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        self.model2 = ModelSequence('YCU-RG', Name='rna2',\
            Alphabet=RNA.Alphabets.DegenGapped)
        self.model3 = ModelSequence('CAA-NR', Name='rna3',\
            Alphabet=RNA.Alphabets.DegenGapped)

        self.aln = Alignment([self.rna1, self.rna2, self.rna3], MolType=RNA)
        self.da = DenseAlignment([self.model1, self.model2, self.model3],\
            MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)

        # seqs no name
        self.nn_rna1 = RnaSequence('UCAGGG')
        self.nn_rna2 = RnaSequence('YCU-RG')
        self.nn_rna3 = RnaSequence('CAA-NR')

        self.nn_model1 = ModelSequence('UCAGGG',\
            Alphabet=RNA.Alphabets.DegenGapped)
        self.nn_model2 = ModelSequence('YCU-RG',\
            Alphabet=RNA.Alphabets.DegenGapped)
        self.nn_model3 = ModelSequence('CAA-NR',\
            Alphabet=RNA.Alphabets.DegenGapped)

        self.nn_aln = Alignment([self.nn_rna1, self.nn_rna2, self.nn_rna3],\
            MolType=RNA)
        self.nn_da = DenseAlignment([self.nn_model1, self.nn_model2,\
            self.nn_model3], MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)

    def test_printing_named_seqs(self):
        """Printing named seqs should work the same on Aln and DenseAln"""
        #Note: the newline trailing each sequence is intentional, because
        #we want each FASTA-format record to be separated.
        exp_lines_general = ['>rna1','UCAGGG','>rna2','YCU-RG','>rna3','CAA-NR']
        self.assertEqual(str(self.aln), '\n'.join(exp_lines_general) + '\n')
        self.assertEqual(str(self.da), '\n'.join(exp_lines_general) + '\n')

    def test_printing_unnamed_seqs(self):
        """Printing unnamed sequences should work the same on Aln and DenseAln
        """
        exp_lines_gen = ['>seq_0','UCAGGG','>seq_1','YCU-RG','>seq_2','CAA-NR\n']
        self.assertEqual(str(self.nn_aln),'\n'.join(exp_lines_gen))
        self.assertEqual(str(self.nn_da),'\n'.join(exp_lines_gen))

    def test_DenseAlignment_without_moltype(self):
        """Expect MolType to be picked up from the sequences."""

        m1 = ModelSequence('UCAG',Alphabet=RNA.Alphabets.DegenGapped,\
            Name='rna1')
        m2 = ModelSequence('CCCR',Alphabet=RNA.Alphabets.DegenGapped,\
            Name='rna2')
        da = DenseAlignment([m1, m2])
        exp_lines = ['>rna1','UCAG','>rna2','CCCR']
        self.assertEqual(str(da), '\n'.join(exp_lines) + '\n')

    def test_names(self):
        # Should both alignments handle names the same way?
        self.assertEqual(self.aln.Names, ['rna1','rna2','rna3'])
        self.assertEqual(self.da.Names, ['rna1','rna2','rna3'])
        # On unnamed sequences the behavior is now the same.
        self.assertEqual(self.nn_aln.Names, ['seq_0','seq_1','seq_2'])
        self.assertEqual(self.nn_da.Names, ['seq_0','seq_1','seq_2'])
    
    def test_seqFreqs(self):
        """seqFreqs should work the same on Alignment and DenseAlignment"""
        # Used alphabet: ('U', 'C', 'A', 'G', '-', 'B', 'D', 'H',\
        # 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y')
        exp = [[1,1,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0],\
            [1,1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0],\
            [0,1,2,0,1,0,0,0,0,0,1,0,1,0,0,0,0]]
        # This works
        self.assertEqual(self.da.getSeqFreqs().Data, exp)
        # This used to raise an error, but now works
        self.assertEqual(self.aln.getSeqFreqs().Data, exp)

    def test_subset_positions_DenseAlignment(self):
        model1 = ModelSequence('UCG', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        model2 = ModelSequence('YCG', Name='rna2',\
            Alphabet=RNA.Alphabets.DegenGapped)
        model3 = ModelSequence('CAR', Name='rna3',\
            Alphabet=RNA.Alphabets.DegenGapped)
        sub_da = DenseAlignment([model1, model2, model3],\
            MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)

        full_data = array([[0,1,2,3,3,3],[15,1,0,4,12,3],[1,2,2,4,10,12]])
        sub_data = array([[0,1,3],[15,1,3],[1,2,12]])
        
        # First check some data
        self.assertEqual(self.da.ArraySeqs, full_data)
        self.assertEqual(self.da.ArrayPositions, transpose(full_data))
        self.assertEqual(sub_da.ArraySeqs, sub_data)
        self.assertEqual(sub_da.ArrayPositions, transpose(sub_data))
        
        obs_sub_da_TP = self.da.takePositions([0,1,5])
        obs_sub_da_SA = self.da.getSubAlignment(pos=[0,1,5])
        
        # When using the getSubAlignment method the data is right 
        self.assertEqual(obs_sub_da_SA, sub_da)
        self.failIfEqual(obs_sub_da_SA, self.da)
        self.assertEqual(obs_sub_da_SA.ArraySeqs, sub_data)
        self.assertEqual(obs_sub_da_SA.ArrayPositions, transpose(sub_data))

        # For the takePositions method: Why does this work
        self.assertEqual(obs_sub_da_TP, sub_da)
        self.failIfEqual(obs_sub_da_TP, self.da)
        # If the data doesn't match?
        self.assertEqual(obs_sub_da_TP.ArraySeqs, sub_data)
        self.assertEqual(obs_sub_da_TP.ArrayPositions, transpose(sub_data))
        # Shouldn't the __eq__ method check the data at least?
        
    def test_subset_positions_Alignment(self):
        rna1 = RnaSequence('UCG', Name='rna1')
        rna2 = RnaSequence('YCG', Name='rna2')
        rna3 = RnaSequence('CAR', Name='rna3')
        
        sub_aln = Alignment([rna1, rna2, rna3], MolType=RNA)

        obs_sub_aln = self.aln.takePositions([0,1,5])
        self.assertEqual(obs_sub_aln, sub_aln)
        self.failIfEqual(obs_sub_aln, self.aln)
        # string representations should be the same. This fails right
        # now, because sequence order is not maintained. See separate test.
        self.assertEqual(str(obs_sub_aln), str(sub_aln))
        
    def test_takePositions_sequence_order(self):
        """Alignment takePositions should maintain seq order"""
        #This works        
        self.assertEqual(self.da.Names,['rna1','rna2','rna3'])
        sub_da = self.da.getSubAlignment(pos=[0,1,5])
        self.assertEqual(sub_da.Names,['rna1','rna2','rna3'])
        # seq order not maintained in Alignment
        self.assertEqual(self.aln.Names,['rna1','rna2','rna3']) 
        sub_aln = self.aln.takePositions([0,1,5])
        self.assertEqual(sub_aln.Names,['rna1','rna2','rna3'])

    def test_subset_seqs_Alignment(self):
        rna1 = RnaSequence('UCG', Name='rna1')
        rna2 = RnaSequence('YCG', Name='rna2')
        rna3 = RnaSequence('CAR', Name='rna3')

        sub_aln = Alignment([rna2, rna3], MolType=RNA)
        aln = Alignment([rna1, rna2, rna3], MolType=RNA)
        obs_sub_aln = aln.takeSeqs(['rna2','rna3'])
        
        self.assertEqual(obs_sub_aln, sub_aln)
        self.assertEqual(str(obs_sub_aln), str(sub_aln))

        # Selected sequences should be in specified order?
        obs_sub_aln_1 = self.aln.takeSeqs(['rna3','rna2'])
        obs_sub_aln_2 = self.aln.takeSeqs(['rna2','rna3'])
        self.failIfEqual(str(obs_sub_aln_1), str(obs_sub_aln_2))
    
    def test_subset_seqs_DenseAlignment(self):
        model1 = ModelSequence('UCG', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        model2 = ModelSequence('YCG', Name='rna2',\
            Alphabet=RNA.Alphabets.DegenGapped)
        model3 = ModelSequence('CAR', Name='rna3',\
            Alphabet=RNA.Alphabets.DegenGapped)
        sub_da = DenseAlignment([model1, model2, model3],\
            MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)
       
        # takeSeqs by name should have the same effect as
        # getSubAlignment by seq idx?
        obs_sub_da_TS = self.da.takeSeqs(['rna1'])
        obs_sub_da_SA = self.da.getSubAlignment(seqs=[0])

        # These two are now the same. Fixed mapping of key to char array.
        self.assertEqual(obs_sub_da_TS, obs_sub_da_SA)
        self.assertEqual(str(obs_sub_da_TS), str(obs_sub_da_SA))

    def test_aln_equality(self):
        # When does something compare equal?
        self.assertEqual(self.da == self.da, True)
        # one sequence less
        other_da1 = DenseAlignment([self.model1, self.model2],\
            MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)
        self.assertEqual(self.da == other_da1, False)
        # seqs in different order -- doesn't matter
        other_da2 = DenseAlignment([self.model1, self.model3, self.model2],\
            MolType=RNA, Alphabet=RNA.Alphabets.DegenGapped)
        self.assertEqual(self.da == other_da2, True)
        # seqs in different encoding -- doesn't matter, only looks at data
        other_da3 = DenseAlignment([self.model1, self.model2, self.model3])
        # Should this compare False even though the data is exactly the same?
        # The MolType is different...
        self.assertEqual(self.da == other_da3, True) 
        assert alltrue(map(alltrue,self.da.ArraySeqs == other_da3.ArraySeqs))

    def test_seq_equality(self):
        model1 = ModelSequence('UCG', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        model2 = ModelSequence('UCG', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        # Shouldn't the above two sequences be equal?
        self.assertEqual(model1, model2)
        # string comparison is True
        self.assertEqual(str(model1), str(model2))

    def test_seq_ungapping(self):
        rna1 = RnaSequence('U-C-A-G-', Name='rna1')
        model1 = ModelSequence('U-C-A-G-', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        
        self.assertEqual(rna1, 'U-C-A-G-')
        self.assertEqual(rna1.degap(), 'UCAG')

        # check is produces the right string from the beginning
        self.assertEqual(str(model1), 'U-C-A-G-')
        self.assertEqual(model1._data, [0,4,1,4,2,4,3,4])
        # ModelSequence should maybe have the same degap method as normal Seq
        self.assertEqual(str(model1.degap()), 'UCAG')

    def test_the_rest_of_ModelSequence(self):
        """The class ModelSequence has 14 methods, but only 2 unittests.
        You might want to add some tests there..."""
        #note: mostly these are tested in derived classes, for convenience.
        pass

if __name__ == "__main__":
    main()

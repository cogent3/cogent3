#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from numpy import array, float64 as Float64, zeros
from math import sqrt
from random import choice
from cogent.core.alignment import Alignment
from cogent.core.sequence import DnaSequence, RnaSequence
from cogent.core.moltype import DNA, RNA
from cogent.align.weights.util import Weights, number_of_pseudo_seqs,\
    pseudo_seqs_exact, pseudo_seqs_monte_carlo, row_to_vote, distance_matrix,\
    eigenvector_for_largest_eigenvalue, DNA_ORDER,RNA_ORDER,PROTEIN_ORDER,\
    SeqToProfile,AlnToProfile, distance_to_closest

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

class WeightsTests(TestCase):

    def test_weights(self):
        """Weights: should behave like a normal dict and can be normalized
        """
        w = Weights({'seq1':2, 'seq2':3, 'seq3':10})
        self.assertEqual(w['seq1'],2)
        w.normalize()
        exp = {'seq1':0.1333333, 'seq2':0.2, 'seq3':0.6666666}
        self.assertFloatEqual(w.values(), exp.values())
        
class UtilTests(TestCase):
    
    def setUp(self):
        """Set up for Voronoi tests"""
        self.aln1 = Alignment(['ABC','BCC','BAC'])
        
        self.aln2 = Alignment({'seq1':'GYVGS','seq2':'GFDGF','seq3':'GYDGF',\
            'seq4':'GYQGG'},Names=['seq1','seq2','seq3','seq4'])

        self.aln3 = Alignment({'seq1':'AA', 'seq2':'AA', 'seq3':'BB'},\
            Names=['seq1','seq2','seq3'])

        self.aln4 = Alignment({'seq1':'AA', 'seq2':'AA', 'seq3':'BB',\
        'seq4':'BB','seq5':'CC'},Names=['seq1','seq2','seq3','seq4','seq5'])

        self.aln5 = Alignment(['ABBA','ABCA','CBCB'])

    
    def test_number_of_pseudo_seqs(self):
        """number_of_pseudo_seqs: should return # of pseudo seqs"""
        self.assertEqual(number_of_pseudo_seqs(self.aln1),6)
        self.assertEqual(number_of_pseudo_seqs(self.aln2),18)
        self.assertEqual(number_of_pseudo_seqs(self.aln3),4)
        self.assertEqual(number_of_pseudo_seqs(self.aln4),9)
    
    def test_pseudo_seqs_exact(self):
        """pseudo_seqs_exact: should generate expected pseudo sequences"""
        self.assertEqualItems(pseudo_seqs_exact(self.aln1),\
            ['AAC','ABC','ACC','BAC','BBC','BCC']) 
        self.assertEqualItems(pseudo_seqs_exact(self.aln3),\
            ['AA','AB','BA','BB'])
        self.assertEqual(len(pseudo_seqs_exact(self.aln2)), 18)

    def test_pseudo_seqs_monte_carlo(self):
        """pseudo_seqs_monte_carlo: random sample from all possible pseudo seqs
        """
        self.assertEqual(len(list(pseudo_seqs_monte_carlo(self.aln1,n=100))),\
            100)
        for i in pseudo_seqs_monte_carlo(self.aln3,n=100):
            self.assertContains(['AA','AB','BA','BB'], i)

    def test_row_to_vote(self):
        """row_to_vote: should return correct votes for int and float distances
        """
        self.assertEqual(row_to_vote(array([2,3,4,5])),array([1,0,0,0]))
        self.assertEqual(row_to_vote(array([2,3,2,5])),array([.5,0,0.5,0]))
        self.assertEqual(row_to_vote(array([2.3,3.5,2.1,5.8]))\
            ,array([0,0,1,0]))

    def test_distance_matrix(self):
        """distance_matrix should obey Names of alignment"""
        #Names=None
        aln1_exp = array([[0,2,2],[2,0,1],[2,1,0]])
        self.assertEqual(distance_matrix(self.aln1),aln1_exp)
        
        a = Alignment(self.aln1.NamedSeqs)
        a.Names=['seq_1','seq_2','seq_0']
        a_exp = array([[0,1,2],[1,0,2],[2,2,0]])
        self.assertEqual(distance_matrix(a),a_exp)

    def test_eigenvector_for_largest_eigenvalue(self):
        """eigenvector_for_largest_eigenvalue: No idea how to test this"""
        pass

    def test_distance_to_closest(self):
        """distance_to_closest: should return closest distances"""
        self.assertEqual(distance_to_closest(self.aln1),[2,1,1])
        self.assertEqual(distance_to_closest(self.aln2),[2,1,1,2])

    def test_SeqToProfile(self):
        """SequenceToProfile: should work with different parameter settings
        """
        seq = DnaSequence("ATCGRYN-")

        #Only non-degenerate bases in the char order, all other
        #characters are ignored. In a sequence this means that 
        #several positions will contain only zeros in the profile.
        exp = zeros([len(seq),4],Float64)
        for x,y in zip(range(len(seq)),[2,0,1,3]):
            exp[x,y] = 1
        self.assertEqual(SeqToProfile(seq,char_order="TCAG",\
            split_degenerates=False).Data.tolist(),exp.tolist()) 
       
        #Same thing should work as well when the char order is not passed in
        exp = zeros([len(seq),4],Float64)
        for x,y in zip(range(len(seq)),[2,0,1,3]):
            exp[x,y] = 1
        self.assertEqual(SeqToProfile(seq, split_degenerates=False)\
            .Data.tolist(),exp.tolist()) 

       
        #All symbols in the sequence are in the char order, no row
        #should contain only zeros. Degenerate symbols are not split.
        exp = zeros([len(seq),8],Float64)
        for x,y in zip(range(len(seq)),[2,0,1,3,4,5,6,7]):
            exp[x,y] = 1
        self.assertEqual(SeqToProfile(seq,char_order="TCAGRYN-",\
            split_degenerates=False).Data.tolist(), exp.tolist())
        
        #splitting all degenerate symbols, having only non-degenerate symbols
        #in the character order (and -)
        exp = array([[0,0,1,0,0],[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],
            [0,0,.5,.5,0],[.5,.5,0,0,0],[.25,.25,.25,.25,0],[0,0,0,0,1]])
        self.assertEqual(SeqToProfile(seq,char_order="TCAG-",\
            split_degenerates=True).Data.tolist(),exp.tolist())
        
        #splitting degenerates, but having one of the degenerate
        #symbols in the character order. In that case the degenerate symbol
        #is not split. 
        exp = array([[0,0,1,0,0,0],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0],
            [0,0,.5,.5,0,0],[.5,.5,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
        self.assertEqual(SeqToProfile(seq,char_order="TCAGN-",\
            split_degenerates=True).Data.tolist(),exp.tolist())

    def test_AlignmentToProfile_basic(self):
        """AlignmentToProfile: should work under basic conditions
        """
        #sequences in the alignment are unweighted
        #Alphabet is the alphabet of the sequences (RNA)
        #CharOrder is set explicitly
        #Degenerate bases are split up
        #Gaps are ignored
        #In all of the columns at least one character is in the CharOrder
        a = Alignment({'a':RnaSequence('UCAGRYN-'),'b':RnaSequence('ACUGAAAA')})
        exp =\
        array([[.5,0,.5,0],
         [0,1,0,0],
         [.5,0,.5,0],
         [0,0,0,1],
         [0,0,.75,.25],
         [.25,.25,.5,0],
         [.125,.125,.625,.125],
         [0,0,1,0]])
        self.assertEqual(AlnToProfile(a,alphabet=RNA,\
            split_degenerates=True).Data.tolist(),exp.tolist())

    def test_AlignmentToProfile_ignore(self):
        """AlignmentToProfile: should raise an error if too many chars ignored
        """
        #Same conditions as previous function, but in the last column 
        #there are only gaps, so normalization will fail at that position
        a = Alignment({'a':RnaSequence('UCAGRYN-'),'b':RnaSequence('ACUGAAA-')})
        exp =\
        array([[.5,0,.5,0],
         [0,1,0,0],
         [.5,0,.5,0],
         [0,0,0,1],
         [0,0,.75,.25],
         [.25,.25,.5,0],
         [.125,.125,.625,.125],
         [0,0,1,0]])
        self.assertRaises(ValueError,AlnToProfile,a,alphabet=RNA,\
            split_degenerates=True)


    def test_AlignmentToProfile_weighted(self):
        """AlignmentToProfile: should work when sequences are weighted
        """
        #Alignment: sequences are just strings and don't have an alphabet
        #Weights: a normal dictionary (could be a real Weights object as well)
        a = Alignment({'seq1':'TCAG','seq2':'TAR-','seq3':'YAG-'},\
        Names=['seq1','seq2','seq3'])
        w = {'seq1':0.5,'seq2':.25,'seq3':.25}
        
        #Error will be raised when no Alphabet is given, since the seqs
        #in the alignment are just strings
        self.assertRaises(AttributeError,AlnToProfile,a)
        
        #Basic situation in which all letters in the sequences occur in the
        #CharOrder, None have to be ignored. In that case it doesn't matter
        #whether we set split_degenerates to True or False, because if it's 
        #True it's overwritten by the fact that the char is in the CharOrder.
        exp = array([[0.75,0,0,0,0,.25,0],
            [0,0.5,0.5,0,0,0,0],
            [0,0.5,0,0.25,0.25,0,0],
            [0,0,0,0.5,0,0,0.5]])
        #split_degenerates = False
        self.assertEqual(AlnToProfile(a,DNA, char_order="TACGRY-",\
            weights=w, split_degenerates=False).Data.tolist(),exp.tolist())
        #split_degenerates = True
        self.assertEqual(AlnToProfile(a,DNA, char_order="TACGRY-",\
            weights=w, split_degenerates=True).Data.tolist(),exp.tolist())

        #Only non-degenerate symbols in the CharOrder. Degenerates are split.
        #Gaps are ignored
        exp = array([[0.875,0,0.125,0],
            [0,0.5,0.5,0],
            [0,0.625,0,0.375],
            [0,0,0,1]])
        self.assertEqual(AlnToProfile(a,DNA, char_order="TACG",\
            weights=w, split_degenerates=True).Data.tolist(),exp.tolist())
        
        #An Error is raised if all chars in an alignment column are ignored
        #CharOrder=AT, degenerates are not split.
        self.assertRaises(ValueError,AlnToProfile,a,DNA,\
            char_order="AT",weights=w, split_degenerates=True)


if __name__ == "__main__":
    main()

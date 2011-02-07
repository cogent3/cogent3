#!/usr/bin/env python
import numpy
from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs, DNA, RNA
from cogent.evolve.pairwise_distance import get_moltype_index_array, \
    seq_to_indices, _fill_diversity_matrix, \
    _jc69_from_matrix, JC69Pair, _tn93_from_matrix, TN93Pair
from cogent.evolve._pairwise_distance import \
    _fill_diversity_matrix as pyx_fill_diversity_matrix
import math

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class TestPair(TestCase):
    dna_char_indices = get_moltype_index_array(DNA)
    rna_char_indices = get_moltype_index_array(RNA)
    alignment = LoadSeqs(data=[('s1', 'ACGTACGTAC'),
                             ('s2', 'GTGTACGTAC')], moltype=DNA)
    
    ambig_alignment = LoadSeqs(data=[('s1', 'RACGTACGTACN'),
                             ('s2', 'AGTGTACGTACA')], moltype=DNA)
    
    diff_alignment = LoadSeqs(data=[('s1', 'ACGTACGTTT'),
                             ('s2', 'GTGTACGTAC')], moltype=DNA)
    
    def test_char_to_index(self):
        """should correctly recode a DNA & RNA seqs into indices"""
        seq = 'TCAGRNY?-'
        expected = [0, 1, 2, 3, -9, -9, -9, -9, -9]
        indices = seq_to_indices(seq, self.dna_char_indices)
        self.assertEquals(indices, expected)
        seq = 'UCAGRNY?-'
        indices = seq_to_indices(seq, self.rna_char_indices)
        self.assertEquals(indices, expected)
    
    def test_fill_diversity_matrix_all(self):
        """make correct diversity matrix when all chars valid"""
        s1 = seq_to_indices('ACGTACGTAC', self.dna_char_indices)
        s2 = seq_to_indices('GTGTACGTAC', self.dna_char_indices)
        matrix = numpy.zeros((4,4), float)
        # self-self should just be an identity matrix
        _fill_diversity_matrix(matrix, s1, s1)
        self.assertEquals(matrix.sum(), len(s1))
        self.assertEquals(matrix,
            numpy.array([[2,0,0,0],
                         [0,3,0,0],
                         [0,0,3,0],
                         [0,0,0,2]], float))
        
        # small diffs
        matrix.fill(0)
        _fill_diversity_matrix(matrix, s1, s2)
        self.assertEquals(matrix,
            numpy.array([[2,0,0,0],
                         [1,2,0,0],
                         [0,0,2,1],
                         [0,0,0,2]], float))
    
    def test_fill_diversity_matrix_some(self):
        """make correct diversity matrix when not all chars valid"""
        s1 = seq_to_indices('RACGTACGTACN', self.dna_char_indices)
        s2 = seq_to_indices('AGTGTACGTACA', self.dna_char_indices)
        matrix = numpy.zeros((4,4), float)
        # small diffs
        matrix.fill(0)
        _fill_diversity_matrix(matrix, s1, s2)
        self.assertEquals(matrix,
            numpy.array([[2,0,0,0],
                         [1,2,0,0],
                         [0,0,2,1],
                         [0,0,0,2]], float))
    
    def test_python_vs_cython_fill_matrix(self):
        """python & cython fill_diversity_matrix give same answer"""
        s1 = seq_to_indices('RACGTACGTACN', self.dna_char_indices)
        s2 = seq_to_indices('AGTGTACGTACA', self.dna_char_indices)
        matrix1 = numpy.zeros((4,4), float)
        _fill_diversity_matrix(matrix1, s1, s2)
        matrix2 = numpy.zeros((4,4), float)
        pyx_fill_diversity_matrix(matrix2, s1, s2)
        self.assertFloatEqual(matrix1, matrix2)
    
    def test_jc69_from_matrix(self):
        """compute JC69 from diversity matrix"""
        s1 = seq_to_indices('ACGTACGTAC', self.dna_char_indices)
        s2 = seq_to_indices('GTGTACGTAC', self.dna_char_indices)
        matrix = numpy.zeros((4,4), float)
        _fill_diversity_matrix(matrix, s1, s2)
        total, p, dist, var = _jc69_from_matrix(matrix)
        self.assertEquals(total, 10.0)
        self.assertEquals(p, 0.2)
    
    def test_jc69_from_alignment(self):
        """compute JC69 dists from an alignment"""
        calc = JC69Pair(DNA, alignment=self.alignment)
        calc.run()
        self.assertEquals(calc.Lengths['s1', 's2'], 10)
        self.assertEquals(calc.Proportions['s1', 's2'], 0.2)
        # value from OSX MEGA 5
        self.assertFloatEqual(calc.Dists['s1', 's2'], 0.2326161962)
        # value**2 from OSX MEGA 5
        self.assertFloatEqual(calc.Variances['s1', 's2'],
                                0.029752066125078681)
        # value from OSX MEGA 5
        self.assertFloatEqual(calc.StdErr['s1', 's2'], 0.1724878724)
        
        # same answer when using ambiguous alignment
        calc.run(self.ambig_alignment)
        self.assertFloatEqual(calc.Dists['s1', 's2'], 0.2326161962)
        
        # but different answer if subsequent alignment is different
        calc.run(self.diff_alignment)
        self.assertTrue(calc.Dists['s1', 's2'] != 0.2326161962)
    
    def test_tn93_from_matrix(self):
        """compute TN93 distances"""
        calc = TN93Pair(DNA, alignment=self.alignment)
        calc.run()
        self.assertEquals(calc.Lengths['s1', 's2'], 10)
        self.assertEquals(calc.Proportions['s1', 's2'], 0.2)
        # value from OSX MEGA 5
        self.assertFloatEqual(calc.Dists['s1', 's2'], 0.2554128119)
        # value**2 from OSX MEGA 5
        self.assertFloatEqual(calc.Variances['s1', 's2'], 0.04444444445376601)
        # value from OSX MEGA 5
        self.assertFloatEqual(calc.StdErr['s1', 's2'], 0.2108185107)
        
        # same answer when using ambiguous alignment
        calc.run(self.ambig_alignment)
        self.assertFloatEqual(calc.Dists['s1', 's2'], 0.2554128119)
        
        # but different answer if subsequent alignment is different
        calc.run(self.diff_alignment)
        self.assertTrue(calc.Dists['s1', 's2'] != 0.2554128119)
    
    def test_distance_pair(self):
        """get distances dict"""
        calc = TN93Pair(DNA, alignment=self.alignment)
        calc.run()
        dists = calc.getPairwiseDistances()
        dist = 0.2554128119
        expect = {('s1', 's2'): dist, ('s2', 's1'): dist}
        self.assertEquals(dists.keys(), expect.keys())
        self.assertFloatEqual(dists.values(), expect.values())
    

if __name__ == '__main__':
    main()

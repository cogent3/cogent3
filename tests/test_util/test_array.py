#!/usr/bin/env python
"""Provides tests for array.py
"""
#SUPPORT2425
#from __future__ import with_statement
from cogent.util.unit_test import main, TestCase#, numpy_err
from cogent.util.array import gapped_to_ungapped, unmasked_to_masked, \
    ungapped_to_gapped, masked_to_unmasked, pairs_to_array,\
    ln_2, log2, safe_p_log_p, safe_log, row_uncertainty, column_uncertainty,\
    row_degeneracy, column_degeneracy, hamming_distance, norm,\
    euclidean_distance, \
    count_simple, count_alphabet, \
    is_complex, is_significantly_complex, \
    has_neg_off_diags, has_neg_off_diags_naive, \
    sum_neg_off_diags, sum_neg_off_diags_naive, \
    scale_row_sum, scale_row_sum_naive, scale_trace, \
    abs_diff, sq_diff, norm_diff, \
    cartesian_product, with_diag, without_diag, \
    only_nonzero, combine_dimensions, split_dimension, \
    non_diag, perturb_one_off_diag, perturb_off_diag, \
    merge_samples, sort_merged_samples_by_value, classifiers, \
    minimize_error_count, minimize_error_rate, mutate_array

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
from numpy import array, zeros, transpose, sqrt, reshape, arange, \
    ravel, trace, ones

__author__ = "Rob Knight and Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class arrayTests(TestCase):
    """Tests of top-level functions."""
    def setUp(self):
        """set up some standard sequences and masks"""
        self.gap_state = array('-', 'c')
        self.s1 = array('ACT-G', 'c')
        self.s2 = array('--CT', 'c')
        self.s3 = array('AC--', 'c')
        self.s4 = array('AC', 'c')
        self.s5 = array('--', 'c')
        self.m1 = array([0,0,0,1,0])
        self.m2 = array([1,1,0,0])
        self.m3 = array([0,0,1,1])
        self.m4 = array([0,0])
        self.m5 = array([1,1])
    
    def test_unmasked_to_masked(self):
        """unmasked_to_masked should match hand-calculated results"""
        u2m = unmasked_to_masked
        self.assertEqual(u2m(self.m1), array([0,1,2,4]))
        self.assertEqual(u2m(self.m2), array([2,3]))
        self.assertEqual(u2m(self.m3), array([0,1]))
        self.assertEqual(u2m(self.m4), array([0,1]))
        self.assertEqual(u2m(self.m5), array([]))

    def test_ungapped_to_gapped(self):
        """ungapped_to_gapped should match hand-calculated results"""
        u2g = ungapped_to_gapped
        gap_state = self.gap_state
        self.assertEqual(u2g(self.s1, gap_state), array([0,1,2,4]))
        self.assertEqual(u2g(self.s2, gap_state), array([2,3]))
        self.assertEqual(u2g(self.s3, gap_state), array([0,1]))
        self.assertEqual(u2g(self.s4, gap_state), array([0,1]))
        self.assertEqual(u2g(self.s5, gap_state), array([]))

    def test_masked_to_unmasked(self):
        """masked_to_unmasked should match hand-calculated results"""
        m2u = masked_to_unmasked
        self.assertEqual(m2u(self.m1), array([0,1,2,2,3]))
        self.assertEqual(m2u(self.m1, True), array([0,1,2,-1,3]))
        self.assertEqual(m2u(self.m2), array([-1,-1,0,1]))
        self.assertEqual(m2u(self.m2, True), array([-1,-1,0,1]))
        self.assertEqual(m2u(self.m3), array([0,1,1,1]))
        self.assertEqual(m2u(self.m3, True), array([0,1,-1,-1]))
        self.assertEqual(m2u(self.m4), array([0,1]))
        self.assertEqual(m2u(self.m4, True), array([0,1]))
        self.assertEqual(m2u(self.m5), array([-1,-1]))
        self.assertEqual(m2u(self.m5, True), array([-1,-1]))
        
    def test_gapped_to_ungapped(self):
        """gapped_to_ungapped should match hand-calculated results"""
        g2u = gapped_to_ungapped
        gap_state = self.gap_state
        self.assertEqual(g2u(self.s1, gap_state), array([0,1,2,2,3]))
        self.assertEqual(g2u(self.s1, gap_state, True), array([0,1,2,-1,3]))
        self.assertEqual(g2u(self.s2, gap_state), array([-1,-1,0,1]))
        self.assertEqual(g2u(self.s2, gap_state, True), array([-1,-1,0,1]))
        self.assertEqual(g2u(self.s3, gap_state), array([0,1,1,1]))
        self.assertEqual(g2u(self.s3, gap_state, True), array([0,1,-1,-1]))
        self.assertEqual(g2u(self.s4, gap_state), array([0,1]))
        self.assertEqual(g2u(self.s4, gap_state, True), array([0,1]))
        self.assertEqual(g2u(self.s5, gap_state), array([-1,-1]))
        self.assertEqual(g2u(self.s5, gap_state, True), array([-1,-1]))

    def test_pairs_to_array(self):
        """pairs_to_array should match hand-calculated results"""
        p2a = pairs_to_array
        p1 = [0, 1, 0.5]
        p2 = [2, 3, 0.9]
        p3 = [1, 2, 0.6]
        pairs = [p1, p2, p3]
        self.assertEqual(p2a(pairs), \
            array([[0,.5,0,0],[0,0,.6,0],[0,0,0,.9],[0,0,0,0]]))
        #try it without weights -- should assign 1
        new_pairs = [[0,1],[2,3],[1,2]]
        self.assertEqual(p2a(new_pairs), \
            array([[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]]))
        #try it with explicit array size
        self.assertEqual(p2a(pairs, 5), \
            array([[0,.5,0,0,0],[0,0,.6,0,0],[0,0,0,.9,0],[0,0,0,0,0],\
            [0,0,0,0,0]]))
        #try it when we want to map the indices into gapped coords
        #we're effectively doing ABCD -> -A--BC-D-
        transform = array([1,4,5,7])
        result = p2a(pairs, transform=transform)
        self.assertEqual(result.shape, (8,8))
        exp = zeros((8,8), Float)
        exp[1,4] = 0.5
        exp[4,5] = 0.6
        exp[5,7] = 0.9
        self.assertEqual(result, exp)

        result = p2a(pairs, num_items=9, transform=transform)
        self.assertEqual(result.shape, (9,9))
        exp = zeros((9,9), Float)
        exp[1,4] = 0.5
        exp[4,5] = 0.6
        exp[5,7] = 0.9
        self.assertEqual(result, exp)
         
class ArrayMathTests(TestCase):
    
    def test_ln_2(self):
        """ln_2: should be constant"""
        self.assertFloatEqual(ln_2, 0.693147)

    def test_log2(self):
        """log2: should work fine on positive/negative numbers and zero"""
        self.assertEqual(log2(1),0)
        self.assertEqual(log2(2),1)
        self.assertEqual(log2(4),2)
        self.assertEqual(log2(8),3)

        #SUPPORT2425
        #with numpy_err(divide='ignore'):
        ori_err = numpy.geterr()
        numpy.seterr(divide='ignore')
        try:
            try:
                self.assertEqual(log2(0),float('-inf'))
            except (ValueError, OverflowError):      #platform-dependent
                pass
        finally:
            numpy.seterr(**ori_err)

        #SUPPORT2425
        ori_err = numpy.geterr()
        numpy.seterr(divide='raise')
        try:
        #with numpy_err(divide='raise'):
            self.assertRaises(FloatingPointError, log2, 0)
        finally:
            numpy.seterr(**ori_err)

        #nan is the only thing that's not equal to itself
        try:
            self.assertNotEqual(log2(-1),log2(-1)) #now nan
        except ValueError:
            pass

    def test_safe_p_log_p(self):
        """safe_p_log_p: should handle pos/neg/zero/empty arrays as expected
        """
        #normal valid array
        a = array([[4,0,8],[2,16,4]])
        self.assertEqual(safe_p_log_p(a),array([[-8,0,-24],[-2,-64,-8]]))
        #just zeros
        a = array([[0,0],[0,0]])
        self.assertEqual(safe_p_log_p(a),array([[0,0],[0,0]]))
        #negative number -- skip
        self.assertEqual(safe_p_log_p(array([-4])), array([0]))
        #integer input, float output
        self.assertFloatEqual(safe_p_log_p(array([3])),array([-4.75488750]))
        #empty array
        self.assertEqual(safe_p_log_p(array([])),array([]))

    def test_safe_log(self):
        """safe_log: should handle pos/neg/zero/empty arrays as expected
        """
        #normal valid array
        a = array([[4,0,8],[2,16,4]])
        self.assertEqual(safe_log(a),array([[2,0,3],[1,4,2]]))
        #input integers, output floats
        self.assertFloatEqual(safe_log(array([1,2,3])),array([0,1,1.5849625]))
        #just zeros
        a = array([[0,0],[0,0]])
        self.assertEqual(safe_log(a),array([[0,0],[0,0]]))
        #negative number
        try:
            self.assertFloatEqual(safe_log(array([0,3,-4]))[0:2], \
                array([0,1.5849625007]))
        except ValueError:      #platform-dependent
            pass
        try:
            self.assertNotEqual(safe_log(array([0,3,-4]))[2],\
                safe_log(array([0,3,-4]))[2])
        except ValueError:      #platform-dependent
            pass
        #empty array
        self.assertEqual(safe_log(array([])),array([]))
        #double empty array
        self.assertEqual(safe_log(array([[]])),array([[]]))

    def test_row_uncertainty(self):
        """row_uncertainty: should handle pos/neg/zero/empty arrays as expected
        """
        #normal valid array
        b = transpose(array([[.25,.2,.45,.25,1],[.25,.2,.45,0,0],\
            [.25,.3,.05,.75,0],[.25,.3,.05,0,0]]))
        self.assertFloatEqual(row_uncertainty(b),[2,1.97,1.47,0.81,0],1e-3)
        #one-dimensional array
        self.assertRaises(ValueError, row_uncertainty,\
            array([.25,.25,.25,.25]))
        #zeros
        self.assertEqual(row_uncertainty(array([[0,0]])),array([0]))
        #empty 2D array
        self.assertEqual(row_uncertainty(array([[]])),array([0]))
        self.assertEqual(row_uncertainty(array([[],[]])),array([0,0]))
        #negative number -- skip
        self.assertEqual(row_uncertainty(array([[-2]])), array([0]))

    def test_col_uncertainty(self):
        """column_uncertainty: should handle pos/neg/zero/empty arrays
        """
        b = array([[.25,.2,.45,.25,1],[.25,.2,.45,0,0],[.25,.3,.05,.75,0],\
            [.25,.3,.05,0,0]])
        self.assertFloatEqual(column_uncertainty(b),[2,1.97,1.47,0.81,0],1e-3)
        #one-dimensional array
        self.assertRaises(ValueError, column_uncertainty,\
            array([.25,.25,.25,.25]))
        #zeros
        self.assertEqual(column_uncertainty(array([[0,0]])),array([0,0]))
        #empty 2D array
        self.assertEqual(column_uncertainty(array([[]])),array([]))
        self.assertEqual(column_uncertainty(array([[],[]])),array([]))
        #negative number -- skip
        self.assertEqual(column_uncertainty(array([[-2]])), array([0]))

    def test_row_degeneracy(self):
        """row_degeneracy: should work with different cutoff values and arrays
        """
        a = array([[.1, .3, .4, .2],[.5, .3, 0, .2],[.8, 0, .1, .1]])
        self.assertEqual(row_degeneracy(a,cutoff=.75),[3,2,1])
        self.assertEqual(row_degeneracy(a,cutoff=.95),[4,3,3])
        #one-dimensional array
        self.assertRaises(ValueError, row_degeneracy,\
            array([.25,.25,.25,.25]))
        #if cutoff value is not found, results are clipped to the
        #number of columns in the array
        self.assertEqual(row_degeneracy(a,cutoff=2), [4,4,4])
        #same behavior on empty array
        self.assertEqual(row_degeneracy(array([[]])),[])

    def test_column_degeneracy(self):
        """column_degeneracy: should work with different cutoff values
        """
        a = array([[.1,.8,.3],[.3,.2,.3],[.6,0,.4]])
        self.assertEqual(column_degeneracy(a,cutoff=.75),[2,1,3])
        self.assertEqual(column_degeneracy(a,cutoff=.45),[1,1,2])
        #one-dimensional array
        self.assertRaises(ValueError, column_degeneracy,\
            array([.25,.25,.25,.25]))
        #if cutoff value is not found, results are clipped to the
        #number of rows in the array
        self.assertEqual(column_degeneracy(a,cutoff=2), [3,3,3])
        #same behavior on empty array
        self.assertEqual(column_degeneracy(array([[]])),[])

    def test_hamming_distance_same_length(self):
        """hamming_distance: should return # of chars different"""
        hd = hamming_distance(array('ABC','c'),array('ABB','c'))
        self.assertEqual(hd,1)
        self.assertEqual(hamming_distance(array('ABC', 'c'),array('ABC', 'c')),0)
        self.assertEqual(hamming_distance(array('ABC', 'c'),array('DDD', 'c')),3)
       
    def test_hamming_distance_diff_length(self):
        """hamming_distance: truncates at shortest sequence"""
        self.assertEqual(hamming_distance(array('ABC', 'c'),array('ABBDDD', 'c')),1)
        self.assertEqual(hamming_distance(array('ABC', 'c'),array('ABCDDD', 'c')),0)
        self.assertEqual(hamming_distance(array('ABC', 'c'),array('DDDDDD', 'c')),3)

    def test_norm(self):
        """norm: should return vector or matrix norm"""
        self.assertFloatEqual(norm(array([2,3,4,5])),sqrt(54))
        self.assertEqual(norm(array([1,1,1,1])),2)
        self.assertFloatEqual(norm(array([[2,3],[4,5]])),sqrt(54))
        self.assertEqual(norm(array([[1,1],[1,1]])),2)

    def test_euclidean_distance(self):
        """euclidean_distance: should return dist between 2 vectors or matrices
        """
        a = array([3,4])
        b = array([8,5])
        c = array([[2,3],[4,5]])
        d = array([[1,5],[8,2]])
        self.assertFloatEqual(euclidean_distance(a,b),sqrt(26))
        self.assertFloatEqual(euclidean_distance(c,d),sqrt(30))

    def test_euclidean_distance_unexpected(self):
        """euclidean_distance: works always when frames are aligned. UNEXPECTED!
        """
        a = array([3,4])
        b = array([8,5])
        c = array([[2,3],[4,5]])
        d = array([[1,5],[8,2]])
        e = array([[4,5],[4,5],[4,5]])
        f = array([1,1,1,1,1])
        self.assertFloatEqual(euclidean_distance(a,c),sqrt(4))
        self.assertFloatEqual(euclidean_distance(c,a),sqrt(4))
        self.assertFloatEqual(euclidean_distance(a,e),sqrt(6))

        #IT DOES RAISE AN ERROR WHEN THE FRAMES ARE NOT ALIGNED
        self.assertRaises(ValueError,euclidean_distance,c,e)
        self.assertRaises(ValueError,euclidean_distance,c,f)

    def test_count_simple(self):
        """count_simple should return correct counts"""
        self.assertEqual(count_simple(array([]), 3), array([0,0,0]))
        self.assertEqual(count_simple(array([1,2,2,1,0]), 3), array([1,2,2]))
        self.assertEqual(count_simple(array([1,1,1,1,1]), 3), array([0,5,0]))
        self.assertEqual(count_simple(array([1,1,1,1,1]), 2), array([0,5]))
        #raises index error if alphabet length is 0
        self.assertRaises(IndexError, count_simple, array([1]), 0)

    def test_count_alphabet(self):
        """count_alphabet should return correct counts"""
        self.assertEqual(count_alphabet(array([]), 3), array([0,0,0]))
        self.assertEqual(count_alphabet(array([1,2,2,1,0]), 3), array([1,2,2]))
        self.assertEqual(count_alphabet(array([1,1,1,1,1]), 3), array([0,5,0]))
        self.assertEqual(count_alphabet(array([1,1,1,1,1]), 2), array([0,5]))
        #raises index error if alphabet length is 0
        self.assertRaises(IndexError, count_alphabet, array([1]), 0)

    def test_is_complex(self):
        """is_complex should return True on matrix with complex values"""
        self.assertEqual(is_complex(array([[1,2],[3,4]])), False)
        self.assertEqual(is_complex(array([[1,2],[3,4.0]])), False)
        self.assertEqual(is_complex(array([[1,2+1j],[3,4]])), True)
        self.assertEqual(is_complex(array([[1,2.0j],[3,4.0]])), True)

    def test_is_significantly_complex(self):
        """is_significantly_complex should return True on complex matrix"""
        isc = is_significantly_complex
        self.assertEqual(isc(array([[1,2],[3,4]])), False)
        self.assertEqual(isc(array([[1,2],[3,4.0]])), False)
        self.assertEqual(isc(array([[1,2+1j],[3,4]])), True)
        self.assertEqual(isc(array([[1,2.0j],[3,4.0]])), True)
        self.assertEqual(isc(array([[1,1e-10j],[3,4.0]])), False)
        self.assertEqual(isc(array([[1,1e-10j],[3,4.0]]), 1e-12), True)

    def test_has_neg_off_diags_naive(self):
        """has_neg_off_diags_naive should return True if any off-diags negative"""
        hnod = has_neg_off_diags_naive
        self.assertEqual(hnod(array([[1,2],[3,4]])), False)
        self.assertEqual(hnod(array([[-1,2],[3,-4]])), False)
        self.assertEqual(hnod(array([[-1,-2],[3,-4]])), True)
        self.assertEqual(hnod(array([[1,-2],[3,4]])), True)

    def test_has_neg_off_diags(self):
        """has_neg_off_diags should be same as has_neg_off_diags_naive"""
        hnod = has_neg_off_diags
        self.assertEqual(hnod(array([[1,2],[3,4]])), False)
        self.assertEqual(hnod(array([[-1,2],[3,-4]])), False)
        self.assertEqual(hnod(array([[-1,-2],[3,-4]])), True)
        self.assertEqual(hnod(array([[1,-2],[3,4]])), True)

    def test_sum_neg_off_diags_naive(self):
        """sum_neg_off_diags_naive should return the sum of negative off-diags"""
        snod = sum_neg_off_diags_naive
        self.assertEqual(snod(array([[1,2],[3,4]])), 0)
        self.assertEqual(snod(array([[-1,2],[3,-4]])), 0)
        self.assertEqual(snod(array([[-1,-2],[3,-4]])), -2)
        self.assertEqual(snod(array([[1,-2],[3,4]])), -2)
        self.assertEqual(snod(array([[1,-2],[-3,4]])), -5)

    def test_sum_neg_off_diags(self):
        """sum_neg_off_diags should return same as sum_neg_off_diags_naive"""
        snod = sum_neg_off_diags
        self.assertEqual(snod(array([[1,2],[3,4]])), 0)
        self.assertEqual(snod(array([[-1,2],[3,-4]])), 0)
        self.assertEqual(snod(array([[-1,-2],[3,-4]])), -2)
        self.assertEqual(snod(array([[1,-2],[3,4]])), -2)
        self.assertEqual(snod(array([[1,-2],[-3,4]])), -5)

    def test_scale_row_sum(self):
        """scale_row_sum should give same result as scale_row_sum_naive"""
        m = array([[1.0,2,3,4],[2,4,4,0],[1,1,1,1],[0,0,0,100]])
        scale_row_sum(m)
        self.assertFloatEqual(m, [[0.1,0.2,0.3,0.4],[0.2,0.4,0.4,0],\
                [0.25,0.25,0.25,0.25],[0,0,0,1.0]])
        scale_row_sum(m,4)
        self.assertFloatEqual(m, [[0.4,0.8,1.2,1.6],[0.8,1.6,1.6,0],\
                [1,1,1,1],[0,0,0,4.0]])
        #if any of the rows sums to zero, an exception will be raised.

        #SUPPORT2425
        ori_err = numpy.geterr()
        numpy.seterr(divide='raise')
        try:
        #with numpy_err(divide='raise'):
            self.assertRaises((ZeroDivisionError, FloatingPointError), \
                scale_row_sum, array([[1,0],[0,0]]))            
        finally:
            numpy.seterr(**ori_err)


    def test_scale_row_sum_naive(self):
        """scale_row_sum_naive should scale rows to correct values"""
        m = array([[1.0,2,3,4],[2,4,4,0],[1,1,1,1],[0,0,0,100]])
        scale_row_sum_naive(m)
        self.assertFloatEqual(m, [[0.1,0.2,0.3,0.4],[0.2,0.4,0.4,0],\
                [0.25,0.25,0.25,0.25],[0,0,0,1.0]])
        scale_row_sum_naive(m,4)
        self.assertFloatEqual(m, [[0.4,0.8,1.2,1.6],[0.8,1.6,1.6,0],\
                [1,1,1,1],[0,0,0,4.0]])
        #if any of the rows sums to zero, an exception will be raised.

        #SUPPORT2425
        ori_err = numpy.geterr()
        numpy.seterr(divide='raise')
        try:
        #with numpy_err(divide='raise'):
            self.assertRaises((ZeroDivisionError, FloatingPointError), \
                scale_row_sum_naive, array([[1,0],[0,0]]))
        finally:
            numpy.seterr(**ori_err)

    def test_scale_trace(self):
        """scale_trace should scale trace to correct values"""
        #should scale to -1 by default
        #WARNING: won't work with integer matrices
        m = array([[-2., 0],[0,-2]])
        scale_trace(m)
        self.assertFloatEqual(m, [[-0.5, 0],[0,-0.5]])
        #should work even with zero rows
        m = array([
                [1.0,2,3,4],
                [2,4,4,0],
                [1,1,0,1],
                [0,0,0,0]
        ])
        m_orig = m.copy()
        scale_trace(m)
        self.assertFloatEqual(m, m_orig / -5)
        #but should fail if trace is zero
        m = array([[0,1,1],[1,0,1],[1,1,0]])

        #SUPPORT2425
        ori_err = numpy.geterr()
        numpy.seterr(divide='raise')
        try:
        #with numpy_err(divide='raise'):
            self.assertRaises((ZeroDivisionError, FloatingPointError), \
                scale_trace, m)
        finally:
            numpy.seterr(**ori_err)

    def test_abs_diff(self):
        """abs_diff should calculate element-wise sum of abs(first-second)"""
        m =  array([[1.0,2,3],[4,5,6], [7,8,9]])
        m2 = array([[1.0,1,4],[2,6,-1],[8,6,-5]])
        #matrix should not be different from itself
        self.assertEqual(abs_diff(m,m), 0.0)
        self.assertEqual(abs_diff(m2,m2), 0.0)
        #difference should be same either direction
        self.assertEqual(abs_diff(m,m2), 29.0)
        self.assertEqual(abs_diff(m2,m), 29.0)
                
    def test_sq_diff(self):
        """sq_diff should calculate element-wise sum square of abs(first-second)"""     
        m =  array([[1.0,2,3],[4,5,6], [7,8,9]])
        m2 = array([[1.0,1,4],[2,6,-1],[8,6,-5]])
        #matrix should not be different from itself
        self.assertEqual(sq_diff(m,m), 0.0)
        self.assertEqual(sq_diff(m2,m2), 0.0)
        #difference should be same either direction
        self.assertEqual(sq_diff(m,m2), 257.0)
        self.assertEqual(sq_diff(m2,m), 257.0)

    def test_norm_diff(self):
        """norm_diff should calculate per-element rms difference"""     
        m =  array([[1.0,2,3],[4,5,6], [7,8,9]])
        m2 = array([[1.0,1,4],[2,6,-1],[8,6,-5]])
        #matrix should not be different from itself
        self.assertEqual(norm_diff(m,m), 0.0)
        self.assertEqual(norm_diff(m2,m2), 0.0)
        #difference should be same either direction
        self.assertEqual(norm_diff(m,m2), sqrt(257.0)/9)
        self.assertEqual(norm_diff(m2,m), sqrt(257.0)/9)

    def test_carteisan_product(self):
        """cartesian_product should return expected results."""
        a = 'abc'
        b = [1,2,3]
        c = [1.0]
        d = [0,1]
        #cartesian_product of list of single list should be same list
        self.assertEqual(cartesian_product([c]), [(1.0,)])
        self.assertEqual(cartesian_product([a]), [('a',),('b',),('c',)])
        #should combine two lists correctly
        self.assertEqual(cartesian_product([a,b]), \
            [('a',1),('a',2),('a',3),('b',1),('b',2),\
             ('b',3),('c',1),('c',2),('c',3)])
        #should combine three lists correctly
        self.assertEqual(cartesian_product([d,d,d]), \
            [(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)])
        self.assertEqual(cartesian_product([c,d,d]), \
            [(1.0,0,0),(1.0,0,1),(1.0,1,0),(1.0,1,1)])
                
    def test_without_diag(self):
        """without_diag should omit diagonal from matrix"""
        a = array([[1,2,3],[4,5,6],[7,8,9]])
        b = without_diag(a)
        self.assertEqual(b, array([[2,3],[4,6],[7,8]]))

    def test_with_diag(self):
        """with_diag should add diagonal to matrix"""
        a = array([[2,3],[4,6],[7,8]])
        b = with_diag(a, array([1,5,9]))
        self.assertEqual(b, array([[1,2,3],[4,5,6],[7,8,9]]))

    def test_only_nonzero(self):
        """only_nonzero should return only items whose first element is nonzero"""
        a = reshape(arange(1,46),(5,3,3))
        a[1,0,0] = 0
        a[3,0,0] = 0
        #expect result to be rows 0, 2 and 3 of a
        result = only_nonzero(a)
        self.assertEqual(result,
            array([[[1,2,3],[4,5,6],[7,8,9]],\
                [[19,20,21],[22,23,24],[25,26,27]],
                [[37,38,39],[40,41,42],[43,44,45]]]))

    def test_combine_dimensions(self):
        """combine_dimensions should aggregate expected dimensions"""
        m = reshape(arange(81), (3,3,3,3))
        a = combine_dimensions(m, 0)
        self.assertEqual(a.shape, (3,3,3,3))
        a = combine_dimensions(m, 1)
        self.assertEqual(a.shape, (3,3,3,3))
        a = combine_dimensions(m, 2)
        self.assertEqual(a.shape, (9,3,3))
        a = combine_dimensions(m, 3)
        self.assertEqual(a.shape, (27,3))
        a = combine_dimensions(m, 4)
        self.assertEqual(a.shape, (81,))
        #should work for negative indices as well, starting at end
        a = combine_dimensions(m, -1)
        self.assertEqual(a.shape, (3,3,3,3))
        a = combine_dimensions(m, -2)
        self.assertEqual(a.shape, (3,3,9))
        a = combine_dimensions(m, -3)
        self.assertEqual(a.shape, (3,27))
        a = combine_dimensions(m, -4)
        self.assertEqual(a.shape, (81,))

    def test_split_dimension(self):
        """split_dimension should unpack specified dimension"""
        m = reshape(arange(12**3), (12,12,12))
        a = split_dimension(m, 0, (4,3))
        self.assertEqual(a.shape, (4,3,12,12))
        a = split_dimension(m, 0, (2,3,2))
        self.assertEqual(a.shape, (2,3,2,12,12))
        a = split_dimension(m, 1, (6,2))
        self.assertEqual(a.shape, (12, 6, 2, 12))
        a = split_dimension(m, 2, (3,4))
        self.assertEqual(a.shape, (12,12,3,4))
        #should work for negative index
        a = split_dimension(m, -1, (3,4))
        self.assertEqual(a.shape, (12,12,3,4))
        a = split_dimension(m, -2, (3,4))
        self.assertEqual(a.shape, (12,3,4,12))
        a = split_dimension(m, -3, (3,4))
        self.assertEqual(a.shape, (3,4,12,12))
        #should fail with IndexError for invalid dimension
        self.assertRaises(IndexError, split_dimension, m, 5, (3,4))
        #should assume even split if not supplied
        m = reshape(arange(16**3), (16,16,16))
        a = split_dimension(m, 0)
        self.assertEqual(a.shape, (4,4,16,16))
        a = split_dimension(m, 1)
        self.assertEqual(a.shape, (16,4,4,16))

    def test_non_diag(self):
        """non_diag should return non-diag elements from flattened matrices"""
        a = reshape(arange(16), (4,4))
        m = non_diag(a)
        self.assertEqual(m, array([[1,2],[5,6],[9,10],[13,14]]))
        a = reshape(arange(27), (3,9))
        m = non_diag(a)
        self.assertEqual(m, array([[1,2,3,5,6,7],[10,11,12,14,15,16],\
            [19,20,21,23,24,25]]))

    def test_perturb_one_off_diag(self):
        """perturb_element should perturb a random off-diagonal element"""
        for i in range(100):
            a = zeros((4,4), Float)
            p = perturb_one_off_diag(a)
            #NOTE: off-diag element and diag element will _both_ change
            self.assertEqual(sum(ravel(p != a)), 2)
            #check that sum is still 0
            self.assertEqual(sum(ravel(p)), 0)
            #check that rrace is negative
            assert trace(p) < 1
        #check that we can pick an element to change
        a = zeros((4,4), Float)
        p = perturb_one_off_diag(a, mean=5, sd=0.1, element_to_change=8)
        #check that row still sums to 0
        self.assertEqual(sum(ravel(p)), 0)
        #set diag in changed row to 0
        p[2][2] = 0
        assert ((4.5 < sum(p)).any() < 5.5).any()
        assert 4.5 < p[2][3] < 5.5
        p[2][3] = 0
        self.assertEqual(sum(ravel(p)), 0)

    def test_perturb_off_diag(self):
        """perturb_off_diag should change all off_diag elements."""
        a = zeros((4,4), Float)
        d = perturb_off_diag(a)
        self.assertFloatEqual(sum(ravel(d)), 0)
        #try it with a valid rate matrix
        a = ones((4,4), Float)
        for i in range(4):
            a[i][i] = -3
        d = perturb_off_diag(a)
        self.assertNotEqual(d, a)
        self.assertFloatEqual(sum(ravel(d)), 0)
        #check that we didn't change it too much
        assert -13 < trace(d) < -11

    def test_merge_samples(self):
        """merge_samples should keep the sample label"""
        self.assertEqual(merge_samples(array([1,2]),array([3,4]),array([5])),
            array([[1,2,3,4,5],[0,0,1,1,2]]))

    def test_sort_merged_samples_by_value(self):
        """sort_merged_samples_by_value should keep label associations"""
        s = merge_samples(array([3,4]), array([5,6]), array([1,2]))
        result = sort_merged_samples_by_value(s)
        self.assertEqual(result, array([[1,2,3,4,5,6],[2,2,0,0,1,1]]))
            
    def test_classifiers(self):
        """classifiers should return all the 1D classifiers of samples"""
        first = array([2,1,5,3,5])
        second = array([2,5,5,4,6,7])
        result = classifiers(first, second)
        self.assertEqual(len(result), 6)
        exp = [(1,False,0,4,1,6),(3,False,1,3,2,5),(4,False,1,2,3,5),\
               (5,False,2,2,3,4),(9,False,4,0,5,2),(10,False,5,0,5,1)]
        self.assertEqual(result, exp)
        #should work in reverse
        result = classifiers(second, first)
        exp = [(1,True,0,4,1,6),(3,True,1,3,2,5),(4,True,1,2,3,5),\
               (5,True,2,2,3,4),(9,True,4,0,5,2),(10,True,5,0,5,1)]

    def test_minimize_error_count(self):
        """minimize_error_count should return correct classifier"""
        first = array([2,1,5,3,5])
        second = array([2,5,5,4,6,7])
        c = classifiers(first, second)
        exp = (4,False,1,2,3,5)
        self.assertEqual(minimize_error_count(c), exp)

    def test_minimize_error_rate(self):
        """minimize_error_rate should return correct classifier"""
        #should be same as error count on example used above
        first = array([2,1,5,3,5])
        second = array([2,5,5,4,6,7])
        c = classifiers(first, second)
        exp = (4,False,1,2,3,5)
        self.assertEqual(minimize_error_rate(c), exp)
        #here's a case where they should differ
        first = array([2,3,11,5])
        second = array([1,4,6,7,8,9,10])
        c = classifiers(first, second)
        self.assertEqual(minimize_error_count(c), (3,False,1,2,2,6))
        self.assertEqual(minimize_error_rate(c), (5,False,2,1,3,5))

    def test_mutate_array(self):
        """mutate_array should return mutated copy"""
        a = arange(5)
        m = mutate_array(a, 1, 2)
        assert a is not m
        self.assertNotEqual(a, m)
        residuals = m - a
        assert min(residuals) > -6
        assert max(residuals) < 6
        

if __name__ == '__main__':
    main()

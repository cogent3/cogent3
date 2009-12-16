#!/usr/bin/env python
"""Tests of the bitvector module.
"""
from cogent.util.unit_test import TestCase, main
from cogent.core.bitvector import is_nonzero_string_char, is_nonzero_char, \
    seq_to_bitstring, is_nonzero_string_int, is_nonzero_int, seq_to_bitlist,\
    num_to_bitstring, bitcount, Bitvector, MutableBitvector, \
    ImmutableBitvector, VectorFromCases, VectorFromMatches, VectorFromRuns, \
    VectorFromSpans, VectorFromPositions, PackedBases, \
    LongBitvector, ShortBitvector
import re

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class bitvectorTests(TestCase):
    """Tests of top-level functions."""

    def test_is_nonzero_string_char(self):
        """is_nonzero_string_char should return '1' for anything but '0', ''"""
        self.assertEqual(is_nonzero_string_char('0'), '0')
        self.assertEqual(is_nonzero_string_char(''), '0')
        for char in "QWERTYUIOPASDFGHJKL:ZXCGHJMK<L>?{|!@#$%^&*()12345678":
            self.assertEqual(is_nonzero_string_char(char), '1')

    def test_is_nonzero_char(self):
        """is_nonzero_char should return '0' for any False item or '0'"""
        zero = ['', 0, '0', [], {}, None, 0L, 0.0, False]
        for z in zero:
            self.assertEqual(is_nonzero_char(z), '0')
        nonzero = ['z', '1', '00', ' ', 1, -1, 1e-30, [''], {'':None}, True]
        for n in nonzero:
            self.assertEqual(is_nonzero_char(n), '1')

    def test_seq_to_bitstring(self):
        """seq_to_bitstring should provide expected results"""
        zero = ['', 0, '0', [], {}, None, 0L, 0.0, False]
        self.assertEqual(seq_to_bitstring(zero), '0'*9)
        nonzero = ['z', '1', '00', ' ', 1, -1, 1e-30, [''], {'':None}, True]
        self.assertEqual(seq_to_bitstring(nonzero), '1'*10)
        self.assertEqual(seq_to_bitstring(''), '')
        self.assertEqual(seq_to_bitstring('305'), '101')
        self.assertEqual(seq_to_bitstring(''), '')

    def test_is_nonzero_string_int(self):
        """is_nonzero_string_int should return 1 for anything but '0', ''"""
        self.assertEqual(is_nonzero_string_int('0'), 0)
        self.assertEqual(is_nonzero_string_int(''), 0)
        for char in "QWERTYUIOPASDFGHJKL:ZXCGHJMK<L>?{|!@#$%^&*()12345678":
            self.assertEqual(is_nonzero_string_int(char), 1)

    def test_is_nonzero_int(self):
        """is_nonzero_int should return 0 for any False item or '0'"""
        zero = ['', 0, '0', [], {}, None, 0L, 0.0, False]
        for z in zero:
            self.assertEqual(is_nonzero_int(z), 0)
        nonzero = ['z', '1', '00', ' ', 1, -1, 1e-30, [''], {'':None}, True]
        for n in nonzero:
            self.assertEqual(is_nonzero_int(n), 1)

    def test_seq_to_bitlist(self):
        """seq_to_bitlist should provide expected results"""
        zero = ['', 0, '0', [], {}, None, 0L, 0.0, False]
        self.assertEqual(seq_to_bitlist(zero), [0]*9)
        nonzero = ['z', '1', '00', ' ', 1, -1, 1e-30, [''], {'':None}, True]
        self.assertEqual(seq_to_bitlist(nonzero), [1]*10)
        self.assertEqual(seq_to_bitlist(''), [])
        self.assertEqual(seq_to_bitlist('305'), [1,0,1])
        self.assertEqual(seq_to_bitlist(''), [])

    def test_number_to_bitstring(self):
        """number_to_bitstring should provide expected results"""
        numbers = [0, 1, 2, 7, 8, 1024, 814715L]
        for n in numbers:
            self.assertEqual(num_to_bitstring(n, 0), '')

        single_results = list('0101001')
        for exp, num in zip(single_results, numbers):
            self.assertEqual(num_to_bitstring(num, 1), exp)

        three_results = ['000','001','010','111','000','000','011']
        for exp, num in zip(three_results, numbers):
            self.assertEqual(num_to_bitstring(num, 3), exp)
        
        #should pad or truncate to the correct length
        self.assertEqual(num_to_bitstring(814715, 20),'11000110111001111011')
        self.assertEqual(num_to_bitstring(814715, 10),'1001111011')
        self.assertEqual(num_to_bitstring(8, 10),'0000001000')

    def test_bitcount(self):
        """bitcount should provide expected results"""
        numbers = [0, 1, 2, 7, 8, 1024, 814715L]

        twenty_results = [0, 1, 1, 3, 1, 1, 13]
        for exp, num in zip(twenty_results, numbers):
            self.assertEqual(bitcount(num, 20), exp)
            self.assertEqual(bitcount(num, 20, 1), exp)
            self.assertEqual(bitcount(num, 20, 0), 20 - exp)

        three_results = [0,1,1,3,0,0,2]
        for exp, num in zip(three_results, numbers):
            self.assertEqual(bitcount(num, 3), exp)
            self.assertEqual(bitcount(num, 3, 1), exp)
            self.assertEqual(bitcount(num, 3, 0), 3 - exp)

        for num in numbers:
            self.assertEqual(bitcount(num, 0), 0)
            self.assertEqual(bitcount(num, 0, 0), 0)
            self.assertEqual(bitcount(num, 0, 1), 0)
        
class BitvectorTests(TestCase):
    """Tests of the (immutable) Bitvector class."""

    def setUp(self):
        """Define a few standard strings and vectors."""
        self.strings = ['', '0', '1', '00', '01', '10', '11']
        self.vectors = map(Bitvector, self.strings)

    def test_init(self):
        """Bitvector init should give expected results."""
        self.assertEqual(Bitvector(), 0)
        self.assertEqual(Bitvector('1001'), 9)
        self.assertEqual(Bitvector(['1','0','0','0']), 8)
        self.assertEqual(Bitvector([]), 0)
        #if passing in non-sequence, must specify length
        self.assertRaises(TypeError, Bitvector, 1024)
        self.assertEqual(Bitvector(1024, 10), 1024)
        bv = Bitvector(10, 3)
        self.assertEqual(bv, 10)
        self.assertEqual(len(bv), 3)
        self.assertEqual(len(Bitvector('1'*1000)), 1000)
        #check that initializing a bv from itself preserves length
        bv2 = Bitvector(bv)
        self.assertEqual(bv2, 10)
        self.assertEqual(len(bv2), 3)

    def test_len(self):
        """Bitvector len should match initialized length"""
        self.assertEqual(len(Bitvector()), 0)
        self.assertEqual(len(Bitvector('010')), 3)
        self.assertEqual(len(Bitvector(1024, 5)), 5)
        self.assertEqual(len(Bitvector(1024, 0)), 0)
        self.assertEqual(len(Bitvector('1'*1000)), 1000)

    def test_str(self):
        """Bitvector str should match expected results"""
        vecs = [Bitvector(i, 0) for i in [0, 1, 2, 7, 8, 1024, 814715L]]
        for v in vecs:
            self.assertEqual(str(v), '')

        vecs = [Bitvector(i, 1) for i in [0, 1, 2, 7, 8, 1024, 814715L,'1'*50]]
        single_results = list('01010011')
        for exp, vec in zip(single_results, vecs):
            self.assertEqual(str(vec), exp)

        vecs = [Bitvector(i, 3) for i in [0, 1, 2, 7, 8, 1024, 814715L,'1'*50]]
        three_results = ['000','001','010','111','000','000','011','111']
        for exp, vec in zip(three_results, vecs):
            self.assertEqual(str(vec), exp)
        
        #should pad or truncate to the correct length
        self.assertEqual(str(Bitvector(814715, 20)),'11000110111001111011')
        self.assertEqual(str(Bitvector(814715, 10)),'1001111011')
        self.assertEqual(str(Bitvector(8, 10)),'0000001000')
        self.assertEqual(str(Bitvector('1'*50)), '1'*50)

    def test_or(self):
        """Bitvector A|B should return 1 for each position that is 1 in A or B"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' or x
                    ['', '0', '1', '0', '0', '1', '1'],         #'0' or x
                    ['', '1', '1', '1', '1', '1', '1'],         #'1' or x
                    ['', '0', '1', '00', '01', '10', '11'],     #'00' or x
                    ['', '0', '1', '01', '01', '11', '11'],     #'01' or x
                    ['', '1', '1', '10', '11', '10', '11'],     #'10' or x
                    ['', '1', '1', '11', '11', '11', '11'],     #'11' or x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first | second),
                                    results[first_pos][second_pos])
        #test chaining
        expected = Bitvector('1110')
        observed = Bitvector('1000') | Bitvector('0100') | Bitvector('0110')
        self.assertEqual(observed, expected)

        #test long
        self.assertEqual(Bitvector('10'*50) | Bitvector('01'*50), \
            Bitvector('11'*50))
        
    def test_and(self):
        """Bitvector A&B should return 0 for each position that is 0 in A and B"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' and x
                    ['', '0', '0', '0', '0', '0', '0'],         #'0' and x
                    ['', '0', '1', '0', '0', '1', '1'],         #'1' and x
                    ['', '0', '0', '00', '00', '00', '00'],     #'00' and x
                    ['', '0', '0', '00', '01', '00', '01'],     #'01' and x
                    ['', '0', '1', '00', '00', '10', '10'],     #'10' and x
                    ['', '0', '1', '00', '01', '10', '11'],     #'11' and x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first & second), 
                                    results[first_pos][second_pos])
        #test chaining
        expected = Bitvector('0110')
        observed = Bitvector('1110') & Bitvector('1111') & Bitvector('0111')
        self.assertEqual(observed, expected)

        #test long
        self.assertEqual(Bitvector('10'*50) & Bitvector('11'*50), \
            Bitvector('10'*50))
        
    def test_xor(self):
        """Bitvector A^B should return 0 for each identical position in A and B"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' xor x
                    ['', '0', '1', '0', '0', '1', '1'],         #'0' xor x
                    ['', '1', '0', '1', '1', '0', '0'],         #'1' xor x
                    ['', '0', '1', '00', '01', '10', '11'],     #'00' xor x
                    ['', '0', '1', '01', '00', '11', '10'],     #'01' xor x
                    ['', '1', '0', '10', '11', '00', '01'],     #'10' xor x
                    ['', '1', '0', '11', '10', '01', '00'],     #'11' xor x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first ^ second),
                                    results[first_pos][second_pos])

        #test chaining
        expected = Bitvector('0110')
        observed = Bitvector('1111') ^ Bitvector('0110') ^ Bitvector('1111')
 
        #test long
        self.assertEqual(Bitvector('11'*50) ^ Bitvector('01'*50), \
            Bitvector('10'*50))
        
    def test_invert(self):
        """Bitvector ~A should return a vector exchanging 1's for 0's"""
        results = map(Bitvector, ['', '1', '0', '11', '10', '01', '00'])
        for data, result in zip(self.vectors, results):
            self.assertEqual(~data, result)
            
            if len(data):
                self.assertNotEqual(data, result)
            else:
                self.assertEqual(data, result)

            #test chaining
            self.assertEqual(~~data, data) #inverting twice should give original
            self.assertEqual(~~~data, ~data)

            #test long
            self.assertEqual(~Bitvector('10'*50), Bitvector('01'*50))
            self.assertEqual(str(~Bitvector('10'*50)), str(Bitvector('01'*50)))
   
    def test_getitem(self):
        """Bitvector getitem should return states at specified position(s)"""

        vec_strings = ['', '0', '1', '10', '10001101', '101'*50]
        vecs = map(Bitvector, vec_strings)
        for vec_string, vec in zip(vec_strings, vecs):
            for char, item in zip(vec_string, vec):
                self.assertEqual(char, str(item))
        #test some 2- and 3-item slices as well
        vec = Bitvector('1001000101001')
        self.assertEqual(vec[3:7], Bitvector('1000'))
        self.assertEqual(vec[:4], Bitvector('1001'))
        self.assertEqual(vec[7:], Bitvector('101001'))
        self.assertEqual(vec[1:11:2], Bitvector('01011'))
        
    def test_bitcount(self):
        """Bitvector bitcount should correctly count 1's or 0's"""
        vec_strings = ['', '0', '1', '10', '10001101', '101'*50]
        vecs = map(Bitvector, vec_strings)
        one_counts = [0, 0, 1, 1, 4, 100]
        zero_counts = [0, 1, 0, 1, 4, 50]
        for v, o, z in zip(vecs, one_counts, zero_counts):
            self.assertEqual(v.bitcount(), o)
            self.assertEqual(v.bitcount(1), o)
            self.assertEqual(v.bitcount(0), z)

    def test_repr(self):
        """Bitvector repr should look like a normal object"""
        v = Bitvector(3, 10)
        v_id = str(hex(id(v)))
        expected = '<cogent.core.bitvector.ShortBitvector object at'
        self.assertTrue(`v`.startswith(expected))

    def test_freeze(self):
        """Bitvector freeze should return same object"""
        v = Bitvector()
        self.assertSameObj(v.freeze(), v)

    def test_thaw(self):
        """Bitvector thaw should return mutable bitvector with same data"""
        b = Bitvector('111')
        c = b.thaw()
        self.assertEqual(c, Bitvector('111'))
        c[1] = 0
        self.assertEqual(c, Bitvector('101'))
        self.assertNotEqual(Bitvector('111'), Bitvector('101'))

    def test_stateChanges(self):
        """Bitvector stateChanges should return indices where state changes"""
        vec_strings = ['', '0', '1', '10', '111', '10001101', '1111100000'*5]
        vecs = map(Bitvector, vec_strings)
        results = [
            [],
            [0,1],
            [0,1],
            [0,1,2],
            [0,3],
            [0,1,4,6,7,8],
            [0,5,10,15,20,25,30,35,40,45,50],
        ]
        for vec, res in zip(vecs, results):
            self.assertEqual(vec.stateChanges(), res)
    
    def test_divideSequence(self):
        """Bitvector divideSequence should cut sequence at state changes"""
        vec_strings = ['', '0', '1', '10', '111', '10001101', '1111100000'*5,
            '0000011111']
        vecs = map(Bitvector, vec_strings)
        seq = 'abc'*30
        results_none = [
            ([], 0),
            (['a'], 0),
            (['a'], 1),
            (['a','b'], 1),
            (['abc'], 1),
            (['a','bca','bc','a','b'], 1),
            (['abcab','cabca','bcabc','abcab','cabca','bcabc','abcab','cabca',
            'bcabc', 'abcab'], 1),
            (['abcab','cabca'], 0),
        ]
        for vec, res in zip(vecs, results_none):
            self.assertEqual(vec.divideSequence(seq), res)

        short_seq = [1,2,3]
        results_short = [
            ([], 0),
            ([[1]],0),
            ([[1]],1),
            ([[1],[2]],1),
            ([[1,2,3]],1),
            ([[1],[2,3]],1),
            ([[1,2,3]],1),
            ([[1,2,3]],0),
        ]
        for vec, res in zip(vecs, results_short):
            self.assertEqual(vec.divideSequence(short_seq), res)
         
        results_zero = [
            ([], 0),
            (['a'], 0),
            ([], 1),
            (['b'], 1),
            ([], 1),
            (['bca','a'], 1),
            (['cabca','abcab','bcabc','cabca', 'abcab'], 1),
            (['abcab'], 0),
        ]
        for vec, res in zip(vecs, results_zero):
            self.assertEqual(vec.divideSequence(seq, 0), res)

        results_one = [
            ([], 0),
            ([], 0),
            (['a'], 1),
            (['a'], 1),
            (['abc'], 1),
            (['a','bc','b'], 1),
            (['abcab','bcabc','cabca','abcab','bcabc'], 1),
            (['cabca'], 0),
        ]
        for vec, res in zip(vecs, results_none):
            self.assertEqual(vec.divideSequence(seq), res)

class MutableBitvectorTests(TestCase):
    """Tests of the MutableBitvector class."""

    def setUp(self):
        """Define a few standard strings and vectors."""
        self.strings = ['', '0', '1', '00', '01', '10', '11']
        self.vectors = map(MutableBitvector, self.strings)

    def test_init(self):
        """MutableBitvector init should give expected results."""
        self.assertEqual(MutableBitvector(), 0)
        self.assertEqual(MutableBitvector('1001'), 9)
        self.assertEqual(MutableBitvector(['1','0','0','0']), 8)
        self.assertEqual(MutableBitvector([]), 0)
        #if passing in non-sequence, must specify length
        self.assertRaises(TypeError, MutableBitvector, 1024)
        self.assertEqual(MutableBitvector(1024, 10), 1024)
        bv = MutableBitvector(10, 5)
        self.assertEqual(bv, 10)
        self.assertEqual(len(bv), 5)
        self.assertEqual(len(MutableBitvector('1'*1000)), 1000)
        #check that initializing a bv from itself preserves length
        bv2 = MutableBitvector(bv)
        self.assertEqual(bv2, 10)
        self.assertEqual(len(bv2), 5)

    def test_len(self):
        """MutableBitvector len should match initialized length"""
        self.assertEqual(len(MutableBitvector()), 0)
        self.assertEqual(len(MutableBitvector('010')), 3)
        self.assertEqual(len(MutableBitvector(1024, 5)), 5)
        self.assertEqual(len(MutableBitvector(1024, 0)), 0)
        self.assertEqual(len(MutableBitvector('1'*1000)), 1000)

    def test_str(self):
        """MutableBitvector str should match expected results"""
        vecs = [MutableBitvector(i, 0) for i in [0, 1, 2, 7, 8, 1024, 814715L]]
        for v in vecs:
            self.assertEqual(str(v), '')

        vecs = [MutableBitvector(i, 1) for i in [0, 1, 2, 7, 8, 1024, 814715L,'1'*50]]
        single_results = list('01010011')
        for exp, vec in zip(single_results, vecs):
            self.assertEqual(str(vec), exp)

        vecs = [MutableBitvector(i, 3) for i in [0, 1, 2, 7, 8, 1024, 814715L,'1'*50]]
        three_results = ['000','001','010','111','000','000','011','111']
        for exp, vec in zip(three_results, vecs):
            self.assertEqual(str(vec), exp)
        
        #should pad or truncate to the correct length
        self.assertEqual(str(MutableBitvector(814715, 20)),'11000110111001111011')
        self.assertEqual(str(MutableBitvector(814715, 10)),'1001111011')
        self.assertEqual(str(MutableBitvector(8, 10)),'0000001000')
        self.assertEqual(str(MutableBitvector('1'*50)), '1'*50)

    def test_or(self):
        """MutableBitvector A|B should return 1 at 1 positions in A or B"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' or x
                    ['', '0', '1', '0', '0', '1', '1'],         #'0' or x
                    ['', '1', '1', '1', '1', '1', '1'],         #'1' or x
                    ['', '0', '1', '00', '01', '10', '11'],     #'00' or x
                    ['', '0', '1', '01', '01', '11', '11'],     #'01' or x
                    ['', '1', '1', '10', '11', '10', '11'],     #'10' or x
                    ['', '1', '1', '11', '11', '11', '11'],     #'11' or x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first | second),
                                    results[first_pos][second_pos])
        #test chaining
        expected = MutableBitvector('1110')
        observed = MutableBitvector('1000') | \
                   MutableBitvector('0100') | \
                   MutableBitvector('0110')
        self.assertEqual(observed, expected)

        #test long
        self.assertEqual(MutableBitvector('10'*50) | MutableBitvector('01'*50), \
            MutableBitvector('11'*50))
        
    def test_and(self):
        """MutableBitvector A&B should return 0 for each position 0 in both"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' and x
                    ['', '0', '0', '0', '0', '0', '0'],         #'0' and x
                    ['', '0', '1', '0', '0', '1', '1'],         #'1' and x
                    ['', '0', '0', '00', '00', '00', '00'],     #'00' and x
                    ['', '0', '0', '00', '01', '00', '01'],     #'01' and x
                    ['', '0', '1', '00', '00', '10', '10'],     #'10' and x
                    ['', '0', '1', '00', '01', '10', '11'],     #'11' and x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first & second),
                                    results[first_pos][second_pos])
        #test chaining
        expected = MutableBitvector('0110')
        observed = MutableBitvector('1110') & \
                   MutableBitvector('1111') & \
                   MutableBitvector('0111')
        self.assertEqual(observed, expected)
        
        #test long
        self.assertEqual(MutableBitvector('10'*50) & MutableBitvector('11'*50), \
            MutableBitvector('10'*50))
        
    def test_xor(self):
        """MutableBitvector A^B should return 0 at identical positions in A, B"""

        results = [
                    ['', '', '', '', '', '', ''],               #'' xor x
                    ['', '0', '1', '0', '0', '1', '1'],         #'0' xor x
                    ['', '1', '0', '1', '1', '0', '0'],         #'1' xor x
                    ['', '0', '1', '00', '01', '10', '11'],     #'00' xor x
                    ['', '0', '1', '01', '00', '11', '10'],     #'01' xor x
                    ['', '1', '0', '10', '11', '00', '01'],     #'10' xor x
                    ['', '1', '0', '11', '10', '01', '00'],     #'11' xor x
                ]
        vectors = self.vectors

        for first_pos, first in enumerate(vectors):
            for second_pos, second in enumerate(vectors):
                self.assertEqual(   str(first ^ second),
                                    results[first_pos][second_pos])

        #test chaining
        expected = MutableBitvector('0110')
        observed = MutableBitvector('1111') ^ \
                   MutableBitvector('0110') ^ \
                   MutableBitvector('1111')
        self.assertEqual(observed, expected)
 
        #test long
        self.assertEqual(MutableBitvector('11'*50) ^ MutableBitvector('01'*50), \
            MutableBitvector('10'*50))
        
    def test_invert(self):
        """MutableBitvector ~A should return a vector exchanging 1's for 0's"""
        results = map(MutableBitvector, ['', '1', '0', '11', '10', '01', '00'])
        for data, result in zip(self.vectors, results):
            self.assertEqual(~data, result)
            
            if len(data):
                self.assertNotEqual(data, result)
            else:
                self.assertEqual(data, result)

            #test chaining
            self.assertEqual(~~data, data) #inverting twice should give original
            self.assertEqual(~~~data, ~data)

            #test long
            self.assertEqual(~MutableBitvector('10'*50), MutableBitvector('01'*50))
            self.assertEqual(str(~MutableBitvector('10'*50)), str(MutableBitvector('01'*50)))
   
    def test_getitem(self):
        """MutableBitvector getitem should return states at specified position(s)"""

        vec_strings = ['', '0', '1', '10', '10001101', '101'*50]
        vecs = map(MutableBitvector, vec_strings)
        for vec_string, vec in zip(vec_strings, vecs):
            for char, item in zip(vec_string, vec):
                self.assertEqual(char, str(item))
        #test some 2- and 3-item slices as well
        vec = MutableBitvector('1001000101001')
        self.assertEqual(vec[3:7], MutableBitvector('1000'))
        self.assertEqual(vec[:4], MutableBitvector('1001'))
        self.assertEqual(vec[7:], MutableBitvector('101001'))
        self.assertEqual(vec[1:11:2], MutableBitvector('01011'))
       

    def test_setitem(self):
        """MutableBitvector setitem should change positions correctly"""
        self.assertRaises(IndexError, MutableBitvector().__setitem__, 1, 1)
        vec_strings = ['0', '1', '10', '10001101', '101'*50]
        vecs = map(MutableBitvector, vec_strings)
        results_1 = ['1', '1', '10', '10001101', '101'*50]
        for vec, res  in zip(vecs, results_1):
            vec[0] = 1
            self.assertEqual(vec, Bitvector(res))

        results_0 = ['0', '0', '10', '10001100', '101'*49+'100']
        for vec, res  in zip(vecs, results_0):
            vec[-1] = 0
            self.assertEqual(vec, Bitvector(res))
            
        vec = MutableBitvector('1001000101001')
        vec[3:7] = '1111'
        self.assertEqual(vec,  Bitvector('1001111101001'))
        vec[:4] = '0101'
        self.assertEqual(vec,  Bitvector('0101111101001'))
        vec[7:] = '000111'
        self.assertEqual(vec,  Bitvector('0101111000111'))
        vec[1:11:2] = '11011'
        self.assertEqual(vec,  Bitvector('0101101101111'))

    def test_bitcount(self):
        """MutableBitvector bitcount should correctly count 1's or 0's"""
        vec_strings = ['', '0', '1', '10', '10001101', '101'*50]
        vecs = map(MutableBitvector, vec_strings)
        one_counts = [0, 0, 1, 1, 4, 100]
        zero_counts = [0, 1, 0, 1, 4, 50]
        for v, o, z in zip(vecs, one_counts, zero_counts):
            self.assertEqual(v.bitcount(), o)
            self.assertEqual(v.bitcount(1), o)
            self.assertEqual(v.bitcount(0), z)

    def test_repr(self):
        """MutableBitvector repr should look like a normal object"""
        v = MutableBitvector(3, 10)
        v_id = str(hex(id(v)))
        expected = '<cogent.core.bitvector.MutableBitvector object at'
        self.assertTrue(`v`.startswith(expected))

    def test_thaw(self):
        """MutableBitvector thaw should return same object"""
        v = MutableBitvector()
        self.assertSameObj(v.thaw(), v)

    def test_freeze(self):
        """MutableBitvector freeze should return immutable bv with same data"""
        b = MutableBitvector('111')
        b[1] = 0
        c = b.freeze()
        self.assertSameObj(c, b._handler)
        self.assertEqual(c, Bitvector('101'))
        try:
            c[1] = 1
        except TypeError:
            pass
        else:
            raise AssertionError, \
                "MutableBitvector.freeze() returned mutable object."

    def test_stateChanges(self):
        """MutableBitvector stateChanges should return indices where state changes"""
        vec_strings = ['', '0', '1', '10', '111', '10001101', '1111100000'*5]
        vecs = map(MutableBitvector, vec_strings)
        results = [
            [],
            [0,1],
            [0,1],
            [0,1,2],
            [0,3],
            [0,1,4,6,7,8],
            [0,5,10,15,20,25,30,35,40,45,50],
        ]
        for vec, res in zip(vecs, results):
            self.assertEqual(vec.stateChanges(), res)
    
    def test_divideSequence(self):
        """MutableBitvector divideSequence should cut sequence at state changes"""
        vec_strings = ['', '0', '1', '10', '111', '10001101', '1111100000'*5,
            '0000011111']
        vecs = map(MutableBitvector, vec_strings)
        seq = 'abc'*30
        results_none = [
            ([], 0),
            (['a'], 0),
            (['a'], 1),
            (['a','b'], 1),
            (['abc'], 1),
            (['a','bca','bc','a','b'], 1),
            (['abcab','cabca','bcabc','abcab','cabca','bcabc','abcab','cabca',
            'bcabc', 'abcab'], 1),
            (['abcab','cabca'], 0),
        ]
        for vec, res in zip(vecs, results_none):
            self.assertEqual(vec.divideSequence(seq), res)

        short_seq = [1,2,3]
        results_short = [
            ([], 0),
            ([[1]],0),
            ([[1]],1),
            ([[1],[2]],1),
            ([[1,2,3]],1),
            ([[1],[2,3]],1),
            ([[1,2,3]],1),
            ([[1,2,3]],0),
        ]
        for vec, res in zip(vecs, results_short):
            self.assertEqual(vec.divideSequence(short_seq), res)
         
        results_zero = [
            ([], 0),
            (['a'], 0),
            ([], 1),
            (['b'], 1),
            ([], 1),
            (['bca','a'], 1),
            (['cabca','abcab','bcabc','cabca', 'abcab'], 1),
            (['abcab'], 0),
        ]
        for vec, res in zip(vecs, results_zero):
            self.assertEqual(vec.divideSequence(seq, 0), res)

        results_one = [
            ([], 0),
            ([], 0),
            (['a'], 1),
            (['a'], 1),
            (['abc'], 1),
            (['a','bc','b'], 1),
            (['abcab','bcabc','cabca','abcab','bcabc'], 1),
            (['cabca'], 0),
        ]
        for vec, res in zip(vecs, results_none):
            self.assertEqual(vec.divideSequence(seq), res)

class BitvectorClassTests(TestCase):
    """Bit operations should work correctly for different bitvector types."""
    def setUp(self):
        """Define a few standard vectors"""
        self.s = ShortBitvector('00000000000000101111000000001011000')
        self.l = LongBitvector('11111111111111111111111111111111111')
        self.tiny = ShortBitvector('1')
        self.huge = LongBitvector('1'*1000)

    def test_and(self):
        """Bitwise and should work between short and long bitvectors"""
        s, l, tiny, huge = self.s, self.l, self.tiny, self.huge
        self.assertTrue(isinstance(l, long))
        self.assertEqual(l & s, s)
        self.assertEqual(s & l, s)
        self.assertEqual(huge & tiny, tiny)
        self.assertEqual(tiny & huge, tiny)

    def test_or(self):
        """Bitwise or should work between short and long bitvectors"""
        s, l, tiny, huge = self.s, self.l, self.tiny, self.huge
        self.assertEqual(l | s, l)
        self.assertEqual(s | l, l)
        self.assertEqual(huge | tiny, tiny)
        self.assertEqual(tiny | huge, tiny)

class VectorFromCasesTests(TestCase):
    """Tests of the VectorFromCases factory fuction."""

    def test_init(self):
        """VectorFromCases should return vector with 1 where string is ucase"""
        valid_strings = ['', 'a', 'X', 'aBc', 'Acb', 'abC', 'AAA', 'aaa', '@']
        results =       ['', '0', '1', '010', '100', '001', '111', '000', '0']

        for data, result in zip(valid_strings, results):
            self.assertEqual(VectorFromCases(data), LongBitvector(result))
            self.assertTrue(isinstance(VectorFromCases(data), LongBitvector))
            self.assertEqual(VectorFromCases(data, ShortBitvector), \
                ShortBitvector(result))
            self.assertTrue(isinstance(VectorFromCases(data, ShortBitvector), \
                                       ShortBitvector))

        v = VectorFromCases('aBC')
        self.assertTrue(isinstance(v, ImmutableBitvector))

        w = VectorFromCases('aBC', MutableBitvector)
        self.assertTrue(isinstance(w, MutableBitvector))
        self.assertEqual(w, Bitvector('011'))
        w[0] = 1
        self.assertEqual(w, Bitvector('111'))

class VectorFromMatchesTests(TestCase):
    """Tests of the VectorFromMatches factory function.
    
    Need to check all combinations of the following cases:
        1. Pattern is:
            -empty
            -single-character string
            -single-character regex
            -multi-character string
            -multi-character regex
            -regex that includes alternation

        2. String is:
            -empty
            -match at every position
            -match at some positions
            -match at no position

        3. Match is:
            -overlapping
            -non-overlapping
    """
    def testBothEmpty(self):
        """VectorFromMatches empty string/pattern should return empty vector"""
        self.assertEqual(str(VectorFromMatches('', '')), '')
        
    def testEmptyPattern(self):
        """VectorFromMatches empty pattern should return zeroes for len(string)"""
        sequences = ['', 'a', 'aa', 'aaaaaaaaaa']
        for s in sequences:
            vec = VectorFromMatches(s, '')
            self.assertEqual(str(vec), '0' * len(s))

    def testSingleBasePattern(self):
        """VectorFromMatches should match every matching char in string"""
        sequences = ['', 'a', 'b', 'aaa', 'bbb', 'aba', 'bab']
        a_matches = ['', '1', '0', '111', '000', '101', '010']
        b_matches = ['', '0', '1', '000', '111', '010', '101']
        
        for i, s in enumerate(sequences):
            vec = VectorFromMatches(s, 'a')
            self.assertEqual(str(vec), a_matches[i])
            vec = VectorFromMatches(s, 'b')
            self.assertEqual(str(vec), b_matches[i])

    def testMultiBasePattern(self):
        """VectorFromMatches should match multi-char string matches"""
        pattern = 'aba'
        sequences = ['','a', 'aba', 'abab', 'ababa', 'ababab', 'abaaba', 'aaba']
        overlap =   ['','0', '111', '1110', '11111', '111110', '111111', '0111']
        discrete =  ['','0', '111', '1110', '11100', '111000', '111111', '0111']
        for i, s in enumerate(sequences):
            vec = VectorFromMatches(s, pattern)
            self.assertEqual(str(vec), overlap[i])
            vec = VectorFromMatches(s, pattern, 1)
            self.assertEqual(str(vec), overlap[i])
            vec = VectorFromMatches(s, pattern, 0)
            self.assertEqual(str(vec), discrete[i])

    def testSingleBaseRegex(self):
        """VectorFromMatches should match every matching character in regex"""
        sequences  = ['', 'a', 'b', 'aaa', 'bbb', 'aba', 'bab', 'axb', 'xxx']
        a_matches  = ['', '1', '0', '111', '000', '101', '010', '100', '000']
        b_matches  = ['', '0', '1', '000', '111', '010', '101', '001', '000']
        ab_matches = ['', '1', '1', '111', '111', '111', '111', '101', '000']
        
        a = re.compile('a')
        b = re.compile('b')
        ab = re.compile('a|b')
        
        for i, s in enumerate(sequences):
            #test that a works as regex or list
            vec = VectorFromMatches(s, a)
            self.assertEqual(str(vec), a_matches[i])
            vec = VectorFromMatches(s, ['a'])
            self.assertEqual(str(vec), a_matches[i])
            
            #test that b works as regex or list
            vec = VectorFromMatches(s, b)
            self.assertEqual(str(vec), b_matches[i])
            vec = VectorFromMatches(s, ['b'])
            self.assertEqual(str(vec), b_matches[i])

            #test that [a or b] works as regex or list
            vec = VectorFromMatches(s, ab)
            self.assertEqual(str(vec), ab_matches[i])
            vec = VectorFromMatches(s, ['a', 'b'])
            self.assertEqual(str(vec), ab_matches[i])

    def testMultiBaseRegex(self):
        """VectorFromMatches should match every matching combination of chars"""
        sequences = ['aaabbb', 'aaaxbbb', 'ababab', 'abaabaabab']
        patterns =  ['aaa', 'bbb', 'aba', 'aaa|bbb', 'aaa|aab']

        overlap = { 'aaa'   :   [
                                    '111000',    #aaabbb
                                    '1110000',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'bbb'   :   [
                                    '000111',    #aaabbb
                                    '0000111',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'aba'   :   [
                                    '000000',    #aaabbb
                                    '0000000',   #aaaxbbb
                                    '111110',    #ababab
                                    '1111111110',#abaabaabab
                                ],
                    'aaa|bbb'   :   [
                                    '111111',    #aaabbb
                                    '1110111',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'aaa|aab'   :   [
                                    '111100',    #aaabbb
                                    '1110000',   #aaaxbbb
                                    '000000',    #ababab
                                    '0011111100',#abaabaabab
                                ]
                    }

        no_overlap = { 'aaa'   :   [
                                    '111000',    #aaabbb
                                    '1110000',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'bbb'   :   [
                                    '000111',    #aaabbb
                                    '0000111',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'aba'   :   [
                                    '000000',    #aaabbb
                                    '0000000',   #aaaxbbb
                                    '111000',    #ababab
                                    '1111111110',#abaabaabab
                                ],
                    'aaa|bbb'   :   [
                                    '111111',    #aaabbb
                                    '1110111',   #aaaxbbb
                                    '000000',    #ababab
                                    '0000000000',#abaabaabab
                                ],
                    'aaa|aab'   :   [
                                    '111000',    #aaabbb
                                    '1110000',   #aaaxbbb
                                    '000000',    #ababab
                                    '0011111100',#abaabaabab
                                ],
                    }
        for i, s in enumerate(sequences):
            for pat in patterns:
                regex = re.compile(pat)
                vec = VectorFromMatches(s, regex, 1)  #overlapping
                self.assertEqual(str(vec), overlap[pat][i])
                vec = VectorFromMatches(s, regex, 0)  #non-overlapping
                self.assertEqual(str(vec), no_overlap[pat][i])
    
    def test_type(self):
        """VectorFromMatches should return correct type of vector"""
        v = VectorFromMatches('a', 'a')
        self.assertEqual(v, Bitvector('1'))
        self.assertTrue(isinstance(v, ImmutableBitvector))
        v = VectorFromMatches('a', 'a', constructor=MutableBitvector)
        self.assertEqual(v, Bitvector('1'))
        v[0] = 0
        self.assertEqual(v, Bitvector('0'))
        self.assertTrue(isinstance(v, MutableBitvector))

class VectorFromRunsTests(TestCase):
    """Tests of the VectorFromRuns factory fuction."""

    def test_init(self):
        """VectorFromRuns should return vector with 1 where string is ucase"""
        empty = []
        empty_run = [[5,0]]
        one_run = [[5,3]]
        two_runs = [[5,1], [11,4]]
        overlap = [[5,3], [6,4]]

        tests = [empty, empty_run, one_run, two_runs, overlap]
        exp_15 = ['0'*15, '0'*15, '0'*5+'1'*3+7*'0', '0'*5+'1'+'0'*5+'1'*4, 
            '0'*5+'1'*5+'0'*5]
        for test, exp in zip(tests, exp_15):
            self.assertEqual(VectorFromRuns(test, 15), Bitvector(exp))
            self.assertEqual(VectorFromRuns(test, 0), Bitvector(''))
        #should fail if vec too short
        self.assertRaises(IndexError, VectorFromRuns, two_runs, 5)
        
class VectorFromSpansTests(TestCase):
    """Tests of the VectorFromSpans factory fuction."""

    def test_init(self):
        """VectorFromSpans hould return vector with 1 in correct spans"""
        empty = []
        empty_run = [[5,5]]
        one_run = [[5,8]]
        two_runs = [[5,6], [11,15]]
        overlap = [[5,8], [6,10]]

        tests = [empty, empty_run, one_run, two_runs, overlap]
        exp_15 = ['0'*15, '0'*15, '0'*5+'1'*3+7*'0', '0'*5+'1'+'0'*5+'1'*4,
            '0'*5+'1'*5+'0'*5]
        for test, exp in zip(tests, exp_15):
            self.assertEqual(VectorFromSpans(test, 15), Bitvector(exp))
            self.assertEqual(VectorFromSpans(test, 0), Bitvector(''))
        #should fail if vec too short
        self.assertRaises(IndexError, VectorFromSpans, two_runs, 5)

class VectorFromPositionsTests(TestCase):
    """Tests of the VectorFromPositions factory function."""
    def test_init(self):
        """VectorFromPositions should return correct vector"""
        empty = []
        first = [0]
        fifth = [4]
        several = [1,4,6,9]

        tests = [empty, first, fifth, several]
        exp_0 = ['']*4
        for test, exp in zip(tests, exp_0):
            self.assertEqual(VectorFromPositions(test, 0), Bitvector(exp))

        exp_5 = ['00000', '10000', '00001']
        for test, exp in zip(tests, exp_5):
            self.assertEqual(VectorFromPositions(test, 5), Bitvector(exp))

        exp_10 = ['0'*10, '1'+'0'*9, '00001'+'0'*5, '0100101001']
        for test, exp in zip(tests, exp_10):
            self.assertEqual(VectorFromPositions(test, 10), Bitvector(exp))
        self.assertRaises(IndexError, VectorFromPositions, [10], 5)

class PackedBasesTests(TestCase):
    """Tests of the PackedBases class."""
    
    def test_init(self):
        """PackedBases init should allow sequence recovery"""
        p = PackedBases()
        self.assertEqual(str(p), '')
        self.assertEqual(len(p), 0)
        p = PackedBases('uCaGGCAU')
        self.assertEqual(len(p), 16)
        self.assertEqual(str(p), 'UCAGGCAU')
        p = PackedBases('aaaAaaAaa')
        self.assertEqual(str(p), 'AAAAAAAAA')
        self.assertEqual(len(p), 18)
        for base in 'ucagUCAGtT':
            p = PackedBases(base)
            self.assertEqual(len(p), 2)
            self.assertEqual(str(p), base.upper().replace('T','U'))

    def test_str(self):
        """PackedBases str should respect RNA/DNA"""
        p = PackedBases('UCAGu')
        self.assertEqual(str(p), 'UCAGU')
        p.Rna = False
        self.assertEqual(str(p), 'TCAGT')
        p.Rna = True
        self.assertEqual(str(p), 'UCAGU')
        q = PackedBases('ucagu', Rna=False)
        self.assertEqual(str(q), 'TCAGT')
        p = PackedBases('T'*100, Rna=False)
        self.assertEqual(str(p), 'T'*100)
        p.Rna = True
        self.assertEqual(str(p), 'U'*100)

    def test_bit_ops(self):
        """PackedBases xor and bitcount should give correct # differences"""
        p = PackedBases('UCAGU')
        q = PackedBases('CAGGA')
        self.assertEqual((p^q).bitcount(), 5)
        r = PackedBases('CAAG')
        self.assertEqual((r^q).bitcount(), 1)
        
    
#main loop of program
if __name__ == "__main__":
    main()

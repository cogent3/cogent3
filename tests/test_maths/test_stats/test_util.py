#!/usr/bin/env python
"""Tests Numbers and Freqs objects, and their Unsafe versions.
"""

from math import sqrt
from cogent.util.unit_test import TestCase, main
from cogent.maths.stats.util import SummaryStatistics, SummaryStatisticsError,\
        Numbers, UnsafeNumbers, Freqs, UnsafeFreqs, NumberFreqs, \
        UnsafeNumberFreqs
from cogent.util.misc import ConstraintError
from operator import add, sub, mul

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


class SummaryStatisticsTests(TestCase):
    """Tests of summary stats functions."""
    def test_init(self):
        """SummaryStatistics should initialize correctly."""
        #check empty init -- can access private vars, but can't get
        #properties.
        s = SummaryStatistics()
        self.assertEqual(s._count, None)
        self.assertRaises(SummaryStatisticsError, getattr, s, 'Count')
        #check init with one positional parameter
        s = SummaryStatistics(1)
        self.assertEqual(s.Count, 1)
        #check init with all positional parameters.
        #note that inconsistent data can sneak in (c.f. sd vs var)
        s = SummaryStatistics(1,2,3,4,5,6)
        self.assertEqual(s.Count, 1)
        self.assertEqual(s.Sum, 2)
        self.assertEqual(s.Mean, 3)
        self.assertEqual(s.StandardDeviation, 4)
        self.assertEqual(s.Variance, 5)
        self.assertEqual(s.SumSquares, 6)

    def test_str(self):
        """SummaryStatistics str should print known fields."""
        s = SummaryStatistics()
        self.assertEqual(str(s), '')
        #note that additional fields will fill in if they can be calculated.
        s = SummaryStatistics(Mean=3, StandardDeviation=2)
        #now expect to print as table
        self.assertEqual(str(s), '==========================\n        Statistic    Value\n--------------------------\n             Mean        3\nStandardDeviation        2\n         Variance        4\n--------------------------')

    def test_Count(self):
        """SummaryStatistics Count should work if Count or Sum and Mean ok"""
        s = SummaryStatistics(Count=3)
        self.assertEqual(s.Count, 3)
        s = SummaryStatistics(Sum=10, Mean=5)
        self.assertEqual(s.Count, 2)
        #if inconsistent, believes Count
        s = SummaryStatistics(Count=3, Sum=2, Mean=5)
        self.assertEqual(s.Count, 3)
        #doesn't work with just sum or mean
        s = SummaryStatistics(Mean=3)
        self.assertRaises(SummaryStatisticsError, getattr, s, 'Count')

    def test_Sum(self):
        """SummaryStatistics Sum should work if Sum or Count and Mean ok"""
        s = SummaryStatistics(Sum=3)
        self.assertEqual(s.Sum, 3)
        s = SummaryStatistics(Count=3, Mean=5)
        self.assertEqual(s.Sum, 15)
        
    def test_Mean(self):
        """SummaryStatistics Mean should work if Mean or Count and Sum ok"""
        s = SummaryStatistics(Mean=3)
        self.assertEqual(s.Mean, 3)
        s = SummaryStatistics(Count=3, Sum=15)
        self.assertEqual(s.Mean, 5)

    def test_StandardDeviation(self):
        """SummaryStatistics StandardDeviation should work if it or variance ok"""
        s = SummaryStatistics(StandardDeviation=3)
        self.assertEqual(s.StandardDeviation, 3)
        self.assertEqual(s.Variance, 9)
        s = SummaryStatistics(Variance=9)
        self.assertEqual(s.StandardDeviation, 3)
        
    def test_Variance(self):
        """SummaryStatistics Variance should work if it or std dev ok"""
        s = SummaryStatistics(StandardDeviation=3)
        self.assertEqual(s.StandardDeviation, 3)
        self.assertEqual(s.Variance, 9)
        s = SummaryStatistics(Variance=9)
        self.assertEqual(s.StandardDeviation, 3)
 
    def test_SumSquares(self):
        """SummaryStatistics SumSquares should work if set"""
        s = SummaryStatistics(SumSquares=3)
        self.assertEqual(s.SumSquares, 3)
        s = SummaryStatistics(Sum=3)
        self.assertRaises(SummaryStatisticsError, getattr, s, 'SumSquares')

    def test_cmp(self):
        """SummaryStatistics should sort by count, then sum, then variance"""
        a = SummaryStatistics(Count=3)
        b = SummaryStatistics(Count=4)
        c = SummaryStatistics(Count=3, Sum=5)
        d = SummaryStatistics(Count=3, Sum=10)
        e = SummaryStatistics(Sum=10)

        assert a < b
        assert b > a
        assert a == a
        assert a < b
        assert c < d
        assert e < a

        all = [c,a,d,b,e]
        all.sort()
        self.assertEqual(all, [e,a,c,d,b])

class NumbersTestsI(object):
    """Abstract class with tests for Numbers objects.

    Inherited by safe and unsafe versions to test polymorphism.
    """
    ClassToTest = None

    def test_init_empty(self):
        """Numbers should initialize OK with empty list"""
        self.assertEqual(self.ClassToTest([]), [])
        
    def test_init_single(self):
        """Numbers should initialize OK with single number"""
        self.assertEqual(self.ClassToTest([5.0]), [5.0])

    def test_init_list(self):
        """Numbers should initialize OK with list of numbers"""
        self.assertEqual(self.ClassToTest([1, 5.0, 3.2]), [1, 5.0, 3.2])

    def test_init_bad_type(self):
        """Numbers should fail with TypeError if input not iterable"""
        self.assertRaises(TypeError, self.ClassToTest, 34)

    def test_add_nonempty(self):
        """Numbers should allow addition of two nonempty Numbers"""
        #test that addition works in the right direction
        self.assertFloatEqual(self.integers + self.floats,
            Numbers([1, 2, 3, 4, 5, 1.5, 2.7]))
        #test that neither of the things that were added was changed
        self.assertFloatEqual(self.integers, [1,2,3,4,5])
        self.assertFloatEqual(self.floats, [1.5, 2.7])

    def test_add_empty(self):
        """Numbers should be unchanged on addition of empty list"""
        #test that addition of an empty list works
        self.assertFloatEqual(self.integers + self.empty, self.integers)
        self.assertFloatEqual(self.empty + self.floats, self.floats)
    
    def test_add_repeated(self):
        """Numbers should support repeated addition, a+b+c"""
        self.assertFloatEqual(self.floats + self.floats + self.floats,
            [1.5, 2.7]*3)

    def test_iadd(self):
        """Numbers should support in-place addition"""
        self.floats += [4]
        self.assertFloatEqual(self.floats, [1.5, 2.7, 4.0])

    def test_setitem(self):
        """Numbers should support assignment to positive index"""
        self.floats[0] = 1
        self.assertFloatEqual(self.floats, [1.0, 2.7])

    def test_setitem_negative_index(self):
        """Numbers should support assignment to negative index"""
        self.floats[-1] = 2
        self.assertFloatEqual(self.floats, [1.5, 2.0])

    def test_setslice(self):
        """Numbers should support slice assignment"""
        self.floats[0:1] = [1, 2, 3]
        self.assertFloatEqual(self.floats, [1, 2, 3, 2.7])
    
    def test_append_good(self):
        """Numbers should support append of a number"""
        self.floats.append(1)
        self.assertFloatEqual(self.floats,
            [1.5, 2.7, 1.0])

    def test_extend(self):
        """Numbers should support extend with a sequence"""
        self.floats.extend([5,5,5])
        self.assertFloatEqual(self.floats,
            [1.5, 2.7, 5.0, 5.0, 5.0])

    def test_items(self):
        """Numbers should support items() method"""
        self.assertFloatEqual(self.floats.items()[0], (1.5, 1))
        self.assertFloatEqual(self.floats.items()[1], (2.7, 1))

    def test_isValid(self):
        """Numbers isValid should return True if all items numbers"""
        for i in [self.empty, self.integers, self.floats, self.mixed]:
            assert i.isValid()

    def test_toFixedWidth(self):
        """Numbers should be able to convert items to fixed-width string"""
        self.assertEqual(self.floats.toFixedWidth(), " +1.50e+00 +2.70e+00")

    def test_toFixedWidth_empty(self):
        """Numbers should return empty string when converting no items"""
        self.assertEqual(self.empty.toFixedWidth(), '')

    def test_toFixedWidth_mixed(self):
        """Numbers should convert all kinds of floats to fixed precision"""
        self.assertEqual(self.mixed.toFixedWidth(), ''.join([
                                            ' +0.00e+00',
                                            ' +1.00e+00',
                                            ' -1.00e+00',
                                            ' +1.23e+00',
                                            ' -1.24e+00',
                                            '+1.23e+302',
                                            '+1.23e-298',
                                            '-1.23e+302',
                                            '-1.23e-298',
                                            ]))
 
    def test_toFixedWidth_specified_precision(self):
        """Numbers should convert all kinds of floats to specified precision"""
        self.assertEqual(self.mixed.toFixedWidth(7), ''.join([
                                            ' +0e+00',
                                            ' +1e+00',
                                            ' -1e+00',
                                            ' +1e+00',
                                            ' -1e+00',
                                            '+1e+302',
                                            '+1e-298',
                                            '-1e+302',
                                            '-1e-298',
                                            ]))

        self.assertEqual(self.mixed.toFixedWidth(8), ''.join([
                                            '  +0e+00',
                                            '  +1e+00',
                                            '  -1e+00',
                                            '  +1e+00',
                                            '  -1e+00',
                                            ' +1e+302',
                                            ' +1e-298',
                                            ' -1e+302',
                                            ' -1e-298',
                                            ]))
 
        self.assertEqual(self.mixed.toFixedWidth(12), ''.join([
                                            ' +0.0000e+00',
                                            ' +1.0000e+00',
                                            ' -1.0000e+00',
                                            ' +1.2346e+00',
                                            ' -1.2368e+00',
                                            '+1.2340e+302',
                                            '+1.2340e-298',
                                            '-1.2340e+302',
                                            '-1.2340e-298',
                                            ]))

    def test_normalize(self):
        """Numbers normalize should return items summing to 1 by default"""
        first = self.ints
        second = self.fracs
        first.normalize()
        second.normalize()
        self.assertFloatEqual(first, second)
        self.assertFloatEqual(first.Sum, 1)
        self.assertFloatEqual(second.Sum, 1)
        empty = self.empty
        empty.normalize()
        self.assertEqual(empty, [])
        zero = self.zero
        zero.normalize()
        self.assertEqual(zero, [0,0,0,0,0])
       
    def test_normalize_parameter(self):
        """Numbers normalize(x) should divide items by x"""
        first = self.ClassToTest([0, 1, 2, 3, 4])
        first.normalize(max(first))
        self.assertFloatEqual(first, [0, 1.0/4, 2.0/4, 3.0/4, 4.0/4])
        second = self.ClassToTest([0, 1, 2])
        second.normalize(0.5)
        self.assertFloatEqual(second, [0, 2, 4])

    def test_accumulate(self):
        """Numbers accumulate should do cumulative sum in place"""
        nl = self.ClassToTest([0, 1, 2, 3, 4])
        nl.accumulate()
        self.assertEqual(nl, [0, 1, 3, 6, 10])
        nl = self.ClassToTest()
        nl.accumulate()
        self.assertEqual(nl, [])

    def test_firstIndexLessThan(self):
        """Numbers firstIndexLessThan should return first index less than val"""
        nl = self.ints
        f = nl.firstIndexLessThan
        self.assertEqual(f(-50), None)
        self.assertEqual(f(100), 0)
        self.assertEqual(f(3), 0)
        self.assertEqual(f(1), None)
        self.assertEqual(f(1, inclusive=True), 0)
        self.assertEqual(f(-50, stop_at_ends=True), 4)

    def test_firstIndexGreaterThan(self):
        """Numbers firstIndexGreaterThan should return first index less than val"""
        nl = self.ints
        f = nl.firstIndexGreaterThan
        self.assertEqual(f(-50), 0)
        self.assertEqual(f(100), None)
        self.assertEqual(f(3), 3)
        self.assertEqual(f(1), 1)
        self.assertEqual(f(1, inclusive=True), 0)
        self.assertEqual(f(2), 2)
        self.assertEqual(f(2, inclusive=True), 1)
        self.assertEqual(f(100, stop_at_ends=True), 4)

        #compatibility tests with old choose()
        """Numbers choose should return correct index"""
        nl = self.ClassToTest([1, 2, 3, 4, 5])
        nl.normalize()
        nl.accumulate()
        known_values = [
                        (-50, 0),
                        (0, 0),
                        (0.001, 0),
                        (1/15.0 - 0.001, 0),
                        (1/15.0 + 0.001, 1),
                        (3/15.0 + 0.001, 2),
                        (1, 4),
                        (10, 4),
                     ]
        for test, result in known_values:
            self.assertFloatEqual(nl.firstIndexGreaterThan(test, inclusive=True, stop_at_ends=True), result)

    def test_lastIndexGreaterThan(self):
        """Numbers lastIndexGreaterThan should return last index > val"""
        nl = self.ints
        f = nl.lastIndexGreaterThan
        self.assertEqual(f(-50), 4)
        self.assertEqual(f(100), None)
        self.assertEqual(f(3), 4)
        self.assertEqual(f(1), 4)
        self.assertEqual(f(1, inclusive=True), 4)
        self.assertEqual(f(100, stop_at_ends=True), 0)

    def test_lastIndexLessThan(self):
        """Numbers lastIndexLessThan should return last index < val"""
        nl = self.ints
        f = nl.lastIndexLessThan
        self.assertEqual(f(-50), None)
        self.assertEqual(f(100), 4)
        self.assertEqual(f(3), 1)
        self.assertEqual(f(1), None)
        self.assertEqual(f(1, inclusive=True), 0)
        self.assertEqual(f(-50, stop_at_ends=True), 0)

    def test_Sum(self):
        """Numbers Sum should be the same as sum()"""
        self.assertEqual(self.ints.Sum, 15)
        self.assertEqual(self.empty.Sum, 0)

    def test_Count(self):
        """Numbers Count should be the same as len()"""
        self.assertEqual(self.ints.Count, 5)
        self.assertEqual(self.empty.Count, 0)

    def test_SumSquares(self):
        """Numbers SumSquares should be sum of squares"""
        self.assertEqual(self.ints.SumSquares, (1*1+2*2+3*3+4*4+5*5))
        self.assertEqual(self.empty.SumSquares, 0)

    def test_Variance(self):
        """Numbers Variance should be variance of individual numbers"""
        self.assertEqual(self.empty.Variance, None)
        self.assertEqual(self.zero.Variance, 0)
        self.assertFloatEqual(self.ints.Variance, 2.5)

    def test_StandardDeviation(self):
        """Numbers StandardDeviation should be sd of individual numbers"""
        self.assertEqual(self.empty.StandardDeviation, None)
        self.assertEqual(self.zero.StandardDeviation, 0)
        self.assertFloatEqual(self.ints.StandardDeviation, sqrt(2.5))

    def test_Mean(self):
        """Numbers Mean should be mean of individual numbers"""
        self.assertEqual(self.empty.Mean, None)
        self.assertEqual(self.zero.Mean, 0)
        self.assertEqual(self.ints.Mean, 3)

    def test_summarize(self):
        """Numbers summarize should return SummaryStatistics object"""
        self.assertEqual(self.ints.summarize(), SummaryStatistics(Mean=3,\
            Variance=2.5, Count=5))

    def test_choice(self):
        """Numbers choice should return random element from self"""
        nums = [self.ints.choice() for i in range(10)]
        self.assertEqual(len(nums), 10)
        for n in nums:
            assert n in self.ints
        v = Numbers(nums).Variance
        self.assertNotEqual(v, 0)

    def test_randomSequence(self):
        """Numbers randomSequence should return random sequence from self"""
        nums = self.ints.randomSequence(10)
        nums = [self.ints.choice() for i in range(10)]
        self.assertEqual(len(nums), 10)
        for n in nums:
            assert n in self.ints
        v = Numbers(nums).Variance
        self.assertNotEqual(v, 0)
        
    def test_subset(self):
        """Numbers subset should delete (or keep) selected items"""
        odd = [5,1,3]
        nums = self.ints
        nums.extend([1,1,1])
        new_nums = nums.copy()
        new_nums.subset(odd)
        self.assertEqual(new_nums, [1,3,5,1,1,1])
        new_nums = nums.copy()
        new_nums.subset(odd, keep=False)
        self.assertEqual(new_nums, [2,4])

    def test_copy(self):
        """Numbers copy should leave class intact (unlike slice)"""
        c = self.ints.copy()
        self.assertEqual(c, self.ints)
        self.assertEqual(c.__class__, self.ints.__class__)

    def test_round(self):
        """Numbers round should round numbers in-place"""
        self.floats.round()
        self.assertEqual(self.floats, [2.0,3.0])
        for i, f in enumerate(self.floats):
            self.floats[i] = self.floats[i] + 0.101
        self.assertNotEqual(self.floats, [2.0,3.0])
        self.assertNotEqual(self.floats, [2.1,3.1])
        self.floats.round(1)
        self.assertEqual(self.floats, [2.1,3.1])

    def test_Uncertainty(self):
        """Numbers Uncertainty should act via Freqs"""
        self.assertEqual(self.floats.Uncertainty, \
                Freqs(self.floats).Uncertainty)
        self.assertNotEqual(self.floats.Uncertainty, None)

    def test_Mode(self):
        """Numbers Mode should return most common element"""
        self.assertEqual(self.empty.Mode, None)
        self.assertEqual(self.zero.Mode, 0)
        self.ints.extend([1,2,2,3,3,3])
        self.assertEqual(self.ints.Mode, 3)

class NumbersTests(TestCase, NumbersTestsI):
    """Tests of the (safe) Numbers class."""
    ClassToTest = Numbers

    def setUp(self):
        """define some standard lists"""
        self.empty = self.ClassToTest([])
        self.integers = self.ClassToTest([1,2,3,4,5])
        self.floats = self.ClassToTest([1.5, 2.7])
        self.mixed = self.ClassToTest([
                            0,
                            1,
                            -1,
                            1.234567890,
                            -1.2367890,
                            123.4e300,
                            123.4e-300,
                            -123.4e300,
                            -123.4e-300,
                        ])
        self.zero = self.ClassToTest([0,0,0,0,0])
        self.ints = self.ClassToTest([1,2,3,4,5])
        self.fracs = self.ClassToTest([0.1,0.2,0.3,0.4,0.5])

    def test_init_string(self):
        """Numbers should initialize by treating string as list of digits"""
        self.assertEqual(self.ClassToTest('102'), [1.0, 0.0, 2.0])

    def test_init_bad_string(self):
        """Numbers should raise ValueError if float() can't convert string"""
        self.assertRaises(ValueError, self.ClassToTest, '102a')

    def test_append_bad(self):
        """Numbers should reject append of a non-number"""
        self.assertRaises(ValueError, self.floats.append, "abc")

class UnsafeNumbersTests(TestCase, NumbersTestsI):
    """Tests of the UnsafeNumbers class."""
    ClassToTest = UnsafeNumbers

    def setUp(self):
        """define some standard lists"""
        self.empty = self.ClassToTest([])
        self.integers = self.ClassToTest([1,2,3,4,5])
        self.floats = self.ClassToTest([1.5, 2.7])
        self.mixed = self.ClassToTest([
                            0,
                            1,
                            -1,
                            1.234567890,
                            -1.2367890,
                            123.4e300,
                            123.4e-300,
                            -123.4e300,
                            -123.4e-300,
                        ])
        self.zero = self.ClassToTest([0,0,0,0,0])
        self.ints = self.ClassToTest([1,2,3,4,5])
        self.fracs = self.ClassToTest([0.1,0.2,0.3,0.4,0.5])

    def test_init_string(self):
        """UnsafeNumbers should treat string as list of chars"""
        self.assertEqual(self.ClassToTest('102'), ['1','0','2'])

    def test_init_bad_string(self):
        """UnsafeNumbers should silently incorporate unfloatable string"""
        self.assertEqual(self.ClassToTest('102a'), ['1','0','2','a'])

    def test_append_bad(self):
        """UnsafeNumbers should allow append of a non-number"""
        self.empty.append('abc')
        self.assertEquals(self.empty, ['abc'])

    def test_isValid_bad(self):
        """UnsafeNumbers should return False if invalid"""
        assert self.mixed.isValid()
        self.mixed.append('abc')
        assert not self.mixed.isValid()

class StaticFreqsTestsI(object):
    """Tests of the interface shared by Freqs and UnsafeFreqs (static keys).
    
    All of these tests assume keys on the alphabet 'abcde'.

    These tests were added to ensure that array-based objects that implement
    a fixed set of keys maintain the appropriate portion of the freqs interface.
    """
    ClassToTest = None

    def setUp(self):
        """Standard cases to test."""
        self.Alphabetic = self.ClassToTest({'a':3,'b':2,'c':1,'d':1,'e':1})
        self.Empty = self.ClassToTest({})
        self.Constant = self.ClassToTest({'a':5})

    #The following test various ways of constructing the objects

    def test_fromTuples(self):
        """Freqs fromTuples should add from key, count pairs w/ repeated keys"""
        ct = self.ClassToTest
        f = ct()
        self.assertEqual(f.fromTuples([('a',4),('b',3),('a',2)]), \
            ct({'a':6,'b':3}))
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative.
        f.fromTuples([('b',-1),('c',4.5)])
        self.assertEqual(f, ct({'a':6,'b':2,'c':4.5}))
        #should work with a different operator
        f.fromTuples([('b',7)], op=mul)
        self.assertEqual(f, ct({'a':6, 'b':14, 'c':4.5}))
        #check that it works with something that depends on the key
        def func(key, first, second):
            if key == 'a':
                return first + second
            else:
                return max(second, first * second)
        f = ct()
        self.assertEqual(f.fromTuples([('a',4),('b',3),('a',2), ('b',4)], \
            func,uses_key=True), ct({'a':6,'b':12}))
        
    def test_newFromTuples(self):
        """Freqs newFromTuples should work as expected."""
        ct = self.ClassToTest
        self.assertEqual(ct.newFromTuples([('a',4),('b',3),('a',2)]), \
            ct({'a':6,'b':3}))
        
    def test_fromDict(self):
        """Freqs fromDict should add from dict of {key:count}"""
        ct = self.ClassToTest
        f = ct()
        self.assertEqual(f.fromDict({'a':6,'b':3}), ct({'a':6,'b':3}))
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative.
        f.fromDict({'b':-1, 'c':4.5})
        self.assertEqual(f, ct({'a':6,'b':2,'c':4.5}))
        #should work with a different operator
        f.fromDict({'b':7}, op=mul)
        self.assertEqual(f, ct({'a':6, 'b':14, 'c':4.5}))

    def test_newFromDict(self):
        """Freqs newFromDict should work as expected."""
        ct = self.ClassToTest
        self.assertEqual(ct.newFromDict({'a':6,'b':3}), ct({'a':6,'b':3}))
     
    def test_fromDicts(self):
        """Freqs fromDicts should add from list of dicts of {key:count}"""
        ct = self.ClassToTest
        f = ct()
        self.assertEqual(f.fromDicts([{'a':6},{'b':3}]), ct({'a':6,'b':3}))
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative. Also tests add of 1-item dict (note: must be in list)
        f.fromDicts([{'b':-1, 'c':4.5}])
        self.assertEqual(f, ct({'a':6,'b':2,'c':4.5}))
        #should work with a different operator
        f.fromDicts([{'b':2},{'b':3}], op=mul)
        self.assertEqual(f, ct({'a':6, 'b':12, 'c':4.5}))

    def test_newFromDicts(self):
        """Freqs newFromDicts should work as expected."""
        ct = self.ClassToTest
        self.assertEqual(ct.newFromDicts([{'a':6},{'b':3}]), ct({'a':6,'b':3}))
 
    def test_fromSeq(self):
        """Freqs fromSeq should add items from sequence, according to weight"""
        ct = self.ClassToTest
        f = self.ClassToTest()
        self.assertEqual(f.fromSeq('aaabbbaaa'), ct({'a':6,'b':3}))
        #should be able to change the operator...
        self.assertEqual(f.fromSeq('aab', sub), ct({'a':4,'b':2}))
        #...or change the weight
        self.assertEqual(f.fromSeq('acc', weight=3.5),\
            ct({'a':7.5,'b':2,'c':7}))

    def test_newFromSeq(self):
        """Freqs newFromSeq should work as expected."""
        ct = self.ClassToTest
        self.assertEqual(ct.newFromSeq('aaabbbaaa'), ct({'a':6,'b':3}))

    def test_fromSeqs(self):
        """Freqs fromSeqs should add items from sequences, according to weight"""
        ct = self.ClassToTest
        f = ct()
        self.assertEqual(f.fromSeqs(['aaa','bbbaaa']), ct({'a':6,'b':3}))
        #should be able to change the operator...
        self.assertEqual(f.fromSeqs(list('aab'), sub), ct({'a':4,'b':2}))
        #...or change the weight. Note that a string counts as a seq of seqs.
        self.assertEqual(f.fromSeqs('acc', weight=3.5), \
            ct({'a':7.5,'b':2,'c':7}))

    def test_newFromSeqs(self):
        """Freqs newFromSeqs should work as expected."""
        ct = self.ClassToTest
        self.assertEqual(ct.newFromSeqs(['aaa','bbbaaa']), ct({'a':6,'b':3}))

    def test_isValid(self):
        """Freqs isValid should return True if valid"""
        d =self.ClassToTest()
        assert d.isValid()
        d.fromSeq('aaaaaaaaaaaaabb')
        assert d.isValid()

    def test_find_conversion_function(self):
        """Freqs _find_conversion_function should return correct value."""
        d = self.ClassToTest()
        f = d._find_conversion_function
        #should always return None if data empty
        for i in [None, 0, False, {}, [], tuple()]:
            self.assertEqual(f(i), None)
        #should return fromDict for non-empty dict
        self.assertEqual(f({3:4}), d.fromDict)
        #should return fromSeq for string or list of scalars or strings
        for i in ['abc', [1,2,3], (1,2,3), ['ab','bb','cb']]:
            self.assertEqual(f(i), d.fromSeq)
        #should return fromSeqs for sequence of sequences
        for i in [[[1,2,3],[3,4,4]], ([1,2,4],[3,4,4]), [(1,2),[3],[], [4]]]:
            self.assertEqual(f(i), d.fromSeqs)
        #should return fromTuples if possibly key-value pairs
        for i in [[('a',3),('b',-1)], [(1,2),(3,4)]]:
            self.assertEqual(f(i), d.fromTuples)
        #should not be fooled by 2-item seqs that can't be key-value pairs
        self.assertEqual(f(['ab','cd']), d.fromSeq)

    #The following test inheritance of dict properties/methods

    def test_setitem_good(self):
        """Freqs should allow non-negative values to be set"""
        ct = self.ClassToTest
        self.Empty['a'] = 0
        self.assertEqual(self.Empty, ct({'a':0}))
        self.Empty['b'] = 5
        self.assertEqual(self.Empty, ct({'a':0, 'b':5}))

    def test_delitem(self):
        """delitem not applicable to freqs w/ constant keys: not tested"""
        pass

    def test_setdefault_good(self):
        """Freqs setdefault should work with positive values if key present"""
        a = self.Alphabetic.setdefault('a', 200)
        self.assertEqual(a, 3)
        self.assertEqual(self.Alphabetic['a'], 3)

    def test_iadd(self):
        """Freqs += should add in place from any known data type"""
        ct = self.ClassToTest
        f = ct({'a':3, 'b':4})
        f += 'aca'
        self.assertEqual(f, ct({'a':5, 'b':4, 'c':1}))
        f += ['b','b']
        self.assertEqual(f, ct({'a':5, 'b':6, 'c':1}))
        f += {'c':10, 'a':-3}
        self.assertEqual(f, ct({'a':2, 'b':6, 'c':11}))
        f += (('a',3),('b',-2))
        self.assertEqual(f, ct({'a':5, 'b':4, 'c':11}))
        f += [['a', 'b', 'b'],['c', 'c', 'c']]
        self.assertEqual(f, ct({'a':6, 'b':6, 'c':14}))
        #note that list of strings will give implementation-dependent result

    def test_add(self):
        """Freqs + should make new object, adding from any known data type"""
        ct = self.ClassToTest
        orig = {'a':3, 'b':4}
        f = ct(orig)
        r = f + 'aca'
        self.assertEqual(r, ct({'a':5, 'b':4, 'c':1}))
        self.assertEqual(f, orig)
        r = f + ['b','b']
        self.assertEqual(r, ct({'a':3, 'b':6}))
        self.assertEqual(f, orig)
        r = f + {'c':10, 'a':-3}
        self.assertEqual(r, ct({'a':0, 'b':4, 'c':10}))
        self.assertEqual(f, orig)
        r = f + (('a',3),('b',-2))
        self.assertEqual(r, ct({'a':6, 'b':2}))
        self.assertEqual(f, orig)
        r = f + [['a', 'b', 'b'],['c', 'c', 'c']]
        self.assertEqual(r, ct({'a':4, 'b':6, 'c':3}))
        self.assertEqual(f, orig)
        #note that list of strings will give implementation-dependent result
        
    def test_isub(self):
        """Freqs -= should subtract in place using any known data type"""
        ct = self.ClassToTest
        f = ct({'a':5, 'b':4})
        f -= 'aba'
        self.assertEqual(f, ct({'a':3, 'b':3}))
        f -= ['b','b']
        self.assertEqual(f, ct({'a':3, 'b':1}))
        f -= {'c':-2, 'a':-3}
        self.assertEqual(f, ct({'a':6, 'b':1, 'c':2}))
        f -= (('a',3),('b',-2))
        self.assertEqual(f, ct({'a':3, 'b':3, 'c':2}))
        f -= [['a', 'b', 'b'],['c', 'c']]
        self.assertEqual(f, ct({'a':2, 'b':1, 'c':0}))
        #note that list of strings will give implementation-dependent result

    def test_sub(self):
        """Freqs - should make new object, subtracting using any known data type"""
        orig = {'a':3, 'b':4}
        ct = self.ClassToTest
        f = self.ClassToTest(orig)
        r = f - 'aba'
        self.assertEqual(r, ct({'a':1, 'b':3}))
        self.assertEqual(f, orig)
        r = f - ['b','b']
        self.assertEqual(r, ct({'a':3, 'b':2}))
        self.assertEqual(f, orig)
        r = f - {'c':-10, 'a':3}
        self.assertEqual(r, ct({'a':0, 'b':4, 'c':10}))
        self.assertEqual(f, orig)
        r = f - (('a',3),('b',-2))
        self.assertEqual(r, ct({'a':0, 'b':6}))
        self.assertEqual(f, orig)
        r = f - [['a', 'b', 'b'],['a','a']]
        self.assertEqual(r, ct({'a':0, 'b':2}))
        self.assertEqual(f, orig)
        #note that list of strings will give implementation-dependent results
        
    def test_copy(self):
        """Freqs copy should preserve class of original"""
        d = {'a':4, 'b':3, 'c':6}
        f = self.ClassToTest(d)
        g = f.copy()
        self.assertEqual(f, g)
        self.assertEqual(f.__class__, g.__class__)

    def test_str(self):
        """Freqs abstract interface doesn't specify string result"""
        pass

    def test_delitem(self):
        """Freqs delitem is implementation-dependent"""
        pass

    #The following test custom methods
        
    def test_rekey(self):
        """Freqs rekey should map the results onto new keys"""
        #note that what happens to unmapped keys is implementation-dependent
        ct = self.ClassToTest
        d = ct({'a':3, 'b':5, 'c':6, 'd':7, 'e':1})
        #should work with simple rekeying
        f = d.rekey({'a':'d', 'b':'e'})
        self.assertEqual(f['d'], d['a'])
        self.assertEqual(f['e'], d['b'])
        #remaining keys might be absent or 0
        for i in 'abc':
            if i in f:
                self.assertEqual(f[i], 0)
        #should work if many old keys map to the same new key
        f = d.rekey({'a':'d', 'b':'e', 'c':'e'})
        self.assertEqual(f['d'], d['a'])
        self.assertEqual(f['e'], d['b'] + d['c'])
        #remaining keys might be absent or 0
        for i in 'abc':
            if i in f:
                self.assertEqual(f[i], 0)
        
        #check with explicit constructor and default
        d = self.ClassToTest({'a':3, 'b':5, 'c':6, 'd':7, 'e':1})
        f = d.rekey({'a':'+', 'b':'-', 'c':'+'}, default='x', constructor=dict)
        self.assertEqual(f, {'+':9, '-':5, 'x':8})
        self.assertEqual(f.__class__, dict)

    def test_purge(self):
        """Freqs purge should have no effect if keys are fixed"""
        ct = self.ClassToTest
        orig = {'a':3, 'b':2, 'c':1, 'd':3, 'e':4}
        d = ct(orig)
        d.purge()
        self.assertEqual(d, ct(orig))
       
    def test_normalize(self):
        """Freqs should allow normalization"""
        ct = self.ClassToTest
        self.Empty.normalize()
        self.assertEqual(self.Empty, ct({}))
        a = self.Alphabetic.copy()
        a.normalize()
        expected = {'a':0.375, 'b':0.25, 'c':0.125, 'd':0.125, 'e':0.125}
        for key, val in expected.items():
            self.assertFloatEqual(a[key], val)

    def test_choice(self):
        """Freqs choice should work as expected"""
        self.Alphabetic.normalize()
        keys = self.Alphabetic.keys()
        vals = Numbers(self.Alphabetic.values())
        vals.accumulate()
        #test first item
        self.assertEqual(self.Alphabetic.choice(-1), keys[0])
        self.assertEqual(self.Alphabetic.choice(-0.0001), keys[0])
        self.assertEqual(self.Alphabetic.choice(-1e300), keys[0])
        self.assertEqual(self.Alphabetic.choice(0), keys[0])
        #test last item
        last_val = vals.pop()
        self.assertEqual(self.Alphabetic.choice(last_val), keys[-1])
        self.assertEqual(self.Alphabetic.choice(1000), keys[-1])
        #test remaining items
        for index, value in enumerate(vals):
            self.assertEqual(self.Alphabetic.choice(value-0.01),keys[index])
            self.assertEqual(self.Alphabetic.choice(value+0.01), keys[index+1])
            
    def test_randomSequence_good(self):
        """Freqs randomSequence should give correct counts"""
        self.Alphabetic.normalize()
        total = self.Alphabetic.Sum
        keys = self.Alphabetic.keys()
        probs = [float(i)/total for i in self.Alphabetic.values()]
        
        rand_seq = self.Alphabetic.randomSequence(10000)
        observed = [rand_seq.count(key) for key in keys]
        expected = [prob*10000 for prob in probs]
        self.assertSimilarFreqs(observed, expected)

    def test_randomSequence_bad(self):
        """Empty Freqs should raise error on randomSequence"""
        self.assertRaises(IndexError, self.Empty.randomSequence, 5)

    def test_randomSequence_one_item(self):
        """Freqs randomSequence should work with one key"""
        self.Constant.normalize()
        rand = self.Constant.randomSequence(1000)
        self.assertEqual(rand.count('a'), 1000)
        self.assertEqual(len(rand), 1000)

    def test_subset_preserve(self):
        """Freqs subset should preserve wanted items"""
        ct = self.ClassToTest
        self.Constant.subset('bc')
        self.assertEqual(self.Constant, self.Empty)
        self.Alphabetic.subset('abx')
        self.assertEqual(self.Alphabetic, ct({'a':3,'b':2}))

    def test_subset_remove(self):
        """Freqs subset should delete unwanted items if asked"""
        ct = self.ClassToTest
        self.Alphabetic.subset('abx', keep=False)
        self.assertEqual(self.Alphabetic, ct({'c':1,'d':1,'e':1}))
        self.Constant.subset('bx', keep=False)
        self.assertEqual(self.Constant, ct({'a':5}))

    def test_scale(self):
        """Freqs scale should multiply all values with the given scale"""
        ct = self.ClassToTest
        f = ct({'a':0.25,'b':0.25})
        f.scale(10)
        self.assertEqual(f,ct({'a':2.5,'b':2.5}))
        f.scale(100)
        self.assertEqual(f,ct({'a':250, 'b':250}))
        f.scale(0.001)
        self.assertEqual(f,ct({'a':0.25,'b':0.25}))
    
    def test_round(self):
        """Freqs round should round all frequencies to integers"""
        ct = self.ClassToTest
        f = ct({'a':23.1, 'b':12.5, 'c':56.7})
        f.round()
        self.assertEqual(f,ct({'a':23, 'b':13, 'c':57}))
        g = ct({'a':23.1356, 'b':12.5731})
        g.round(3)
        self.assertEqual(g,ct({'a':23.136, 'b':12.573}))
 
    def test_expand(self):
        """Freqs expand should give expected results"""
        ct = self.ClassToTest
        f = ct({'a':3, 'c':5, 'b':2})
        self.assertEqual(f.expand(order='acb'), list('aaacccccbb'))
        self.assertEqual(f.expand(order='dba'), list('bbaaa'))
        self.assertEqual(f.expand(order='cba',convert_to=''.join),'cccccbbaaa')
        f['c'] = 0
        self.assertEqual(f.expand(order='acb'), list('aaabb'))
        f.normalize()
        self.assertEqual(f.expand(order='cba'), ['a'])
        self.assertEqual(f.expand(convert_to=''.join), 'a')
        f.normalize(total=1.0/20)
        self.assertEqual(f.expand(order='abc'), list('a'*12 + 'b'*8))
        #test expand with scaling
        g = ct({'c':0.5,'d':0.5})
        self.assertEqual(g.expand(order='cd'),['c','d'])
        self.assertEqual(g.expand(order='cd',scale=10),list(5*'c'+5*'d'))
        self.assertRaises(ValueError,g.expand,scale=33)

    def test_Count(self):
        """Freqs Count should return correct count (number of categories)"""
        self.assertEqual(self.Alphabetic.Count, 5)
        #WARNING: Count of empty categories is implementation-dependent

    def test_Sum(self):
        """Freqs Sum should return sum of item counts in all categories"""
        self.assertEqual(self.Alphabetic.Sum, 8)
        self.assertEqual(self.Empty.Sum, 0)
        self.assertEqual(self.Constant.Sum, 5)

    def test_SumSquares(self):
        """Freqs SumSquared should return sum of squared freq of each category"""
        self.assertEqual(self.Alphabetic.SumSquares, 16)
        self.assertEqual(self.Empty.SumSquares, 0)
        self.assertEqual(self.Constant.SumSquares, 25)

    def test_Variance(self):
        """Freqs Variance should return variance of counts in categories"""
        self.assertFloatEqual(self.Alphabetic.Variance, 0.8)
        self.assertFloatEqual(self.Empty.Variance, None)
        #WARNING: Variance with empty categories is implementation-dependent

    def test_StandardDeviation(self):
        """Freqs StandardDeviation should return stdev of counts in categories"""
        self.assertFloatEqual(self.Alphabetic.StandardDeviation, sqrt(0.8))
        self.assertFloatEqual(self.Empty.StandardDeviation, None)
        #WARNING: Standard deviation with empty categories is implementation-
        #dependent
       
    def test_Mean(self):
        """Freqs Mean should return mean of counts in categories"""
        self.assertEqual(self.Alphabetic.Mean, 8/5.0)
        self.assertEqual(self.Empty.Mean, None)
        #WARNING: Mean with empty categories is implementation-dependent

    def test_Uncertainty(self):
        """Freqs Shannon uncertainty values should match spreadsheet"""
        self.assertEqual(self.Empty.Uncertainty, 0)
        self.assertFloatEqual(self.Alphabetic.Uncertainty, 2.1556, eps=1e-4)
        #WARNING: Uncertainty with empty categories is implementation-dependent

    def test_mode(self):
        """Freqs mode should return most frequent item"""
        self.assertEqual(self.Empty.Mode, None)
        self.assertEqual(self.Alphabetic.Mode, 'a')
        self.assertEqual(self.Constant.Mode, 'a')

    def test_summarize(self):
        """Freqs summarize should return Summary: Count, Sum, SumSquares, Var"""
        s = self.Alphabetic.summarize()
        self.assertEqual(s.Sum, 8)
        self.assertEqual(s.Count, 5)
        self.assertEqual(s.SumSquares, 16)
        self.assertFloatEqual(s.Variance, 0.8)
        self.assertFloatEqual(s.StandardDeviation, sqrt(0.8))
        self.assertFloatEqual(s.Mean, 8.0/5)
 
    def test_getSortedList(self):
        """Freqs getSortedList should return sorted list of key, val tuples"""
        #behavior is implementation-defined with empty list, so skip tests.
        a = self.Alphabetic
        a['b'] = 5
        self.assertEqual(a.getSortedList(), \
            [('b',5),('a',3),('e',1),('d',1),('c',1)])
        self.assertEqual(a.getSortedList(descending=True), \
            [('b',5),('a',3),('e',1),('d',1),('c',1)])
        self.assertEqual(a.getSortedList(descending=False), \
            [('c',1),('d',1),('e',1),('a',3),('b',5)])
        self.assertEqual(a.getSortedList(by_val=False), \
            [('e',1),('d',1),('c',1),('b',5),('a',3)])
        self.assertEqual(a.getSortedList(by_val=False, descending=False), \
            [('a',3),('b',5),('c',1),('d',1),('e',1)])
 
class FreqsStaticTests(StaticFreqsTestsI, TestCase):
    ClassToTest = Freqs

class UnsafeFreqsStaticTests(StaticFreqsTestsI, TestCase):
    ClassToTest = UnsafeFreqs
        
class FreqsTestsI(object):
    """Tests of the interface shared by Freqs and UnsafeFreqs."""
    ClassToTest = None

    #The following test various ways of constructing the objects

    def test_fromTuples(self):
        """Freqs fromTuples should add from key, count pairs w/ repeated keys"""
        f = self.ClassToTest()
        self.assertEqual(f.fromTuples([('a',4),('b',3),('a',2)]), {'a':6,'b':3})
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative.
        f.fromTuples([('b',-1),('c',4.5)])
        self.assertEqual(f, {'a':6,'b':2,'c':4.5})
        #should work with a different operator
        f.fromTuples([('b',7)], op=mul)
        self.assertEqual(f, {'a':6, 'b':14, 'c':4.5})
        #check that it works with something that depends on the key
        def func(key, first, second):
            if key == 'a':
                return first + second
            else:
                return max(second, first * second)
        f = self.ClassToTest()
        self.assertEqual(f.fromTuples([('a',4),('b',3),('a',2), ('b',4)], \
            func,uses_key=True), {'a':6,'b':12})
        
    
    def test_fromDict(self):
        """Freqs fromDict should add from dict of {key:count}"""
        f = self.ClassToTest()
        self.assertEqual(f.fromDict({'a':6,'b':3}), {'a':6,'b':3})
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative.
        f.fromDict({'b':-1, 'c':4.5})
        self.assertEqual(f, {'a':6,'b':2,'c':4.5})
        #should work with a different operator
        f.fromDict({'b':7}, op=mul)
        self.assertEqual(f, {'a':6, 'b':14, 'c':4.5})
    
    def test_fromDicts(self):
        """Freqs fromDicts should add from list of dicts of {key:count}"""
        f = self.ClassToTest()
        self.assertEqual(f.fromDicts([{'a':6},{'b':3}]), {'a':6,'b':3})
        #note: should be allowed to subtract, as long as freq doesn't go
        #negative. Also tests add of 1-item dict (note: must be in list)
        f.fromDicts([{'b':-1, 'c':4.5}])
        self.assertEqual(f, {'a':6,'b':2,'c':4.5})
        #should work with a different operator
        f.fromDicts([{'b':2},{'b':3}], op=mul)
        self.assertEqual(f, {'a':6, 'b':12, 'c':4.5})

    def test_fromSeq(self):
        """Freqs fromSeq should add items from sequence, according to weight"""
        f = self.ClassToTest()
        self.assertEqual(f.fromSeq('aaabbbaaa'), {'a':6,'b':3})
        #should be able to change the operator...
        self.assertEqual(f.fromSeq('aab', sub), {'a':4,'b':2})
        #...or change the weight
        self.assertEqual(f.fromSeq('acc', weight=3.5), {'a':7.5,'b':2,'c':7})

    def test_fromSeqs(self):
        """Freqs fromSeqs should add items from sequences, according to weight"""
        f = self.ClassToTest()
        self.assertEqual(f.fromSeqs(['aaa','bbbaaa']), {'a':6,'b':3})
        #should be able to change the operator...
        self.assertEqual(f.fromSeqs(list('aab'), sub), {'a':4,'b':2})
        #...or change the weight. Note that a string counts as a seq of seqs.
        self.assertEqual(f.fromSeqs('acc', weight=3.5), {'a':7.5,'b':2,'c':7})

    def test_isValid(self):
        """Freqs isValid should return True if valid"""
        d =self.ClassToTest()
        assert d.isValid()
        d.fromSeq('aaaaaaaaaaaaabb')
        assert d.isValid()

    def test_find_conversion_function(self):
        """Freqs _find_conversion_function should return correct value."""
        d = self.ClassToTest()
        f = d._find_conversion_function
        #should always return None if data empty
        for i in [None, 0, False, {}, [], tuple()]:
            self.assertEqual(f(i), None)
        #should return fromDict for non-empty dict
        self.assertEqual(f({3:4}), d.fromDict)
        #should return fromSeq for string or list of scalars or strings
        for i in ['abc', [1,2,3], (1,2,3), ['ab','bb','cb']]:
            self.assertEqual(f(i), d.fromSeq)
        #should return fromSeqs for sequence of sequences
        for i in [[[1,2,3],[3,4,4]], ([1,2,4],[3,4,4]), [(1,2),[3],[], [4]]]:
            self.assertEqual(f(i), d.fromSeqs)
        #should return fromTuples if possibly key-value pairs
        for i in [[('a',3),('b',-1)], [(1,2),(3,4)]]:
            self.assertEqual(f(i), d.fromTuples)
        #should not be fooled by 2-item seqs that can't be key-value pairs
        self.assertEqual(f(['ab','cd']), d.fromSeq)

    #The following test inheritance of dict properties/methods

    def test_setitem_good(self):
        """Freqs should allow non-negative values to be set"""
        self.Empty[3] = 0
        self.assertEqual(self.Empty, {3:0})
        self.Empty['xyz'] = 5
        self.assertEqual(self.Empty, {3:0, 'xyz':5})

    def test_delitem(self):
        """Freqs should delete all counts of item with del"""
        del self.Alphabetic['a']
        del self.Alphabetic['b']
        del self.Alphabetic['c']
        self.assertEqual(self.Alphabetic, {'d':1,'e':1})

    def test_setdefault_good(self):
        """Freqs setdefault should work with positive values"""
        a = self.Alphabetic.setdefault('a', 200)
        self.assertEqual(a, 3)
        self.assertEqual(self.Alphabetic['a'], 3)

        f = self.Alphabetic.setdefault('f', 0)
        self.assertEqual(f, 0)
        self.assertEqual(self.Alphabetic['f'], 0)

        g = self.Alphabetic.setdefault('g', 1000)
        self.assertEqual(g, 1000)
        self.assertEqual(self.Alphabetic['g'], 1000)

    #The following test overridden operators and methods

    def test_iadd(self):
        """Freqs += should add in place from any known data type"""
        f = self.ClassToTest({'a':3, 'b':4})
        f += 'aca'
        self.assertEqual(f, {'a':5, 'b':4, 'c':1})
        f += ['b','b']
        self.assertEqual(f, {'a':5, 'b':6, 'c':1})
        f += {'c':10, 'a':-3}
        self.assertEqual(f, {'a':2, 'b':6, 'c':11})
        f += (('a',3),('b',-2))
        self.assertEqual(f, {'a':5, 'b':4, 'c':11})
        f += [['a', 'b', 'b'],['c', 'c', 'c']]
        self.assertEqual(f, {'a':6, 'b':6, 'c':14})
        #note that list of strings will use the strings as keys
        f += ['abc', 'def', 'abc']
        self.assertEqual(f, {'a':6, 'b':6, 'c':14, 'abc':2, 'def':1})

    def test_add(self):
        """Freqs + should make new object, adding from any known data type"""
        orig = {'a':3, 'b':4}
        f = self.ClassToTest(orig)
        self.assertEqual(f, orig)
        r = f + 'aca'
        self.assertEqual(r, {'a':5, 'b':4, 'c':1})
        self.assertEqual(f, orig)
        r = f + ['b','b']
        self.assertEqual(r, {'a':3, 'b':6})
        self.assertEqual(f, orig)
        r = f + {'c':10, 'a':-3}
        self.assertEqual(r, {'a':0, 'b':4, 'c':10})
        self.assertEqual(f, orig)
        r = f + (('a',3),('b',-2))
        self.assertEqual(r, {'a':6, 'b':2})
        self.assertEqual(f, orig)
        r = f + [['a', 'b', 'b'],['c', 'c', 'c']]
        self.assertEqual(r, {'a':4, 'b':6, 'c':3})
        self.assertEqual(f, orig)
        #note that list of strings will use the strings as keys
        r = f + ['abc', 'def', 'abc']
        self.assertEqual(r, {'a':3, 'b':4, 'abc':2, 'def':1})
        self.assertEqual(f, f)
        
    def test_isub(self):
        """Freqs -= should subtract in place using any known data type"""
        f = self.ClassToTest({'a':5, 'b':4})
        f -= 'aba'
        self.assertEqual(f, {'a':3, 'b':3,})
        f -= ['b','b']
        self.assertEqual(f, {'a':3, 'b':1})
        f -= {'c':-2, 'a':-3}
        self.assertEqual(f, {'a':6, 'b':1, 'c':2})
        f -= (('a',3),('b',-2))
        self.assertEqual(f, {'a':3, 'b':3, 'c':2})
        f -= [['a', 'b', 'b'],['c', 'c']]
        self.assertEqual(f, {'a':2, 'b':1, 'c':0})
        f['abc'] = 3
        f['def'] = 10
        #note that list of strings will use the strings as keys
        f -= ['abc', 'def', 'abc']
        self.assertEqual(f, {'a':2, 'b':1, 'c':0, 'abc':1, 'def':9})

    def test_sub(self):
        """Freqs - should make new object, subtracting using any known data type"""
        orig = {'a':3, 'b':4}
        f = self.ClassToTest(orig)
        self.assertEqual(f, orig)
        r = f - 'aba'
        self.assertEqual(r, {'a':1, 'b':3})
        self.assertEqual(f, orig)
        r = f - ['b','b']
        self.assertEqual(r, {'a':3, 'b':2})
        self.assertEqual(f, orig)
        r = f - {'c':-10, 'a':3}
        self.assertEqual(r, {'a':0, 'b':4, 'c':10})
        self.assertEqual(f, orig)
        r = f - (('a',3),('b',-2))
        self.assertEqual(r, {'a':0, 'b':6})
        self.assertEqual(f, orig)
        r = f - [['a', 'b', 'b'],['a','a']]
        self.assertEqual(r, {'a':0, 'b':2})
        self.assertEqual(f, orig)
        #note that list of strings will use the strings as keys
        orig['abc'] = 5
        orig['def'] = 10
        f['abc'] = 5
        f['def'] = 10
        r = f - ['abc', 'def', 'abc']
        self.assertEqual(r, {'a':3, 'b':4, 'abc':3, 'def':9})
        self.assertEqual(f, orig)
        
    def test_copy(self):
        """Freqs copy should preserve class of original"""
        d = {'a':4, 'b':3, 'c':6}
        f = self.ClassToTest(d)
        g = f.copy()
        self.assertEqual(d, f)
        self.assertEqual(d, g)
        self.assertEqual(f, g)
        self.assertEqual(f.__class__, g.__class__)

    def test_str(self):
        """Freqs should print as tab-delimited table, or 'Empty'"""
        #should work with empty freq distribution
        self.assertEqual(str(self.ClassToTest([])), \
                "Empty frequency distribution")
        #should work with single element
        self.assertEqual(str(self.ClassToTest({'X':1.0})), \
                "Value\tCount\nX\t1.0")
        #should work with multiples of same key
        self.assertEqual(str(self.ClassToTest({1.0:5.0})), \
                "Value\tCount\n1.0\t5.0")
        #should work with different keys
        self.assertEqual(str(self.ClassToTest({0:3.0,1:2.0})), \
                "Value\tCount\n0\t3.0\n1\t2.0")

    def test_delitem(self):
        """Freqs delitem should refuse to delete a required key"""
        a = self.Alphabetic
        del a['a']
        self.assertEqual(a, {'b':2, 'c':1, 'd':1, 'e':1})
        #can't delete RequiredKeys once set
        a.RequiredKeys = 'bcd'
        del a['e']
        self.assertEqual(a, {'b':2,'c':1,'d':1})
        for k in 'bcd':
            self.assertRaises(KeyError, a.__delitem__, k)
        #when RequiredKeys is removed, can delete them again
        a.RequiredKeys = None
        del a['b']
        self.assertEqual(a, {'c':1,'d':1})

    #The following test custom methods
        
    def test_rekey(self):
        """Freqs rekey should map the results onto new keys."""
        d = self.ClassToTest({'a':3, 'b':5, 'c':6, 'd':7, 'e':1})
        f = d.rekey({'a':'+', 'b':'-', 'c':'+'})
        self.assertEqual(f, {'+':9, '-':5, None:8})
        self.assertEqual(f.__class__, d.__class__)
        #check with explicit constructor and default
        d = self.ClassToTest({'a':3, 'b':5, 'c':6, 'd':7, 'e':1})
        f = d.rekey({'a':'+', 'b':'-', 'c':'+'}, default='x', constructor=dict)
        self.assertEqual(f, {'+':9, '-':5, 'x':8})
        self.assertEqual(f.__class__, dict)

    def test_purge(self):
        """Freqs purge should have no effect unless RequiredKeys set."""
        working = self.PosNeg.copy()
        self.assertEqual(working, self.PosNeg)
        working.purge()
        self.assertEqual(working, self.PosNeg)
        working.RequiredKeys=(-2,-1)
        working[-2] = 3
        working.purge()
        self.assertEqual(working, {-2:3,-1:1})
        #should have no effect if repeated
        working.purge()
        self.assertEqual(working, {-2:3,-1:1})
       
    def test_normalize(self):
        """Freqs should allow normalization on any type"""
        self.Empty.normalize()
        self.assertEqual(self.Empty, {})
        a = self.Alphabetic.copy()
        a.normalize()
        expected = {'a':0.375, 'b':0.25, 'c':0.125, 'd':0.125, 'e':0.125}
        for key, val in expected.items():
            self.assertFloatEqual(a[key], val)
        self.PosNeg.normalize()
        expected = {-2:0.25, -1:0.25, 1:0.25, 2:0.25}
        for key, val in expected.items():
            self.assertFloatEqual(self.PosNeg[key], val)
        #check that it works when we specify a total
        self.PosNeg.normalize(total=0.2)
        expected = {-2:1.25, -1:1.25, 1:1.25, 2:1.25}
        for key, val in expected.items():
            self.assertFloatEqual(self.PosNeg[key], val)
        #check that purging works
        a = self.Alphabetic.copy()
        a.RequiredKeys = 'ac'
        a.normalize()
        self.assertEqual(a, {'a':0.75, 'c':0.25})
        a = self.Alphabetic.copy()
        a.RequiredKeys = 'ac'
        a.normalize(purge=False)
        self.assertEqual(a, \
            {'a':0.375, 'b':0.25, 'c':0.125, 'd':0.125, 'e':0.125})
        #normalize should also create keys when necessary
        a.RequiredKeys = 'bdex'
        a.normalize(purge=True)
        self.assertEqual(a, {'b':0.5, 'd':0.25, 'e':0.25, 'x':0})

    def test_choice(self):
        """Freqs choice should work as expected"""
        self.Alphabetic.normalize()
        keys = self.Alphabetic.keys()
        vals = Numbers(self.Alphabetic.values())
        vals.accumulate()
        #test first item
        self.assertEqual(self.Alphabetic.choice(-1), keys[0])
        self.assertEqual(self.Alphabetic.choice(-0.0001), keys[0])
        self.assertEqual(self.Alphabetic.choice(-1e300), keys[0])
        self.assertEqual(self.Alphabetic.choice(0), keys[0])
        #test last item
        last_val = vals.pop()
        self.assertEqual(self.Alphabetic.choice(last_val), keys[-1])
        self.assertEqual(self.Alphabetic.choice(1000), keys[-1])
        #test remaining items
        for index, value in enumerate(vals):
            self.assertEqual(self.Alphabetic.choice(value-0.01),keys[index])
            self.assertEqual(self.Alphabetic.choice(value+0.01), keys[index+1])
            
    def test_randomSequence_good(self):
        """Freqs randomSequence should give correct counts"""
        self.Alphabetic.normalize()
        total = self.Alphabetic.Sum
        keys = self.Alphabetic.keys()
        probs = [float(i)/total for i in self.Alphabetic.values()]
        
        rand_seq = self.Alphabetic.randomSequence(10000)
        observed = [rand_seq.count(key) for key in keys]
        expected = [prob*10000 for prob in probs]
        self.assertSimilarFreqs(observed, expected)

    def test_randomSequence_bad(self):
        """Empty Freqs should raise error on randomSequence"""
        self.assertRaises(IndexError, self.Empty.randomSequence, 5)

    def test_randomSequence_one_item(self):
        """Freqs randomSequence should work with one key"""
        self.Constant.normalize()
        rand = self.Constant.randomSequence(1000)
        self.assertEqual(rand.count(1), 1000)
        self.assertEqual(len(rand), 1000)

    def test_subset_preserve(self):
        """Freqs subset should preserve unwanted items"""
        self.Constant.subset('abc')
        self.assertEqual(self.Constant, self.Empty)
        self.Alphabetic.subset('abx')
        self.assertEqual(self.Alphabetic, Freqs('aaabb'))

    def test_subset_remove(self):
        """Freqs subset should delete unwanted items if asked"""
        self.Alphabetic.subset('abx', keep=False)
        self.assertEqual(self.Alphabetic, Freqs('cde'))
        self.Constant.subset('abx', keep=False)
        self.assertEqual(self.Constant, Freqs([1]*5))

    def test_scale(self):
        """Freqs scale should multiply all values with the given scale"""
        f = self.ClassToTest({'a':0.25,'b':0.25})
        f.scale(10)
        self.assertEqual(f,{'a':2.5,'b':2.5})
        f.scale(100)
        self.assertEqual(f,{'a':250, 'b':250})
        f.scale(0.001)
        self.assertEqual(f,{'a':0.25,'b':0.25})
    
    def test_round(self):
        """Freqs round should round all frequencies to integers"""
        f = self.ClassToTest({'a':23.1, 'b':12.5, 'c':56.7})
        f.round()
        self.assertEqual(f,{'a':23, 'b':13, 'c':57})
        g = Freqs({'a':23.1356, 'b':12.5731})
        g.round(3)
        self.assertEqual(g,{'a':23.136, 'b':12.573})
 
    def test_expand(self):
        """Freqs expand should give expected results"""
        f = self.ClassToTest({'U':3, 'A':5, 'C':2})
        self.assertEqual(f.expand(order='UAC'), list('UUUAAAAACC'))
        self.assertEqual(f.expand(order='GCU'), list('CCUUU'))
        self.assertEqual(f.expand(order='ACU',convert_to=''.join),'AAAAACCUUU')
        del f['A']
        self.assertEqual(f.expand(order='UAC'), list('UUUCC'))
        f.normalize()
        self.assertEqual(f.expand(order='ACU'), ['U'])
        self.assertEqual(f.expand(convert_to=''.join), 'U')
        f.normalize(total=1.0/20)
        self.assertEqual(f.expand(order='UCA'), list('U'*12 + 'C'*8))
        #test expand with scaling
        g = self.ClassToTest({'A':0.5,'G':0.5})
        self.assertEqual(g.expand(order='AG'),['A','G'])
        self.assertEqual(g.expand(order='AG',scale=10),list(5*'A'+5*'G'))
        self.assertRaises(ValueError,g.expand,scale=33)

    def test_Count(self):
        """Freqs Count should return correct count (number of categories)"""
        self.assertEqual(self.Alphabetic.Count, 5)
        self.assertEqual(self.NumericDuplicated.Count, 3)
        self.assertEqual(self.Empty.Count, 0)
        self.assertEqual(self.Constant.Count, 1)

    def test_Sum(self):
        """Freqs Sum should return sum of item counts in all categories"""
        self.assertEqual(self.Alphabetic.Sum, 8)
        self.assertEqual(self.NumericUnique.Sum, 5)
        self.assertEqual(self.NumericDuplicated.Sum, 4)
        self.assertEqual(self.Empty.Sum, 0)
        # WARNING: For numeric keys, the value of the key is not taken into 
        # account (i.e. each key counts as a separate category)
        self.assertEqual(self.PosNeg.Sum, 4)
        self.assertEqual(self.Constant.Sum, 5)

    def test_SumSquares(self):
        """Freqs SumSquared should return sum of squared freq of each category"""
        self.assertEqual(self.Alphabetic.SumSquares, 16)
        self.assertEqual(self.NumericUnique.SumSquares, 5)
        self.assertEqual(self.NumericDuplicated.SumSquares, 6)
        self.assertEqual(self.Empty.SumSquares, 0)
        self.assertEqual(self.PosNeg.SumSquares, 4)
        self.assertEqual(self.Constant.SumSquares, 25)

    def test_Variance(self):
        """Freqs Variance should return variance of counts in categories"""
        self.assertFloatEqual(self.Alphabetic.Variance, 0.8)
        self.assertFloatEqual(self.NumericUnique.Variance, 0)
        self.assertFloatEqual(self.NumericDuplicated.Variance, 1.0/3)
        self.assertFloatEqual(self.Empty.Variance, None)
        self.assertFloatEqual(self.PosNeg.Variance, 0)
        self.assertEqual(self.Constant.Variance, 0)

    def test_StandardDeviation(self):
        """Freqs StandardDeviation should return stdev of counts in categories"""
        self.assertFloatEqual(self.Alphabetic.StandardDeviation, sqrt(0.8))
        self.assertFloatEqual(self.NumericUnique.StandardDeviation, 0)
        self.assertFloatEqual(self.NumericDuplicated.StandardDeviation, \
                sqrt(1.0/3))
        self.assertFloatEqual(self.Empty.StandardDeviation, None)
        self.assertFloatEqual(self.PosNeg.StandardDeviation, 0)
        self.assertEqual(self.Constant.StandardDeviation, 0)
       
    def test_Mean(self):
        """Freqs Mean should return mean of counts in categories"""
        self.assertEqual(self.Alphabetic.Mean, 8/5.0)
        self.assertEqual(self.NumericUnique.Mean, 1)
        self.assertEqual(self.NumericDuplicated.Mean, 4/3.0)
        self.assertEqual(self.Empty.Mean, None)
        self.assertEqual(self.PosNeg.Mean, 1)
        self.assertEqual(self.Constant.Mean, 5)

    def test_Uncertainty(self):
        """Freqs Shannon uncertainty values should match spreadsheet"""
        self.assertEqual(self.Empty.Uncertainty, 0)
        self.assertFloatEqual(self.Alphabetic.Uncertainty, 2.1556, eps=1e-4)
        self.assertFloatEqual(self.PosNeg.Uncertainty, 2)
        self.assertFloatEqual(self.NumericDuplicated.Uncertainty, 1.5)
        self.assertFloatEqual(self.NumericUnique.Uncertainty, 2.3219, eps=1e-4)
        self.assertFloatEqual(self.Constant.Uncertainty, 0)

    def test_mode(self):
        """Freqs mode should return most frequent item"""
        self.assertEqual(self.Empty.Mode, None)
        self.assertEqual(self.Alphabetic.Mode, 'a')
        assert(self.PosNeg.Mode in self.PosNeg)
        assert(self.NumericUnique.Mode in self.NumericUnique)
        self.assertEqual(self.NumericDuplicated.Mode, 1.5)
        self.assertEqual(self.Constant.Mode, 1)

    def test_summarize(self):
        """Freqs summarize should return Summary: Count, Sum, SumSquares, Var"""
        self.assertEqual(self.Empty.summarize(), SummaryStatistics())
        s = self.Alphabetic.summarize()
        self.assertEqual(s.Sum, 8)
        self.assertEqual(s.Count, 5)
        self.assertEqual(s.SumSquares, 16)
        self.assertFloatEqual(s.Variance, 0.8)
        self.assertFloatEqual(s.StandardDeviation, sqrt(0.8))
        self.assertFloatEqual(s.Mean, 8.0/5)
 
    def test_getSortedList(self):
        """Freqs getSortedList should return sorted list of key, val tuples"""
        e = self.Empty
        self.assertEqual(e.getSortedList(), [])
        self.assertEqual(e.getSortedList(descending=True), [])
        self.assertEqual(e.getSortedList(descending=False), [])
        self.assertEqual(e.getSortedList(by_val=True), [])
        self.assertEqual(e.getSortedList(by_val=False), [])
        a = self.Alphabetic
        a['b'] = 5
        self.assertEqual(a.getSortedList(), \
            [('b',5),('a',3),('e',1),('d',1),('c',1)])
        self.assertEqual(a.getSortedList(descending=True), \
            [('b',5),('a',3),('e',1),('d',1),('c',1)])
        self.assertEqual(a.getSortedList(descending=False), \
            [('c',1),('d',1),('e',1),('a',3),('b',5)])
        self.assertEqual(a.getSortedList(by_val=False), \
            [('e',1),('d',1),('c',1),('b',5),('a',3)])
        self.assertEqual(a.getSortedList(by_val=False, descending=False), \
            [('a',3),('b',5),('c',1),('d',1),('e',1)])
 
class FreqsTests(FreqsTestsI, TestCase):
    """Tests of Freqs-specific behavior, mostly validation."""
    ClassToTest = Freqs
    def setUp(self):
        """defines some standard frequency distributions to check"""
        self.Alphabetic = self.ClassToTest('abcdeaab')
        self.NumericUnique = self.ClassToTest([1,2,3,4,5])
        self.NumericDuplicated = self.ClassToTest([1, 1.5, 1.5, 3.5])
        self.Empty = self.ClassToTest('')
        self.PosNeg = self.ClassToTest([-2, -1, 1, 2])
        self.Constant = self.ClassToTest([1]*5)

    def test_isValid_bad(self):
        """Freqs should reject invalid data, so isValid() always True"""
        self.assertRaises(ConstraintError, self.ClassToTest, {'a':3, 'b':-10})

    def test_init_empty(self):
        """Freqs should initialize OK with empty list"""
        self.assertEqual(self.ClassToTest([]), {})
                         
    def test_init_single(self):
        """Freqs should initialize OK with single item"""
        self.assertEqual(self.ClassToTest(['X']), {'X':1.0})
                        
    def test_init_same_key(self):
        """Freqs should initialize OK with duplicate items"""
        self.assertEqual(self.ClassToTest([1]*5), {1:5})
    
    def test_init_two_keys(self):
        """Freqs should initialize OK with distinct items"""
        self.assertEqual(self.ClassToTest([0,1,0,0,1]), {1:2,0:3})
        
    def test_init_strings(self):
        """Freqs should initialize OK with characters in string"""
        self.assertEqual(self.ClassToTest('zabcz'), {'z':2,'a':1,'b':1,'c':1})
    
    def test_init_fails_negative(self):
        """Freqs init should fail on negative frequencies"""
        self.assertRaises(ConstraintError, self.ClassToTest, {'a':3, 'b':-3})

    def test_init_from_dict(self):
        """Freqs should init OK from dictionary"""
        self.assertEqual(self.ClassToTest({'a':3,'b':2}), {'a':3, 'b':2})

    def test_init_from_dicts(self):
        """Freqs should init OK from list of dicts"""
        self.assertEqual(self.ClassToTest([{'a':1,'b':1}, {'a':2,'b':1}]), \
            {'a':3, 'b':2})

    def test_init_from_strings(self):
        """Freqs should init OK from list of strings"""
        self.assertEqual(self.ClassToTest(['abc','def','abc']), \
            {'abc':2,'def':1})

    def test_init_from_tuples(self):
        """Freqs should init OK from list of key-value pairs"""
        self.assertEqual(self.ClassToTest([('a',3),('b',10),('a',2)]), \
            {'a':5,'b':10})

    def test_init_alphabet_success(self):
        """Freqs should init ok with keys matching alphabet"""
        fd = self.ClassToTest('abc', Constraint='abcd')
        self.assertEqual(fd, {'a':1,'b':1,'c':1})
        self.assertRaises(ConstraintError, fd.setdefault, 'x', 1)
        self.assertRaises(ConstraintError, fd.__setitem__, 'x', 1)
    
    def test_init_alphabet_failure(self):
        """Freqs should fail if keys don't match alphabet"""
        try:
            f = Freqs('abcd', Constraint='abc')
        except ConstraintError:
            pass
        else:
            self.fail()
    
    def test_setitem_bad(self):
        """Freqs should not allow negative values"""
        self.assertRaises(ConstraintError, self.Empty.__setitem__, 'xyz', -0.01)

    def test_setdefault_bad(self):
        """Freqs setdefault should fail if default < 0"""
        self.assertRaises(ConstraintError, self.Empty.setdefault, 'a', -1)
        self.assertRaises(ConstraintError, self.Empty.setdefault, 'a', -0.00001)
        self.assertRaises(ConstraintError, self.Empty.setdefault, 'a', "-1")
        self.assertRaises(ConstraintError, self.Empty.setdefault, 'a', "xxxx")

class UnsafeFreqsTests(FreqsTestsI, TestCase):
    """Tests of UnsafeFreqs-specific behavior, mostly validation."""
    ClassToTest = UnsafeFreqs

    def setUp(self):
        """defines some standard frequency distributions to check"""
        self.Alphabetic = self.ClassToTest({'a':3,'b':2,'c':1,'d':1,'e':1})
        self.NumericUnique = self.ClassToTest({'1':1,'2':1,'3':1,'4':1,'5':1})
        self.NumericDuplicated = self.ClassToTest({1:1, 1.5:2, 3.5:1})
        self.Empty = self.ClassToTest({})
        self.PosNeg = self.ClassToTest({-2:1, -1:1, 1:1, 2:1})
        self.Constant = self.ClassToTest({1:5})

    def test_isValid_bad(self):
        """UnsafeFreqs should allow invalid data, returning False for isValid"""
        d = self.ClassToTest({'a':3, 'b':'x'})
        assert not d.isValid()

    def test_init_empty(self):
        """UnsafeFreqs should initialize OK with empty list"""
        self.assertEqual(self.ClassToTest([]), {})
                         
    def test_init_single(self):
        """UnsafeFreqs init FAILS with single item"""
        self.assertRaises(ValueError, self.ClassToTest, ['X'])
                        
    def test_init_same_key(self):
        """UnsafeFreqs init FAILS with list of items"""
        self.assertRaises(TypeError, self.ClassToTest, [1]*5)
    
    def test_init_strings(self):
        """UnsafeFreqs init FAILS with string"""
        self.assertRaises(ValueError, self.ClassToTest, 'zabcz')
    
    def test_init_negative(self):
        """UnsafeFreqs init should SUCCEED on negative frequencies"""
        self.assertEqual(self.ClassToTest({'a':3, 'b':-3}), {'a':3,'b':-3})

    def test_init_from_dict(self):
        """UnsafeFreqs should init OK from dictionary"""
        self.assertEqual(self.ClassToTest({'a':3,'b':2}), {'a':3, 'b':2})

    def test_init_from_dicts(self):
        """UnsafeFreqs init should init LIKE A DICT from list of dicts"""
        # WARNING: Note the difference between this and Freqs init!
        self.assertEqual(self.ClassToTest([{'a':1,'b':1}, {'a':2,'b':1}]), \
            {'a':'b'})

    def test_init_from_strings(self):
        """UnsafeFreqs init should FAIL from list of strings"""
        self.assertRaises(ValueError, self.ClassToTest, ['abc','def','abc'])

    def test_init_from_tuples(self):
        """UnsafeFreqs should init LIKE A DICT from list of key-value pairs"""
        # WARNING: Note the difference between this and Freqs init!
        self.assertEqual(self.ClassToTest([('a',3),('b',10),('a',2)]), \
            {'a':2,'b':10})

class FreqsSubclassTests(TestCase):
    """Freqs subclassing should work correctly, esp. with RequiredKeys."""
    class BaseFreqs(Freqs):
        RequiredKeys = 'UCAG'

    def test_init(self):
        """Freqs subclass init should add RequiredKeys"""
        b = self.BaseFreqs()
        self.assertEqual(b, {'U':0.0,'C':0.0,'A':0.0,'G':0.0})
        self.assertEqual(self.BaseFreqs('UUCCCCAAAabc'), \
            {'U':2, 'C':4, 'A':3, 'a':1, 'b':1, 'c':1, 'G':0})

    def test_delitem(self):
        """Freqs subclass delitem shouldn't allow deletion of RequiredKeys"""
        b = self.BaseFreqs('AAGCg')
        self.assertEqual(b, {'A':2,'G':1,'C':1,'U':0,'g':1})
        del b['g']
        self.assertEqual(b, {'A':2,'G':1,'C':1,'U':0})
        self.assertRaises(KeyError, b.__delitem__, 'A')

    def test_purge(self):
        """Freqs subclass purge should eliminate anything not in RequiredKeys"""
        b = self.BaseFreqs('AjaknadjkAjnjndfjndCnjdjsfnfdsjkC32478737&#^&@GGGG')
        b.purge()
        self.assertEqual(b, {'A':2,'C':2,'G':4, 'U':0})
        b.purge()
        self.assertEqual(b, {'A':2,'C':2,'G':4, 'U':0})

    def test_normalize(self):
        """Freqs subclass normalize should optionally elminate non-required keys"""
        b = self.BaseFreqs('UUUCX')
        b.normalize(purge=False)
        self.assertEqual(b, {'U':0.6, 'C':0.2, 'X':0.2, 'A':0, 'G':0})
        b.normalize(purge=True)
        self.assertFloatEqual(b, {'U':0.75, 'C':0.25, 'A':0, 'G':0})
        b = self.BaseFreqs()
        b.normalize()
        self.assertEqual(b, {'U':0, 'C':0, 'A':0, 'G':0})


class NumberFreqsTestsI(object):
    """Interface for tests of safe and unsafe NumberFreqs classes."""
    ClassToTest = None
    
    def setUp(self):
        """defines some standard frequency distributions to check"""
        self.NumericUnique = self.ClassToTest([1,2,3,4,5])
        self.NumericDuplicated = self.ClassToTest([1, 1.5, 1.5, 3.5])
        self.Empty = self.ClassToTest('')
        self.PosNeg = self.ClassToTest([-2, -1, 1, 2])
        self.Constant = self.ClassToTest([1]*5)

    def test_setitem_good(self):
        """NumberFreqs should allow non-negative values to be set"""
        self.Empty[3] = 0
        self.assertEqual(self.Empty, {3:0})

    def test_add_good(self):
        """NumberFreqs should allow addition of counts or items"""
        self.Empty += [1]*5
        self.assertEqual(self.Empty, {1:5})

    def test_Mean(self):
        """NumberFreqs means should match hand-calculated values"""
        self.assertEqual(self.Empty.Mean, None)
        self.assertFloatEqual(self.NumericUnique.Mean, 15.0/5)
        self.assertFloatEqual(self.NumericDuplicated.Mean, 7.5/4)
        self.assertFloatEqual(self.PosNeg.Mean, 0.0)
        self.assertFloatEqual(self.Constant.Mean, 1.0)

    def test_Variance(self):
        """NumberFreqs variance should match values from R."""
        self.assertEqual(None, self.Empty.Variance)
        self.assertFloatEqual(2.5, self.NumericUnique.Variance)
        self.assertFloatEqual(1.229167, self.NumericDuplicated.Variance)
        self.assertFloatEqual(10/3.0, self.PosNeg.Variance)
        self.assertFloatEqual(0, self.Constant.Variance)

    def test_Sum(self):
        """NumberFreqs sums should match hand-calculated values"""
        self.assertEqual(self.Empty.Sum, None)
        self.assertFloatEqual(self.NumericUnique.Sum, 15)
        self.assertFloatEqual(self.NumericDuplicated.Sum, 7.5)
        self.assertFloatEqual(self.PosNeg.Sum, 0.0)
        self.assertFloatEqual(self.Constant.Sum, 5.0)

    def test_Count(self):
        """NumberFreqs counts should match hand-calculated values"""
        self.assertEqual(self.NumericUnique.Count, 5)
        self.assertEqual(self.NumericDuplicated.Count, 4)
        self.assertEqual(self.Empty.Count, 0)
        self.assertEqual(self.PosNeg.Count, 4)
        self.assertEqual(self.Constant.Count, 5)

    def test_Sumsq(self):
        """NumberFreqs sum of squares should match spreadsheet"""
        self.assertEqual(self.Empty.SumSquares, None)
        self.assertFloatEqual(self.NumericUnique.SumSquares, 55.0)
        self.assertFloatEqual(self.NumericDuplicated.SumSquares, 17.75)
        self.assertFloatEqual(self.PosNeg.SumSquares, 10.0)
        self.assertFloatEqual(self.Constant.SumSquares, 5.0)

    def test_Stdev(self):
        """NumberFreqs stdev should match spreadsheet"""
        self.assertEqual(self.Empty.StandardDeviation, None)
        self.assertFloatEqual(self.NumericUnique.StandardDeviation,1.581139)
        self.assertFloatEqual(self.NumericDuplicated.StandardDeviation,1.108678)
        self.assertFloatEqual(self.PosNeg.StandardDeviation, 1.825742)
        self.assertFloatEqual(self.Constant.StandardDeviation, 0.0)

    def test_normalize(self):
        """NumberFreqs should allow normalization on any type"""
        self.Empty.normalize()
        self.assertEqual(self.Empty, {})
        # will refuse to normalize if sum is 0
        orig = self.PosNeg.copy()
        self.PosNeg.normalize()
        self.assertEqual(self.PosNeg, orig)
        # will normalize OK if total passed in
        self.PosNeg.normalize(4) #self.PosNeg.Count)
        expected = {-2:0.25, -1:0.25, 1:0.25, 2:0.25}
        for key, val in expected.items():
            self.assertFloatEqual(self.PosNeg[key], val)
        
    def test_Uncertainty(self):
        """NumberFreqs Shannon entropy values should match spreadsheet"""
        self.assertEqual(self.Empty.Uncertainty, 0)
        self.assertEqual(self.PosNeg.Uncertainty, 2)
        self.assertEqual(self.NumericDuplicated.Uncertainty, 1.5)
        self.assertEqual("%6.4f" % self.NumericUnique.Uncertainty, '2.3219')
        self.assertEqual(self.Constant.Uncertainty, 0)

    def test_Mode(self):
        """NumberFreqs mode should return most frequent item"""
        self.assertEqual(self.Empty.Mode, None)
        assert(self.PosNeg.Mode in self.PosNeg)
        assert(self.NumericUnique.Mode in self.NumericUnique)
        self.assertEqual(self.NumericDuplicated.Mode, 1.5)
        self.assertEqual(self.Constant.Mode, 1)

    def test_randomSequence_one_item(self):
        """NumberFreqs randomSequence should work with one key"""
        self.Constant.normalize()
        rand = self.Constant.randomSequence(1000)
        self.assertEqual(rand.count(1), 1000)
        self.assertEqual(len(rand), 1000)

class NumberFreqsTests(NumberFreqsTestsI, TestCase):
    """Tests of (safe) NumberFreqs classes."""
    ClassToTest = NumberFreqs
    
    def setUp(self):
        """defines some standard frequency distributions to check"""
        self.NumericUnique = self.ClassToTest([1,2,3,4,5])
        self.NumericDuplicated = self.ClassToTest([1, 1.5, 1.5, 3.5])
        self.Empty = self.ClassToTest()
        self.PosNeg = self.ClassToTest([-2, -1, 1, 2])
        self.Constant = self.ClassToTest([1]*5)

    def test_setitem_bad(self):
        """NumberFreqs should not allow non-numeric values"""
        self.assertRaises(ValueError, self.Empty.__setitem__, 'xyz', -0.01)

    def test_add_bad(self):
        """NumberFreqs add should fail if key not numeric"""
        self.assertRaises(ValueError, self.Empty.__iadd__, {'a':-1})

class UnsafeNumberFreqsTests(NumberFreqsTestsI, TestCase):
    """Tests of UnsafeNumberFreqs classes."""
    ClassToTest = UnsafeNumberFreqs
    
    def setUp(self):
        """defines some standard frequency distributions to check"""
        self.NumericUnique = self.ClassToTest({1:1,2:1,3:1,4:1,5:1})
        self.NumericDuplicated = self.ClassToTest({1:1,1.5:2,3.5:1})
        self.Empty = self.ClassToTest()
        self.PosNeg = self.ClassToTest({-2:1,-1:1,1:1,2:1})
        self.Constant = self.ClassToTest({1:5})


#execute tests if called from command line
if __name__ == '__main__':
    main()




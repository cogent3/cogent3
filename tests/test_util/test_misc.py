#!/usr/bin/env python

"""Unit tests for utility functions and classes.
"""
from copy import copy, deepcopy
from os import remove
from os.path import exists
from cogent.app.util import get_tmp_filename
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import iterable, max_index, min_index, \
    flatten, is_iterable, is_char, is_char_or_noniterable,\
    is_str_or_noniterable, not_list_tuple, list_flatten,\
    recursive_flatten, unflatten, unzip, select, sort_order, find_all, \
    find_many, unreserve,\
    extract_delimited, caps_from_underscores,\
    add_lowercase, InverseDict, InverseDictMulti, DictFromPos, DictFromFirst, \
    DictFromLast, DistanceFromMatrix, PairsFromGroups, \
    ClassChecker, Delegator, FunctionWrapper, \
    ConstraintError, ConstrainedContainer,\
    ConstrainedString, ConstrainedList, ConstrainedDict, \
    MappedString, MappedList, MappedDict, \
    generateCombinations, makeNonnegInt, \
    NonnegIntError, revComp, not_none, get_items_except,\
    NestedSplitter, curry, app_path, remove_files
from numpy import array

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Amanda Birmingham", "Sandra Smit",
                    "Zongzhi Liu", "Peter Maxwell", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class UtilsTests(TestCase):
    """Tests of individual functions in utils"""

    def test_iterable(self):
        """iterable(x) should return x or [x], always an iterable result"""
        self.assertEqual(iterable('x'), 'x')
        self.assertEqual(iterable(''), '')
        self.assertEqual(iterable(3), [3])
        self.assertEqual(iterable(None), [None])
        self.assertEqual(iterable({'a':1}), {'a':1})
        self.assertEqual(iterable(['a','b','c']), ['a', 'b', 'c'])

    def test_max_index(self):
        """max_index should return index of largest item, last if tie"""
        self.assertEqual(max_index('abcde'), 4)
        self.assertEqual(max_index('ebcda'), 0)
        self.assertRaises(ValueError, max_index, '')
        self.assertEqual(max_index('ebcde'), 4)
        self.assertEqual(max_index([0, 0, 1, 0]), 2)

    def test_min_index(self):
        """min_index should return index of smallest item, first if tie"""
        self.assertEqual(min_index('abcde'), 0)
        self.assertEqual(min_index('ebcda'), 4)
        self.assertRaises(ValueError, min_index, '')
        self.assertEqual(min_index('ebcde'), 1)
        self.assertEqual(min_index([0,0,1,0]), 0)

    def test_flatten_no_change(self):
        """flatten should not change non-nested sequences (except to list)"""
        self.assertEqual(flatten('abcdef'), list('abcdef')) #test identities
        self.assertEqual(flatten([]), []) #test empty sequence
        self.assertEqual(flatten(''), []) #test empty string

    def test_flatten(self):
        """flatten should remove one level of nesting from nested sequences"""
        self.assertEqual(flatten(['aa', 'bb', 'cc']), list('aabbcc'))
        self.assertEqual(flatten([1,[2,3], [[4, [5]]]]), [1, 2, 3, [4,[5]]])

    def test_is_iterable(self):
        """is_iterable should return True for iterables"""
        #test str
        self.assertEqual(is_iterable('aa'), True)
        #test list
        self.assertEqual(is_iterable([3,'aa']), True)
        #test Number, expect False
        self.assertEqual(is_iterable(3), False)
    
    def test_is_char(self):
        """is_char(obj) should return True when obj is a char"""
        self.assertEqual(is_char('a'), True)
        self.assertEqual(is_char('ab'), False)
        self.assertEqual(is_char(''), True)
        self.assertEqual(is_char([3]), False)
        self.assertEqual(is_char(3), False)

    def test_is_char_or_noniterable(self):
        """is_char_or_noniterable should return True or False"""
        self.assertEqual(is_char_or_noniterable('a'), True)
        self.assertEqual(is_char_or_noniterable('ab'), False)
        self.assertEqual(is_char_or_noniterable(3), True)
        self.assertEqual(is_char_or_noniterable([3]), False)

    def test_is_str_or_noniterable(self):
        """is_str_or_noniterable should return True or False"""
        self.assertEqual(is_str_or_noniterable('a'), True)
        self.assertEqual(is_str_or_noniterable('ab'), True)
        self.assertEqual(is_str_or_noniterable(3), True)
        self.assertEqual(is_str_or_noniterable([3]), False)

    def test_recursive_flatten(self):
        """recursive_flatten should remove all nesting from nested sequences"""
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]]), [1,2,3,4,5])

        #test default behavior on str unpacking
        self.assertEqual(recursive_flatten(
            ['aa',[8,'cc','dd'], ['ee',['ff','gg']]]),
            ['a', 'a', 8, 'c', 'c', 'd', 'd', 'e', 'e', 'f', 'f', 'g', 'g'])

        #test str untouched flattening using is_leaf=is_str_or_noniterable
        self.assertEqual(recursive_flatten(
            ['aa',[8,'cc','dd'], ['ee',['ff','gg']]],
            is_leaf=is_str_or_noniterable),
            ['aa',8,'cc','dd','ee','ff','gg'])

    def test_not_list_tuple(self):
        """not_list_tuple(obj) should return False when obj is list or tuple"""
        self.assertEqual(not_list_tuple([8,3]), False)
        self.assertEqual(not_list_tuple((8,3)), False)
        self.assertEqual(not_list_tuple('34'), True)

    def test_list_flatten(self):
        """list_flatten should remove all nesting, str is untouched """
        self.assertEqual(list_flatten(
            ['aa',[8,'cc','dd'], ['ee',['ff','gg']]], ),
            ['aa',8,'cc','dd','ee','ff','gg'])

    def test_recursive_flatten_max_depth(self):
        """recursive_flatten should not remover more than max_depth levels"""
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]]), [1,2,3,4,5])
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]], 0), \
            [1,[2,3], [[4, [5]]]])
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]], 1), \
            [1,2,3, [4, [5]]])
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]], 2), \
            [1,2,3,4, [5]])
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]], 3), \
            [1,2,3,4,5])
        self.assertEqual(recursive_flatten([1,[2,3], [[4, [5]]]], 5000), \
            [1,2,3,4,5])

    def test_unflatten(self):
       """unflatten should turn a 1D sequence into a 2D list"""
       self.assertEqual(unflatten("abcdef", 1), list("abcdef"))
       self.assertEqual(unflatten("abcdef", 1, True), list("abcdef"))
       self.assertEqual(unflatten("abcdef", 2), ['ab','cd','ef'])
       self.assertEqual(unflatten("abcdef", 3), ['abc','def'])
       self.assertEqual(unflatten("abcdef", 4), ['abcd'])
       #should be able to preserve extra items
       self.assertEqual(unflatten("abcdef", 4, True), ['abcd', 'ef'])
       self.assertEqual(unflatten("abcdef", 10), [])
       self.assertEqual(unflatten("abcdef", 10, True), ['abcdef'])
       #should succeed on empty sequnce
       self.assertEqual(unflatten('',10), [])

    def test_unflatten_bad_row_width(self):
        "unflatten should raise ValueError with row_width < 1"""
        self.assertRaises(ValueError, unflatten, "abcd", 0)
        self.assertRaises(ValueError, unflatten, "abcd", -1)

    def test_unzip(self):
        """unzip(items) should be the inverse of zip(*items)"""
        chars = [list('abcde'), list('ghijk')]
        numbers = [[1,2,3,4,5], [0,0,0,0,0]]
        strings = [["abcde", "fghij", "klmno"], ['xxxxx'] * 3]
        empty = [[]]

        lists = [chars, numbers, strings]
        zipped = [zip(*i) for i in lists]
        unzipped = [unzip(i) for i in zipped]

        for u, l in zip(unzipped, lists):
            self.assertEqual(u, l)

    def test_select_sequence(self):
        """select should work on a sequence with a list of indices"""
        chars = 'abcdefghij'
        strings = list(chars)

        tests = {   (0,):['a'],
                    (-1,):['j'],
                    (0, 2, 4): ['a', 'c', 'e'],
                    (9,8,7,6,5,4,3,2,1,0):list('jihgfedcba'),
                    (-8, 8): ['c', 'i'],
                    ():[],
                }
        for test, result in tests.items():
            self.assertEqual(select(test, chars), result)
            self.assertEqual(select(test, strings), result)

    def test_select_empty(self):
        """select should raise error if indexing into empty sequence"""
        self.assertRaises(IndexError, select, [1], [])

    def test_select_mapping(self):
        """select should return the values corresponding to a list of keys"""
        values = {'a':5, 'b':2, 'c':4, 'd':6, 'e':7}
        self.assertEqual(select('abc', values), [5,2,4])
        self.assertEqual(select(['e','e','e'], values), [7,7,7])
        self.assertEqual(select(('e', 'b', 'a'), values), [7, 2, 5])
        #check that it raises KeyError on anything out of range
        self.assertRaises(KeyError, select, 'abx', values)

    def test_sort_order(self):
        """sort_order should return the ordered indices of items"""
        self.assertEqual(sort_order('abc'), [0, 1, 2])
        self.assertEqual(sort_order('cba'), [2,1,0])
        self.assertEqual(sort_order('bca'), [2,0,1])

    def test_sort_order_cmpfunc(self):
        """sort_order should use cmpfunc if passed"""
        self.assertEqual(sort_order([4, 8, 10], lambda x,y:cmp(y,x)), [2, 1, 0])

    def test_sort_order_empty(self):
        """sort_order should return empty list on empty sequence"""
        self.assertEqual(sort_order([]), [])

    def test_find_all(self):
        """find_all should return list of all occurrences"""
        self.assertEqual(find_all('abc', 'd'), [])
        self.assertEqual(find_all('abc', 'a'), [0])
        self.assertEqual(find_all('abcabca', 'a'), [0,3,6])
        self.assertEqual(find_all('abcabca', 'c'), [2,5])
        self.assertEqual(find_all('abcabca', '3'), [])
        self.assertEqual(find_all('abcabca', 'bc'), [1,4])
        self.assertRaises(TypeError, find_all,'abcabca', 3)

    def test_find_many(self):
        """find_many should return list of all occurrences of all items"""
        #should be same as find_all for single chars
        self.assertEqual(find_many('abc', 'd'), [])
        self.assertEqual(find_many('abc', 'a'), [0])
        self.assertEqual(find_many('abcabca', 'a'), [0,3,6])
        self.assertEqual(find_many('abcabca', 'c'), [2,5])
        self.assertEqual(find_many('abcabca', '3'), [])
        #should sort together the items from the two lists
        self.assertEqual(find_many('abcabca', 'bc'), [1,2,4,5])
        #note difference between 2-char string and 1-string list
        self.assertEqual(find_many('abcabca', ['bc']), [1,4])
        self.assertRaises(TypeError, find_many,'abcabca', [3])
 

    def test_unreserve(self):
        """unreserve should trim trailing underscore if present."""
        for i in (None, [], ['x'], 'xyz', '', 'a', '__abc'):
            self.assertEqual(unreserve(i), i)
        self.assertEqual(unreserve('_'), '')
        self.assertEqual(unreserve('class_'), 'class')

    def test_extract_delimited_bad_delimiters(self):
        """extract_delimited should raise error if delimiters identical"""
        self.assertRaises(TypeError, extract_delimited, '|acb|acx', '|','|')

    def test_extract_delimited_missing_right(self):
        """extract_delimited should raise error if right delimiter missing"""
        self.assertRaises(ValueError, extract_delimited, 'ac[acgsd', '[', ']')

    def test_extract_delimited_normal(self):
        """extract_delimited should return correct field if present, or None"""
        self.assertEqual(extract_delimited('[]', '[', ']'), '')
        self.assertEqual(extract_delimited('asdsad', '[', ']'), None)
        self.assertEqual(extract_delimited('ac[abc]ac', '[', ']'), 'abc')
        self.assertEqual(extract_delimited('[xyz]asd', '[', ']'), 'xyz')
        self.assertEqual(extract_delimited('acg[xyz]', '[', ']'), 'xyz')
        self.assertEqual(extract_delimited('abcdef', 'a', 'e'), 'bcd')

    def test_extract_delimited_indexed(self):
        """extract_delimited should return correct field with starting index"""
        self.assertEqual(extract_delimited('[abc][def]', '[',']', 0), 'abc')
        self.assertEqual(extract_delimited('[abc][def]','[',']',1), 'def')
        self.assertEqual(extract_delimited('[abc][def]', '[',']',5), 'def')

    def test_caps_from_underscores(self):
        """caps_from_underscores should become CapsFromUnderscores"""
        cfu = caps_from_underscores
        #should still convert strings without underscores
        self.assertEqual(cfu('ABCDE  abcde!$'), 'Abcde  Abcde!$')
        self.assertEqual(cfu('abc_def'), 'AbcDef')
        #should read through multiple underscores
        self.assertEqual(cfu('_caps__from_underscores___'), 
            'CapsFromUnderscores')

    def test_add_lowercase(self):
        """add_lowercase should add lowercase version of each key w/ same val"""
        d = {'a':1, 'b':'test', 'A':5, 'C':123, 'D':[], 'AbC':'XyZ', \
            None:'3', '$':'abc', 145:'5'}
        add_lowercase(d)
        assert d['d'] is d['D']
        d['D'].append(3)
        self.assertEqual(d['D'], [3])
        self.assertEqual(d['d'], [3])
        self.assertNotEqual(d['a'], d['A'])
        self.assertEqual(d, {'a':1, 'b':'test', 'A':5, 'C':123, 'c':123, \
            'D':[3], 'd':[3], 'AbC':'XyZ', 'abc':'xyz', None:'3', '$':'abc', \
            145:'5'})

        #should work with strings
        d = 'ABC'
        self.assertEqual(add_lowercase(d), 'ABCabc')
        #should work with tuples
        d = tuple('ABC')
        self.assertEqual(add_lowercase(d), tuple('ABCabc'))
        #should work with lists
        d = list('ABC')
        self.assertEqual(add_lowercase(d), list('ABCabc'))
        #should work with sets
        d = set('ABC')
        self.assertEqual(add_lowercase(d), set('ABCabc'))
        #...even frozensets
        d = frozenset('ABC')
        self.assertEqual(add_lowercase(d), frozenset('ABCabc'))

    def test_add_lowercase_tuple(self):
        """add_lowercase should deal with tuples correctly"""
        d = {('A','B'):'C', ('D','e'):'F', ('b','c'):'H'}
        add_lowercase(d)
        self.assertEqual(d, {
            ('A','B'):'C',
            ('a','b'):'c',
            ('D','e'):'F',
            ('d','e'):'f',
            ('b','c'):'H',
            })

    def test_InverseDict(self):
        """InverseDict should invert dict's keys and values"""
        self.assertEqual(InverseDict({}), {})
        self.assertEqual(InverseDict({'3':4}), {4:'3'})
        self.assertEqual(InverseDict({'a':'x','b':1,'c':None,'d':('a','b')}), \
            {'x':'a',1:'b',None:'c',('a','b'):'d'})
        self.assertRaises(TypeError, InverseDict, {'a':['a','b','c']})
        d = InverseDict({'a':3, 'b':3, 'c':3})
        self.assertEqual(len(d), 1)
        assert 3 in d
        assert d[3] in 'abc'
        
    def test_InverseDictMulti(self):
        """InverseDictMulti should invert keys and values, keeping all keys"""
        self.assertEqual(InverseDictMulti({}), {})
        self.assertEqual(InverseDictMulti({'3':4}), {4:['3']})
        self.assertEqual(InverseDictMulti(\
            {'a':'x','b':1,'c':None,'d':('a','b')}), \
            {'x':['a'],1:['b'],None:['c'],('a','b'):['d']})
        self.assertRaises(TypeError, InverseDictMulti, {'a':['a','b','c']})
        d = InverseDictMulti({'a':3, 'b':3, 'c':3, 'd':'3', 'e':'3'})
        self.assertEqual(len(d), 2)
        assert 3 in d
        d3_items = d[3][:]
        self.assertEqual(len(d3_items), 3)
        d3_items.sort()
        self.assertEqual(''.join(d3_items), 'abc')
        assert '3' in d
        d3_items = d['3'][:]
        self.assertEqual(len(d3_items), 2)
        d3_items.sort()
        self.assertEqual(''.join(d3_items), 'de') 

    def test_DictFromPos(self):
        """DictFromPos should return correct lists of positions"""
        d = DictFromPos
        self.assertEqual(d(''), {})
        self.assertEqual(d('a'), {'a':[0]})
        self.assertEqual(d(['a','a','a']), {'a':[0,1,2]})
        self.assertEqual(d('abacdeeee'), {'a':[0,2],'b':[1],'c':[3],'d':[4], \
            'e':[5,6,7,8]})
        self.assertEqual(d(('abc',None, 'xyz', None, 3)), {'abc':[0],None:[1,3],
            'xyz':[2], 3:[4]})

    def test_DictFromFirst(self):
        """DictFromFirst should return correct first positions"""
        d = DictFromFirst
        self.assertEqual(d(''), {})
        self.assertEqual(d('a'), {'a':0})
        self.assertEqual(d(['a','a','a']), {'a':0})
        self.assertEqual(d('abacdeeee'), {'a':0,'b':1,'c':3,'d':4,'e':5})
        self.assertEqual(d(('abc',None, 'xyz', None, 3)), {'abc':0,None:1,
            'xyz':2, 3:4})

    def test_DictFromLast(self):
        """DictFromLast should return correct last positions"""
        d = DictFromLast
        self.assertEqual(d(''), {})
        self.assertEqual(d('a'), {'a':0})
        self.assertEqual(d(['a','a','a']), {'a':2})
        self.assertEqual(d('abacdeeee'), {'a':2,'b':1,'c':3,'d':4,'e':8})
        self.assertEqual(d(('abc',None, 'xyz', None, 3)), {'abc':0,None:3,
            'xyz':2, 3:4})

    def test_DistanceFromMatrix(self):
        """DistanceFromMatrix should return correct elements"""
        m = {'a':{'3':4, 6:1}, 'b':{'3':5,'6':2}}
        d = DistanceFromMatrix(m)
        self.assertEqual(d('a','3'), 4)
        self.assertEqual(d('a',6), 1)
        self.assertEqual(d('b','3'), 5)
        self.assertEqual(d('b','6'), 2)
        self.assertRaises(KeyError, d, 'c', 1)
        self.assertRaises(KeyError, d, 'b', 3)

    def test_PairsFromGroups(self):
        """PairsFromGroups should return dict with correct pairs"""
        empty = []
        self.assertEqual(PairsFromGroups(empty), {})
        one = ['abc']
        self.assertEqual(PairsFromGroups(one), dict.fromkeys([ \
            ('a','a'), ('a','b'), ('a','c'), \
            ('b','a'), ('b','b'), ('b','c'), \
            ('c','a'), ('c','b'), ('c','c'), \
            ]))

        two = ['xy', 'abc']
        self.assertEqual(PairsFromGroups(two), dict.fromkeys([ \
            ('a','a'), ('a','b'), ('a','c'), \
            ('b','a'), ('b','b'), ('b','c'), \
            ('c','a'), ('c','b'), ('c','c'), \
            ('x','x'), ('x','y'), ('y','x'), ('y','y'), \
            ]))
        #if there's overlap, note that the groups should _not_ be expanded
        #(e.g. in the following case, 'x' is _not_ similar to 'c', even though
        #both 'x' and 'c' are similar to 'a'.
        overlap = ['ax', 'abc']
        self.assertEqual(PairsFromGroups(overlap), dict.fromkeys([ \
            ('a','a'), ('a','b'), ('a','c'), \
            ('b','a'), ('b','b'), ('b','c'), \
            ('c','a'), ('c','b'), ('c','c'), \
            ('x','x'), ('x','a'), ('a','x'), \
            ]))
            
    def test_remove_files(self):
        """Remove files functions as expected """
        # create list of temp file paths
        test_filepaths = \
         [get_tmp_filename(prefix='remove_files_test') for i in range(5)]
        
        # try to remove them with remove_files and verify that an IOError is 
        # raises
        self.assertRaises(OSError,remove_files,test_filepaths)
        # now get no error when error_on_missing=False
        remove_files(test_filepaths,error_on_missing=False)
        
        # touch one of the filepaths so it exists
        open(test_filepaths[2],'w').close()
        # check that an error is raised on trying to remove the files...
        self.assertRaises(OSError,remove_files,test_filepaths)
        # ... but that the existing file was still removed
        self.assertFalse(exists(test_filepaths[2]))
        
        # touch one of the filepaths so it exists
        open(test_filepaths[2],'w').close()
        # no error is raised on trying to remove the files 
        # (although 4 don't exist)...
        remove_files(test_filepaths,error_on_missing=False)
        # ... and the existing file was removed
        self.assertFalse(exists(test_filepaths[2]))
        
            

class _my_dict(dict):
    """Used for testing subclass behavior of ClassChecker"""
    pass

class ClassCheckerTests(TestCase):
    """Unit tests for the ClassChecker class."""

    def setUp(self):
        """define a few standard checkers"""
        self.strcheck = ClassChecker(str)
        self.intcheck = ClassChecker(int, long)
        self.numcheck = ClassChecker(float, int, long)
        self.emptycheck = ClassChecker()
        self.dictcheck = ClassChecker(dict)
        self.mydictcheck = ClassChecker(_my_dict)
        
    def test_init_good(self):
        """ClassChecker should init OK when initialized with classes"""
        self.assertEqual(self.strcheck.Classes, [str])
        self.assertEqual(self.numcheck.Classes, [float, int, long])
        self.assertEqual(self.emptycheck.Classes, [])

    def test_init_bad(self):
        """ClassChecker should raise TypeError if initialized with non-class"""
        self.assertRaises(TypeError, ClassChecker, 'x')
        self.assertRaises(TypeError, ClassChecker, str, None)

    def test_contains(self):
        """ClassChecker should return True only if given instance of class"""
        self.assertEqual(self.strcheck.__contains__('3'), True)
        self.assertEqual(self.strcheck.__contains__('ahsdahisad'), True)
        self.assertEqual(self.strcheck.__contains__(3), False)
        self.assertEqual(self.strcheck.__contains__({3:'c'}), False)

        self.assertEqual(self.intcheck.__contains__('ahsdahisad'), False)
        self.assertEqual(self.intcheck.__contains__('3'), False)
        self.assertEqual(self.intcheck.__contains__(3.0), False)
        self.assertEqual(self.intcheck.__contains__(3), True)
        self.assertEqual(self.intcheck.__contains__(4**60), True)
        self.assertEqual(self.intcheck.__contains__(4**60 * -1), True)

        d = _my_dict()
        self.assertEqual(self.dictcheck.__contains__(d), True)
        self.assertEqual(self.dictcheck.__contains__({'d':1}), True)
        self.assertEqual(self.mydictcheck.__contains__(d), True)
        self.assertEqual(self.mydictcheck.__contains__({'d':1}), False)

        self.assertEqual(self.emptycheck.__contains__('d'), False)

        self.assertEqual(self.numcheck.__contains__(3), True)
        self.assertEqual(self.numcheck.__contains__(3.0), True)
        self.assertEqual(self.numcheck.__contains__(-3), True)
        self.assertEqual(self.numcheck.__contains__(-3.0), True)
        self.assertEqual(self.numcheck.__contains__(3e-300), True)
        self.assertEqual(self.numcheck.__contains__(0), True)
        self.assertEqual(self.numcheck.__contains__(4**1000), True)
        self.assertEqual(self.numcheck.__contains__('4**1000'), False)

    def test_str(self):
        """ClassChecker str should be the same as str(self.Classes)"""
        for c in [self.strcheck, self.intcheck, self.numcheck, self.emptycheck,
            self.dictcheck, self.mydictcheck]:
            self.assertEqual(str(c), str(c.Classes))

    def test_copy(self):
        """copy.copy should work correctly on ClassChecker"""
        c = copy(self.strcheck)
        assert c is not self.strcheck
        assert '3' in c
        assert 3 not in c
        assert c.Classes is self.strcheck.Classes

    def test_deepcopy(self):
        """copy.deepcopy should work correctly on ClassChecker"""
        c = deepcopy(self.strcheck)
        assert c is not self.strcheck
        assert '3' in c
        assert 3 not in c
        assert c.Classes is not self.strcheck.Classes

class modifiable_string(str):
    """Mutable class to allow arbitrary attributes to be set"""
    pass

class _list_and_string(list, Delegator):
    """Trivial class to demonstrate Delegator.
    """
    def __init__(self, items, string):
        Delegator.__init__(self, string)
        self.NormalAttribute = 'default'
        self._x = None
        self._constant = 'c'
        for i in items:
            self.append(i)

    def _get_rand_property(self):
        return self._x

    def _set_rand_property(self, value):
        self._x = value
    prop = property(_get_rand_property, _set_rand_property)

    def _get_constant_property(self):
        return self._constant
    constant = property(_get_constant_property)

class DelegatorTests(TestCase):
    """Verify that Delegator works with attributes and properties."""

    def test_init(self):
        """Delegator should init OK when data supplied"""
        ls = _list_and_string([1,2,3], 'abc')
        self.assertRaises(TypeError, _list_and_string, [123])
        
    def test_getattr(self):
        """Delegator should find attributes in correct places"""
        ls = _list_and_string([1,2,3], 'abcd')
        #behavior as list
        self.assertEqual(len(ls), 3)
        self.assertEqual(ls[0], 1)
        ls.reverse()
        self.assertEqual(ls, [3,2,1])
        #behavior as string
        self.assertEqual(ls.upper(), 'ABCD')
        self.assertEqual(len(ls.upper()), 4)
        self.assertEqual(ls.replace('a', 'x'), 'xbcd')
        #behavior of normal attributes
        self.assertEqual(ls.NormalAttribute, 'default')
        #behavior of properties
        self.assertEqual(ls.prop, None)
        self.assertEqual(ls.constant, 'c')
        #shouldn't be allowed to get unknown properties
        self.assertRaises(AttributeError, getattr, ls, 'xyz')
        #if the unknown property can be set in the forwarder, do it there
        flex = modifiable_string('abcd')
        ls_flex = _list_and_string([1,2,3], flex)
        ls_flex.blah = 'zxc'
        self.assertEqual(ls_flex.blah, 'zxc')
        self.assertEqual(flex.blah, 'zxc')
        #should get AttributeError if changing a read-only property
        self.assertRaises(AttributeError, setattr, ls, 'constant', 'xyz')
        
    
    def test_setattr(self):
        """Delegator should set attributes in correct places"""
        ls = _list_and_string([1,2,3], 'abcd')
        #ability to create a new attribute
        ls.xyz = 3
        self.assertEqual(ls.xyz, 3)
        #modify a normal attribute
        ls.NormalAttribute = 'changed'
        self.assertEqual(ls.NormalAttribute, 'changed')
        #modify a read/write property
        ls.prop = 'xyz'
        self.assertEqual(ls.prop, 'xyz')

    def test_copy(self):
        """copy.copy should work correctly on Delegator"""
        l = ['a']
        d = Delegator(l)
        c = copy(d)
        assert c is not d
        assert c._handler is d._handler

    def test_deepcopy(self):
        """copy.deepcopy should work correctly on Delegator"""
        l = ['a']
        d = Delegator(l)
        c = deepcopy(d)
        assert c is not d
        assert c._handler is not d._handler
        assert c._handler == d._handler

        

class FunctionWrapperTests(TestCase):
    """Tests of the FunctionWrapper class"""
    def test_init(self):
        """FunctionWrapper should initialize with any callable"""
        f = FunctionWrapper(str)
        g = FunctionWrapper(id)
        h = FunctionWrapper(iterable)
        x = 3
        self.assertEqual(f(x), '3')
        self.assertEqual(g(x), id(x))
        self.assertEqual(h(x), [3])

    def test_copy(self):
        """copy should work for FunctionWrapper objects"""
        f = FunctionWrapper(str)
        c = copy(f)
        assert c is not f
        assert c.Function is f.Function

    #NOTE: deepcopy does not work for FunctionWrapper objects because you
    #can't copy a function.


class _simple_container(object):
    """example of a container to constrain"""
    def __init__(self, data):
        self._data = list(data)
    def __getitem__(self, item):
        return self._data.__getitem__(item)

class _constrained_simple_container(_simple_container, ConstrainedContainer):
    """constrained version of _simple_container"""
    def __init__(self, data):
        _simple_container.__init__(self, data)
        ConstrainedContainer.__init__(self, None)

class ConstrainedContainerTests(TestCase):
    """Tests of the generic ConstrainedContainer interface."""
    def setUp(self):
        """Make a couple of standard containers"""
        self.alphabet = _constrained_simple_container('abc')
        self.numbers = _constrained_simple_container([1,2,3])
        self.alphacontainer = 'abcdef'
        self.numbercontainer = ClassChecker(int)
        
    def test_matchesConstraint(self):
        """ConstrainedContainer matchesConstraint should return true if items ok"""
        self.assertEqual(self.alphabet.matchesConstraint(self.alphacontainer), \
            True)
        self.assertEqual(self.alphabet.matchesConstraint(self.numbercontainer),\
            False)
        self.assertEqual(self.numbers.matchesConstraint(self.alphacontainer), \
            False)
        self.assertEqual(self.numbers.matchesConstraint(self.numbercontainer),\
            True)

    def test_otherIsValid(self):
        """ConstrainedContainer should use constraint for checking other"""
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), False)
        self.alphabet.Constraint = list('abcdefghijkl12345678')
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), True)
        self.assertEqual(self.alphabet.otherIsValid('z'), False)

    def test_itemIsValid(self):
        """ConstrainedContainer should use constraint for checking item"""
        self.assertEqual(self.alphabet.itemIsValid(3), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.itemIsValid(3), False)
        self.assertEqual(self.alphabet.itemIsValid('a'), True)

    def test_sequenceIsValid(self):
        """ConstrainedContainer should use constraint for checking sequence"""
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), False)
        self.alphabet.Constraint = list('abcdefghijkl12345678')
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), True)
        self.assertEqual(self.alphabet.sequenceIsValid('z'), False)

    def test_Constraint(self):
        """ConstrainedContainer should only allow valid constraints to be set"""
        try:
            self.alphabet.Constraint = self.numbers
        except ConstraintError:
            pass
        else:
            raise AssertionError, \
            "Failed to raise ConstraintError with invalid constraint."
        self.alphabet.Constraint = 'abcdefghi'
        self.alphabet.Constraint = ['a','b', 'c', 1, 2, 3]
        self.numbers.Constraint = range(20)
        self.numbers.Constraint = xrange(20)
        self.numbers.Constraint = [5,1,3,7,2]
        self.numbers.Constraint = {1:'a',2:'b',3:'c'}
        self.assertRaises(ConstraintError, setattr, self.numbers, \
            'Constraint', '1')
            
class ConstrainedStringTests(TestCase):
    """Tests that ConstrainedString can only contain allowed items."""
    
    def test_init_good_data(self):
        """ConstrainedString should init OK if string matches constraint"""
        self.assertEqual(ConstrainedString('abc', 'abcd'), 'abc')
        self.assertEqual(ConstrainedString('', 'abcd'), '')
        items = [1,2,3.2234, (['a'], ['b'],), 'xyz']
        #should accept anything str() does if no constraint is passed
        self.assertEqual(ConstrainedString(items), str(items))
        self.assertEqual(ConstrainedString(items, None), str(items))
        self.assertEqual(ConstrainedString('12345'), str(12345))
        self.assertEqual(ConstrainedString(12345, '1234567890'), str(12345))
        #check that list is formatted correctly and chars are all there
        test_list = [1,2,3,4,5]
        self.assertEqual(ConstrainedString(test_list, '][, 12345'), str(test_list))

    def test_init_bad_data(self):
        """ConstrainedString should fail init if unknown chars in string"""
        self.assertRaises(ConstraintError, ConstrainedString, 1234, '123')
        self.assertRaises(ConstraintError, ConstrainedString, '1234', '123')
        self.assertRaises(ConstraintError, ConstrainedString, [1,2,3], '123')

    def test_add_prevents_bad_data(self):
        """ConstrainedString should allow addition only of compliant string"""
        a = ConstrainedString('123', '12345')
        b = ConstrainedString('444', '4')
        c = ConstrainedString('45', '12345')
        d = ConstrainedString('x')
        self.assertEqual(a + b, '123444')
        self.assertEqual(a + c, '12345')
        self.assertRaises(ConstraintError, b.__add__, c)
        self.assertRaises(ConstraintError, c.__add__, d)
        #should be OK if constraint removed
        b.Constraint = None
        self.assertEqual(b + c, '44445')
        self.assertEqual(b + d, '444x')
        #should fail if we add the constraint back
        b.Constraint = '4x'
        self.assertEqual(b + d, '444x')
        self.assertRaises(ConstraintError, b.__add__, c)
        #check that added strings retain constraint
        self.assertRaises(ConstraintError, (a+b).__add__, d)
   
    def test_mul(self):
        """ConstrainedString mul amd rmul should retain constraint"""
        a = ConstrainedString('123', '12345')
        b = 3*a
        c = b*2
        self.assertEqual(b, '123123123')
        self.assertEqual(c, '123123123123123123')
        self.assertRaises(ConstraintError, b.__add__, 'x')
        self.assertRaises(ConstraintError, c.__add__, 'x')
   
    def test_getslice(self):
        """ConstrainedString getslice should remember constraint"""
        a = ConstrainedString('123333', '12345')
        b = a[2:4]
        self.assertEqual(b, '33')
        self.assertEqual(b.Constraint, '12345')

    def test_getitem(self):
        """ConstrainedString getitem should handle slice objects"""
        a = ConstrainedString('7890543', '1234567890')
        self.assertEqual(a[0], '7')
        self.assertEqual(a[1], '8')
        self.assertEqual(a[-1], '3')
        self.assertRaises(AttributeError, getattr, a[1], 'Alphabet')
        self.assertEqual(a[1:6:2], '804')
        self.assertEqual(a[1:6:2].Constraint, '1234567890')

    def test_init_masks(self):
        """ConstrainedString should init OK with masks"""
        def mask(x):
            return str(int(x) + 3)
        a = ConstrainedString('12333', '45678', mask)
        self.assertEqual(a, '45666')
        assert 'x' not in a
        self.assertRaises(TypeError, a.__contains__, 1)

class MappedStringTests(TestCase):
    """MappedString should behave like ConstrainedString, but should map items."""
    def test_init_masks(self):
        """MappedString should init OK with masks"""
        def mask(x):
            return str(int(x) + 3)
        a = MappedString('12333', '45678', mask)
        self.assertEqual(a, '45666')
        assert 1 in a
        assert 'x' not in a

    

class ConstrainedListTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedLists."""
    def test_init_good_data(self):
        """ConstrainedList should init OK if list matches constraint"""
        self.assertEqual(ConstrainedList('abc', 'abcd'), list('abc'))
        self.assertEqual(ConstrainedList('', 'abcd'), list(''))
        items = [1,2,3.2234, (['a'], ['b'],), list('xyz')]
        #should accept anything str() does if no constraint is passed
        self.assertEqual(ConstrainedList(items), items)
        self.assertEqual(ConstrainedList(items, None), items)
        self.assertEqual(ConstrainedList('12345'), list('12345'))
        #check that list is formatted correctly and chars are all there
        test_list = list('12345')
        self.assertEqual(ConstrainedList(test_list, '12345'), test_list)

    def test_init_bad_data(self):
        """ConstrainedList should fail init with items not in constraint"""
        self.assertRaises(ConstraintError, ConstrainedList, '1234', '123')
        self.assertRaises(ConstraintError,ConstrainedList,[1,2,3],['1','2','3'])

    def test_add_prevents_bad_data(self):
        """ConstrainedList should allow addition only of compliant data"""
        a = ConstrainedList('123', '12345')
        b = ConstrainedList('444', '4')
        c = ConstrainedList('45', '12345')
        d = ConstrainedList('x')
        self.assertEqual(a + b, list('123444'))
        self.assertEqual(a + c, list('12345'))
        self.assertRaises(ConstraintError, b.__add__, c)
        self.assertRaises(ConstraintError, c.__add__, d)
        #should be OK if constraint removed
        b.Constraint = None
        self.assertEqual(b + c, list('44445'))
        self.assertEqual(b + d, list('444x'))
        #should fail if we add the constraint back
        b.Constraint = {'4':1, 5:2}
        self.assertRaises(ConstraintError, b.__add__, c)
                
    
    def test_iadd_prevents_bad_data(self):
        """ConstrainedList should allow in-place addition only of compliant data"""
        a = ConstrainedList('12', '123')
        a += '2'
        self.assertEqual(a, list('122'))
        self.assertEqual(a.Constraint, '123')
        self.assertRaises(ConstraintError, a.__iadd__, '4')
    
    def test_imul(self):
        """ConstrainedList imul should preserve constraint"""
        a = ConstrainedList('12', '123')
        a *= 3
        self.assertEqual(a, list('121212'))
        self.assertEqual(a.Constraint, '123')

    def test_mul(self):
        """ConstrainedList mul should preserve constraint"""
        a = ConstrainedList('12', '123')
        b = a * 3
        self.assertEqual(b, list('121212'))
        self.assertEqual(b.Constraint, '123')

    def test_rmul(self):
        """ConstrainedList rmul should preserve constraint"""
        a = ConstrainedList('12', '123')
        b = 3 * a
        self.assertEqual(b, list('121212'))
        self.assertEqual(b.Constraint, '123')

    def test_setitem(self):
        """ConstrainedList setitem should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a[0] = '3'
        self.assertEqual(a, list('32'))
        self.assertRaises(ConstraintError, a.__setitem__, 0, 3)
        a = ConstrainedList('1'*20, '123')
        self.assertRaises(ConstraintError, a.__setitem__, slice(0,1,1), [3])
        self.assertRaises(ConstraintError, a.__setitem__, slice(0,1,1), ['111'])
        a[2:9:2] = '3333'
        self.assertEqual(a, list('11313131311111111111'))

    def test_append(self):
        """ConstrainedList append should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a.append('3')
        self.assertEqual(a, list('123'))
        self.assertRaises(ConstraintError, a.append, 3)

    def test_extend(self):
        """ConstrainedList extend should work only if all items in constraint"""
        a = ConstrainedList('12', '123')
        a.extend('321')
        self.assertEqual(a, list('12321'))
        self.assertRaises(ConstraintError, a.extend, ['1','2', 3])

    def test_insert(self):
        """ConstrainedList insert should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a.insert(0, '2')
        self.assertEqual(a, list('212'))
        self.assertRaises(ConstraintError, a.insert, 0, [2])

    def test_getslice(self):
        """ConstrainedList getslice should remember constraint"""
        a = ConstrainedList('123333', '12345')
        b = a[2:4]
        self.assertEqual(b, list('33'))
        self.assertEqual(b.Constraint, '12345')

    def test_setslice(self):
        """ConstrainedList setslice should fail if slice has invalid chars"""
        a = ConstrainedList('123333', '12345')
        a[2:4] = ['2','2']
        self.assertEqual(a, list('122233'))
        self.assertRaises(ConstraintError, a.__setslice__, 2,4, [2,2])
        a[:] = []
        self.assertEqual(a, [])
        self.assertEqual(a.Constraint, '12345')

    def test_setitem_masks(self):
        """ConstrainedList setitem with masks should transform input"""
        a = ConstrainedList('12333', range(5), lambda x: int(x) + 1)
        self.assertEqual(a, [2,3,4,4,4])
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.Mask is a.Mask
        assert '1' not in a
        assert '2' not in a
        assert 2 in a
        assert 'x' not in a

class MappedListTests(TestCase):
    """MappedList should behave like ConstrainedList, but map items."""
    def test_setitem_masks(self):
        """MappedList setitem with masks should transform input"""
        a = MappedList('12333', range(5), lambda x: int(x) + 1)
        self.assertEqual(a, [2,3,4,4,4])
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.Mask is a.Mask
        assert '1' in a
        assert 'x' not in a

class ConstrainedDictTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedDicts."""
    def test_init_good_data(self):
        """ConstrainedDict should init OK if list matches constraint"""
        self.assertEqual(ConstrainedDict(dict.fromkeys('abc'), 'abcd'), \
            dict.fromkeys('abc'))
        self.assertEqual(ConstrainedDict('', 'abcd'), dict(''))
        items = [1,2,3.2234, tuple('xyz')]
        #should accept anything dict() does if no constraint is passed
        self.assertEqual(ConstrainedDict(dict.fromkeys(items)), \
            dict.fromkeys(items))
        self.assertEqual(ConstrainedDict(dict.fromkeys(items), None), \
            dict.fromkeys(items))
        self.assertEqual(ConstrainedDict([(x,1) for x in '12345']), \
            dict.fromkeys('12345', 1))
        #check that list is formatted correctly and chars are all there
        test_dict = dict.fromkeys('12345')
        self.assertEqual(ConstrainedDict(test_dict, '12345'), test_dict)

    def test_init_sequence(self):
        """ConstrainedDict should init from sequence, unlike normal dict"""
        self.assertEqual(ConstrainedDict('abcda'), {'a':2,'b':1,'c':1,'d':1})

    def test_init_bad_data(self):
        """ConstrainedDict should fail init with items not in constraint"""
        self.assertRaises(ConstraintError, ConstrainedDict, \
            dict.fromkeys('1234'), '123')
        self.assertRaises(ConstraintError,ConstrainedDict, \
        dict.fromkeys([1,2,3]),['1','2','3'])
   
    def test_setitem(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        a['1'] = '3'
        self.assertEqual(a, {'1':'3','2':None})
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')

    def test_copy(self):
        """ConstrainedDict copy should retain constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        b = a.copy()
        self.assertEqual(a.Constraint, b.Constraint)
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')
        self.assertRaises(ConstraintError, b.__setitem__, 1, '3')

    def test_fromkeys(self):
        """ConstrainedDict instance fromkeys should retain constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        b = a.fromkeys('23')
        self.assertEqual(a.Constraint, b.Constraint)
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')
        self.assertRaises(ConstraintError, b.__setitem__, 1, '3')
        b['2'] = 5
        self.assertEqual(b, {'2':5, '3':None})

    def test_setdefault(self):
        """ConstrainedDict setdefault shouldn't allow bad keys"""
        a = ConstrainedDict({'1':None, '2': 'xyz'}, '123')
        self.assertEqual(a.setdefault('2', None), 'xyz')
        self.assertEqual(a.setdefault('1', None), None)
        self.assertRaises(ConstraintError, a.setdefault, 'x', 3)
        a.setdefault('3', 12345)
        self.assertEqual(a, {'1':None, '2':'xyz', '3': 12345})

    def test_update(self):
        """ConstrainedDict should allow update only of compliant data"""
        a = ConstrainedDict(dict.fromkeys('123'), '12345')
        b = ConstrainedDict(dict.fromkeys('444'), '4')
        c = ConstrainedDict(dict.fromkeys('45'), '12345')
        d = ConstrainedDict([['x','y']])
        a.update(b)
        self.assertEqual(a, dict.fromkeys('1234'))
        a.update(c)
        self.assertEqual(a, dict.fromkeys('12345'))
        self.assertRaises(ConstraintError, b.update, c)
        self.assertRaises(ConstraintError, c.update, d)
        #should be OK if constraint removed
        b.Constraint = None
        b.update(c)
        self.assertEqual(b, dict.fromkeys('45'))
        b.update(d)
        self.assertEqual(b, {'4':None, '5':None, 'x':'y'})
        #should fail if we add the constraint back
        b.Constraint = {'4':1, 5:2, '5':1, 'x':1}
        self.assertRaises(ConstraintError, b.update, {4:1})
        b.update({5:1})
        self.assertEqual(b, {'4':None, '5':None, 'x':'y', 5:1})
   
    def test_setitem_masks(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        key_mask = str
        val_mask = lambda x: int(x) + 3
        d = ConstrainedDict({1:4, 2:6}, '123', key_mask, val_mask)
        d[1] = '456'
        self.assertEqual(d, {'1':459,'2':9,})
        d['1'] = 234
        self.assertEqual(d, {'1':237,'2':9,})
        self.assertRaises(ConstraintError, d.__setitem__, 4, '3')
        e = d.copy()
        assert e.Mask is d.Mask
        assert '1' in d
        assert not 1 in d

class MappedDictTests(TestCase):
    """MappedDict should work like ConstrainedDict, but map keys."""
    
    def test_setitem_masks(self):
        """MappedDict setitem should work only if key in constraint"""
        key_mask = str
        val_mask = lambda x: int(x) + 3
        d = MappedDict({1:4, 2:6}, '123', key_mask, val_mask)
        d[1] = '456'
        self.assertEqual(d, {'1':459,'2':9,})
        d['1'] = 234
        self.assertEqual(d, {'1':237,'2':9,})
        self.assertRaises(ConstraintError, d.__setitem__, 4, '3')
        e = d.copy()
        assert e.Mask is d.Mask
        assert '1' in d
        assert 1 in d
        assert 1 not in d.keys()
        assert 'x' not in d.keys()

    def test_getitem(self):
        """MappedDict getitem should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        self.assertEqual(d, {'1':5})
        self.assertEqual(d[1], 5)

    def test_get(self):
        """MappedDict get should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        self.assertEqual(d, {'1':5})
        self.assertEqual(d.get(1, 'x'), 5)
        self.assertEqual(d.get(5, 'x'), 'x')

    def test_has_key(self):
        """MappedDict has_key should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        assert d.has_key('1')
        assert d.has_key(1)
        assert not d.has_key('5')

        


class generateCombinationsTests(TestCase):
    """Tests for public generateCombinations function."""
    
    def test_generateCombinations(self):
        """function should return all combinations of given length"""
        
        #test all 3-position combinations of a 2-digit alphabet, since I can
        #work that one out by hand ...
        correct_result = [  "AAA", "AAB", "ABA", "ABB", \
                            "BBB", "BBA", "BAB", "BAA"]
                            
        real_result = generateCombinations("AB", 3)
        
        correct_result.sort()
        real_result.sort()
        self.assertEquals(str(real_result), str(correct_result))
    #end test_generateCombinations
    
    def test_generateCombinations_singleAlphabet(self):
        """function should return correct value when alphabet is one char"""
        
        real_result = generateCombinations("A", 4)
        self.assertEquals(str(real_result), str(["AAAA"]))
    #end test_generateCombinations_singleAlphabet
    
    def test_generateCombinations_singleLength(self):
        """function should return correct values if length is 1"""
        
        real_result = generateCombinations("ABC", 1)
        self.assertEquals(str(real_result), str(["A", "B", "C"]))
    #end test_generateCombinations_singleLength
    
    def test_generateCombinations_emptyAlphabet(self):
        """function should return empty list if alphabet arg is [], "" """
        
        real_result = generateCombinations("", 4)
        self.assertEquals(str(real_result), str([]))
        
        real_result = generateCombinations([], 4)
        self.assertEquals(str(real_result), str([]))
    #end test_generateCombinations_emptyAlphabet
    
    def test_generateCombinations_zeroLength(self):
        """function should return empty list if length arg is 0 """
        
        real_result = generateCombinations("ABC", 0)
        self.assertEquals(str(real_result), str([]))
    #end test_generateCombinations_zeroLength
    
    def test_generateCombinations_badArgs(self):
        """function should error if args are not castable to right type."""
    
        self.assertRaises(RuntimeError, generateCombinations, 12, 4)
        self.assertRaises(RuntimeError, generateCombinations, [], None) 
    #end test_generateCombinations_badArgs
#end generateCombinationsTests

class makeNonnegIntTests(TestCase):
    """Tests of the public makeNonnegInt function"""
    
    def test_makeNonnegInt_unchanged(self):
        """Should return an input nonneg int unchanged"""
        
        self.assertEquals(makeNonnegInt(3), 3)
    #end test_makeNonnegInt_unchanged
    
    def test_makeNonnegInt_castable(self):
        """Should return nonneg int version of a castable input"""
        
        self.assertEquals(makeNonnegInt(-4.2), 4)
    #end test_makeNonnegInt_castable
    
    def test_makeNonnegInt_noncastable(self):
        """Should raise a special NonnegIntError if input isn't castable"""
        
        self.assertRaises(NonnegIntError, makeNonnegInt, "blue")
    #end test_makeNonnegInt_noncastable
#end makeNonnegIntTests

class revCompTests(TestCase):
    """Tests of the public revComp function"""
    
    def test_revComp_DNA(self):
        """revComp should correctly return reverse complement of DNA"""
        
        #input and correct output taken from example at 
        #http://bioweb.uwlax.edu/GenWeb/Molecular/Seq_Anal/
        #Reverse_Comp/reverse_comp.html
        user_input = "ATGCAGGGGAAACATGATTCAGGAC"
        correct_output = "GTCCTGAATCATGTTTCCCCTGCAT"
        real_output = revComp(user_input)
        self.assertEquals(real_output, correct_output)
    #end test_revComp_DNA
    
    def test_revComp_RNA(self):
        """revComp should correctly return reverse complement of RNA"""
        
        #input and correct output taken from test_revComp_DNA test,
        #with all Ts changed to Us
        user_input = "AUGCAGGGGAAACAUGAUUCAGGAC"
        correct_output = "GUCCUGAAUCAUGUUUCCCCUGCAU"
        
        #remember to use False toggle to get RNA instead of DNA
        real_output = revComp(user_input, False)
        self.assertEquals(real_output, correct_output)        
    #end test_revComp_RNA
    
    def test_revComp_caseSensitive(self):
        """revComp should convert bases without changing case"""
        
        user_input = "aCGtAcgT"
        correct_output = "AcgTaCGt"
        real_output = revComp(user_input)
        self.assertEquals(real_output, correct_output) 
    #end test_revComp_caseSensitive
    
    def test_revComp_nonNucleicSeq(self):
        """revComp should just reverse any chars but ACGT/U"""
        
        user_input = "BDeF"
        correct_output = "FeDB"
        real_output = revComp(user_input)
        self.assertEquals(real_output, correct_output)   
    #end test_revComp_nonNucleicSeq
    
    def test_revComp_emptySeq(self):
        """revComp should return empty string if given empty sequence"""
        
        #shouldn't matter whether in DNA or RNA mode
        real_output = revComp("")
        self.assertEquals(real_output, "") 
    #end test_revComp_emptySeq
    
    def test_revComp_noSeq(self):
        """revComp should return error if given no sequence argument"""
        
        self.assertRaises(TypeError, revComp)
    #end test_revComp_noSeq
#end revCompTests

    def test_not_none(self):
        """not_none should return True if none of the items is None"""
        assert not_none([1,2,3,4])
        assert not not_none([1,2,3,None])
        self.assertEqual(filter(not_none,[(1,2),(3,None)]),[(1,2)])
    #end test_not_none

    def test_get_items_except(self):
        """get_items_except should return all items of seq not in indices"""
        self.assertEqual(get_items_except('a-b-c-d',[1,3,5]),'abcd')
        self.assertEqual(get_items_except([0,1,2,3,4,5,6],[1,3,5]),[0,2,4,6])
        self.assertEqual(get_items_except((0,1,2,3,4,5,6),[1,3,5]),(0,2,4,6))
        self.assertEqual(get_items_except('a-b-c-d',[1,3,5],tuple),
            ('a','b','c','d'))
    #end test_get_items_except    
    
    def test_NestedSplitter(self):
        """NestedSplitter should make a function which return expected list"""
        #test delimiters, constructor, filter_
        line='ii=0; oo= 9, 6 5;  ; xx=  8;  '
        cmds = [
            "NestedSplitter(';=,')(line)",
            "NestedSplitter([';', '=', ','])(line)",
            "NestedSplitter([(';'), '=', ','], constructor=None)(line)",
            "NestedSplitter([(';'), '=', ','], filter_=None)(line)",
            "NestedSplitter([(';',1), '=', ','])(line)",
            "NestedSplitter([(';',-1), '=', ','])(line)"
        ]
        results=[    
            [['ii', '0'], ['oo', ['9', '6 5']], '', ['xx', '8'], ''],
            [['ii', '0'], ['oo', ['9', '6 5']], '', ['xx', '8'], ''],
            [['ii', '0'], [' oo', [' 9', ' 6 5']], '  ', [' xx', '  8'], '  '],
            [['ii', '0'], ['oo', ['9', '6 5']], ['xx', '8']],
            [['ii', '0'], ['oo', ['9', '6 5;  ; xx'], '8;']],
            [['ii', '0; oo', ['9', '6 5;  ; xx'], '8'], '']
        ]
        for cmd, result in zip(cmds, results):
            self.assertEqual(eval(cmd), result)

        #test uncontinous level of delimiters
        test = 'a; b,c; d,e:f; g:h;' #g:h should get [[g,h]] instead of [g,h]
        self.assertEqual(NestedSplitter(';,:')(test),
                ['a', ['b', 'c'], ['d', ['e', 'f']], [['g', 'h']], ''])

        #test empty
        self.assertEqual(NestedSplitter(';,:')(''), [''])
        self.assertEqual(NestedSplitter(';,:')('  '), [''])
        self.assertEqual(NestedSplitter(';,:', filter_=None)(' ;, :'), [[[]]])

        

    def test_curry(self):
        """curry should generate the function with parameters setted"""
        curry_test = curry(cmp, 5)
        knowns = ((3, 1),
                (9, -1),
                (5, 0))
        for arg2, result in knowns:
            self.assertEqual (curry_test(arg2), result)

    def test_app_path(self):
        """app_path should return correct paths"""
        self.assertEqual(app_path('ls'), '/bin/ls')
        self.assertEqual(app_path('lsxxyyx'), False)

#run tests on command-line invocation

if __name__ == '__main__':
    main()

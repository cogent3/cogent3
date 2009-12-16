#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.util.dict2d import Dict2D, \
    average, largest, smallest, swap, nonzero, not_0, upper_to_lower, \
    lower_to_upper, Dict2DInitError, Dict2DError, Dict2DSparseError
from cogent.maths.stats.util import Numbers, Freqs

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Development"

class Dict2DTests(TestCase):
    """ Tests of the Dict2DTests class """

    def setUp(self):
        """Define a few standard matrices"""
        self.empty = {}
        self.single_same = {'a':{'a':2}}
        self.single_diff = {'a':{'b':3}}
        self.square = {
            'a':{'a':1,'b':2,'c':3},
            'b':{'a':2,'b':4,'c':6},
            'c':{'a':3,'b':6,'c':9},
            }
        self.top_triangle = {
            'a':{'a':1, 'b':2, 'c':3},
            'b':{'b':4, 'c':6},
            'c':{'c':9}
        }
        self.bottom_triangle = {
            'b':{'a':2},
            'c':{'a':3, 'b':6}
        }
        self.sparse = {
            'a':{'a':1, 'c':3},
            'd':{'b':2},
        }
        self.dense = {
            'a':{'a':1,'b':2,'c':3},
            'b':{'a':2,'b':4,'c':6},
            }

    def test_init(self):
        """Dict2D init should work as expected"""
        #NOTE: currently only tests init from dict of dicts. Other initializers
        #are tested in the test_guess_input* and test_from* methods
       
        #should compare equal to the relevant dict
        for d in [self.empty, self.single_same, self.single_diff, self.dense, \
            self.sparse]:
            d2d = Dict2D(d)
            self.assertEqual(d2d, d)
            self.assertEqual(d2d.__class__, Dict2D)
        #spot-check values
        d2d = Dict2D(self.sparse)
        self.assertEqual(d2d['a']['c'], 3)
        self.assertEqual(d2d['d']['b'], 2)
        self.assertEqual(len(d2d), 2)
        self.assertEqual(len(d2d['a']), 2)
        self.assertRaises(KeyError, d2d.__getitem__, 'c')
        #check truth values
        assert not Dict2D(self.empty)
        assert Dict2D(self.single_same)

    def test_fromDicts(self):
        """Dict2D.fromDicts should construct from dict of dicts"""
        d2d = Dict2D()
        d2d.fromDicts(self.sparse)
        self.assertEqual(d2d['a']['c'], 3)
        self.assertEqual(d2d['d']['b'], 2)
        self.assertEqual(len(d2d), 2)
        self.assertEqual(len(d2d['a']), 2)
        self.assertRaises(KeyError, d2d.__getitem__, 'c')
        self.assertRaises(Dict2DInitError, d2d.fromDicts, [1,2,3])

    def test_fromIndices(self):
        """Dict2D.fromIndices should construct from list of indices"""
        d2d = Dict2D(self.sparse)
        d2d2 = Dict2D()
        self.assertNotEqual(d2d, d2d2)
        d2d2.fromIndices([('a','a',1),('a','c',3),('d','b',2)])
        self.assertEqual(d2d, d2d2)
        self.assertRaises(Dict2DInitError, d2d2.fromIndices, [1,2,3])

    def test_fromLists(self):
        """Dict2D.fromLists should construct from list of lists"""
        #Note that this only works for dense matrices, not sparse ones
        orig = Dict2D(self.dense)
        new = Dict2D(self.dense)    #will overwrite this data
        self.assertEqual(orig, new)
        assert orig is not new
        new.RowOrder = ['b','a']
        new.ColOrder = ['c','a','b']
        new.fromLists([[3,6,9],[1,3,5]])
        self.assertNotEqual(orig, new)
        test = Dict2D({'b':{'c':3,'a':6,'b':9},'a':{'c':1,'a':3,'b':5}})
        self.assertEqual(new, test)

    def test_guess_input_type_fromLists(self):
        """Dict2D init can correctly guess input type: Lists """
        # Will fail if Error is raised
        d = Dict2D(data=[[1,2,3],[4,5,6]], RowOrder=list('ab'), \
                ColOrder=list('def'))
    
    def test_guess_input_type_fromDict(self):
        """Dict2D init can correctly guess input type: Dict """
        # Will fail if error is raised 
        d = Dict2D({})

    def test_guess_input_type_fromIndices(self):
        """Dict2D init can correctly guess input type: Indices """
        # Will fail if error is raised
        d = Dict2D([('a','b',1)])

    def test_init_without_data(self):
        """Dict2D init functions correctly without a data parameter """
        d = Dict2D(RowOrder=['a'],ColOrder=['b'],Pad=True,Default=42, RowConstructor=Freqs)
        self.assertEqual(d.RowOrder,['a'])
        self.assertEqual(d.ColOrder,['b'])
        self.assertEqual(d.Pad,True)
        self.assertEqual(d.Default,42)
        self.assertEqual(d.RowConstructor, Freqs)
        self.assertEqual(d,{'a':{'b':42.}})

    def test_pad(self):
        """Dict2D pad should fill empty slots with default, but not make square"""
        d = Dict2D(self.sparse)
        d.pad()
        self.assertEqual(len(d), 2)
        self.assertEqual(len(d['a']), 3)
        self.assertEqual(len(d['d']), 3)
        self.assertEqual(d['a'].keys(), d['d'].keys())
        self.assertEqual(d['a']['b'], None)

        #check that it works with a different default value
        d = Dict2D(self.sparse, Default='x')
        d.pad()
        self.assertEqual(d['a']['b'], 'x')

        #check that it works with a different constructor
        d = Dict2D(self.sparse, Default=0, RowConstructor=Freqs)
        d.pad()
        self.assertEqual(d['a']['b'], 0)
        assert isinstance(d['a'], Freqs)

    def test_purge(self):
        """Dict2D purge should delete unwanted keys at both levels"""
        d = Dict2D(self.square)
        d.RowOrder = 'ab'
        d.ColOrder = 'bc'
        d.purge()
        self.assertEqual(d, Dict2D({'a':{'b':2,'c':3},'b':{'b':4,'c':6}}))
        #check that a superset of the keys is OK
        d = Dict2D(self.square)
        d.RowOrder = dict.fromkeys('abcd')
        d.ColOrder = dict.fromkeys('abcd')
        d.purge()
        self.assertEqual(d, Dict2D(self.square))
        #check that everything gets deleted if nothing is valid
        d.RowOrder = list('xyz')
        d.ColOrder = list('xyz')
        d.purge()
        self.assertEqual(d, {})

    def test_rowKeys(self):
        """Dict2D rowKeys should find all the keys of component rows"""
        self.assertEqual(Dict2D(self.empty).rowKeys(), [])
        self.assertEqual(Dict2D(self.single_diff).rowKeys(), ['a'])
        #note that keys will be returned in arbitrary order
        self.assertEqualItems(Dict2D(self.dense).rowKeys(), ['a','b',])
        self.assertEqualItems(Dict2D(self.square).rowKeys(), ['a','b','c'])
        self.assertEqualItems(Dict2D(self.sparse).rowKeys(), ['a','d'])
        
    def test_colKeys(self):
        """Dict2D colKeys should find all the keys of component cols"""
        self.assertEqual(Dict2D(self.empty).colKeys(), [])
        self.assertEqual(Dict2D(self.single_diff).colKeys(), ['b'])
        #note that keys will be returned in arbitrary order
        self.assertEqualItems(Dict2D(self.square).colKeys(), ['a','b','c'])
        self.assertEqualItems(Dict2D(self.dense).colKeys(), ['a','b','c'])
        self.assertEqualItems(Dict2D(self.sparse).colKeys(), ['a','b','c'])

    def test_sharedColKeys(self):
        """Dict2D sharedColKeys should find keys shared by all component cols"""
        self.assertEqual(Dict2D(self.empty).sharedColKeys(), [])
        self.assertEqual(Dict2D(self.single_diff).sharedColKeys(), ['b'])
        #note that keys will be returned in arbitrary order
        self.assertEqualItems(Dict2D(self.square).sharedColKeys(),['a','b','c'])
        self.assertEqualItems(Dict2D(self.dense).sharedColKeys(), ['a','b','c'])
        self.assertEqualItems(Dict2D(self.sparse).sharedColKeys(), [])
        self.square['x'] = {'b':3, 'c':5, 'e':7}
        self.assertEqualItems(Dict2D(self.square).colKeys(),['a','b','c','e'])
        self.assertEqualItems(Dict2D(self.square).sharedColKeys(),['b','c'])

    def test_square(self):
        """Dict2D square should ensure that all d[i][j] exist"""
        #will raise exception if rows and cols aren't equal...
        self.assertRaises(Dict2DError, Dict2D(self.sparse).square)
        self.assertRaises(Dict2DError, Dict2D(self.dense).square)
        #...unless reset_order is True
        d = Dict2D(self.sparse)
        d.square(reset_order=True)
        self.assertEqual(d, Dict2D({
            'a':{'a':1,'b':None,'c':3,'d':None},
            'b':{'a':None, 'b':None, 'c':None, 'd':None},
            'c':{'a':None, 'b':None, 'c':None, 'd':None},
            'd':{'a':None, 'b':2, 'c':None, 'd':None},
            }))
        #Check that passing in a default works too
        d = Dict2D(self.sparse)
        d.square(reset_order=True, default='x')
        self.assertEqual(d, Dict2D({
            'a':{'a':1,'b':'x','c':3,'d':'x'},
            'b':{'a':'x', 'b':'x', 'c':'x', 'd':'x'},
            'c':{'a':'x', 'b':'x', 'c':'x', 'd':'x'},
            'd':{'a':'x', 'b':2, 'c':'x', 'd':'x'},
            }))

    def test_rows(self):
        """Dict2D Rows property should return list in correct order"""
        #should work with no data
        self.assertEqual(list(Dict2D(self.empty).Rows), [])
        #should work on square matrix
        sq = Dict2D(self.square, RowOrder='abc', ColOrder='abc')
        self.assertEqual(list(sq.Rows), [[1,2,3],[2,4,6],[3,6,9]])
        #check that it works when we change the row and col order
        sq.RowOrder = 'ba'
        sq.ColOrder = 'ccb'
        self.assertEqual(list(sq.Rows), [[6,6,4],[3,3,2]])
        #check that it doesn't raise an error on sparse matrices...
        sp = Dict2D(self.sparse)
        rows = list(sp.Rows)
        for r in rows:
            r.sort()
        rows.sort()
        self.assertEqual(rows, [[1,3],[2]])
        #...unless self.RowOrder and self.ColOrder are set...
        sp.RowOrder = 'ad'
        sp.ColOrder = 'abc'
        self.assertRaises(Dict2DSparseError, list, sp.Rows)
        #...and then, only if self.Pad is not set
        sp.Pad = True
        sp.Default = 'xxx'
        self.assertEqual(list(sp.Rows), [[1, 'xxx', 3],['xxx',2,'xxx']])

    def test_cols(self):
        """Dict2D Cols property should return list in correct order"""
        #should work with no data
        self.assertEqual(list(Dict2D(self.empty).Cols), [])
        #should work with square matrix
        sq = Dict2D(self.square, RowOrder='abc', ColOrder='abc')
        self.assertEqual(list(sq.Cols), [[1,2,3],[2,4,6],[3,6,9]])
        #check that it works when we change the row and col order
        sq.RowOrder = 'ba'
        sq.ColOrder = 'ccb'
        self.assertEqual(list(sq.Cols), [[6,3],[6,3],[4,2]])
        #check that it _does_ raise an error on sparse matrices...
        sp = Dict2D(self.sparse)
        self.assertRaises(Dict2DSparseError, list, sp.Cols)
        #...especially if self.RowOrder and self.ColOrder are set...
        sp.RowOrder = 'ad'
        sp.ColOrder = 'abc'
        self.assertRaises(Dict2DSparseError, list, sp.Cols)
        #...and then, only if self.Pad is not set
        sp.Pad = True
        sp.Default = 'xxx'
        self.assertEqual(list(sp.Cols), [[1,'xxx'],['xxx',2],[3,'xxx']])

    def test_items(self):
        """Dict2D Items property should return list in correct order"""
        #should work with no data
        self.assertEqual(list(Dict2D(self.empty).Items), [])
        #should work on square matrix
        sq = Dict2D(self.square, RowOrder='abc', ColOrder='abc')
        self.assertEqual(list(sq.Items), [1,2,3,2,4,6,3,6,9])
        #check that it works when we change the row and col order
        sq.RowOrder = 'ba'
        sq.ColOrder = 'ccb'
        self.assertEqual(list(sq.Items), [6,6,4,3,3,2])
        #check that it doesn't raise an error on sparse matrices...
        sp = Dict2D(self.sparse)
        items = list(sp.Items)
        items.sort()
        self.assertEqual(items, [1,2,3])
        #...unless self.RowOrder and self.ColOrder are set...
        sp.RowOrder = 'ad'
        sp.ColOrder = 'abc'
        self.assertRaises(Dict2DSparseError, list, sp.Items)
        #...and then, only if self.Pad is not set
        sp.Pad = True
        sp.Default = 'xxx'
        self.assertEqual(list(sp.Items), [1, 'xxx', 3,'xxx',2,'xxx'])

    def test_getRows(self):
        """Dict2D getRows should get specified rows"""
        self.assertEqual(Dict2D(self.square).getRows(['a','c']), \
            {'a':{'a':1,'b':2,'c':3},'c':{'a':3,'b':6,'c':9}})
        #should work on sparse matrix
        self.assertEqual(Dict2D(self.sparse).getRows(['d']), {'d':{'b':2}})
        #should raise KeyError if row doesn't exist...
        d = Dict2D(self.sparse)
        self.assertRaises(KeyError, d.getRows, ['c'])
        #...unless we're Padding
        d.Pad = True
        self.assertEqual(d.getRows('c'), {'c':{}})
        #should work when we negate it
        self.assertEqual(Dict2D(self.square).getRows(['a','c'], negate=True),
            {'b':{'a':2,'b':4,'c':6}})

    def test_getRowIndices(self):
        """Dict2D getRowIndices should return indices of rows where f(x) True"""
        d = Dict2D(self.square)
        lt_15 = lambda x: sum(x) < 15
        self.assertEqual(d.getRowIndices(lt_15), ['a','b'])
        #should be bound by RowOrder and ColOrder
        d.RowOrder = d.ColOrder = 'ac'
        self.assertEqual(d.getRowIndices(lt_15), ['a','c'])
        #negate should work
        d.RowOrder = d.ColOrder = None
        self.assertEqual(d.getRowIndices(lt_15, negate=True), ['c'])
        

    def test_getRowsIf(self):
        """Dict2D getRowsIf should return object with rows wher f(x) is True"""
        d = Dict2D(self.square)
        lt_15 = lambda x: sum(x) < 15
        self.assertEqual(d.getRowsIf(lt_15), \
            {'a':{'a':1,'b':2,'c':3},'b':{'a':2,'b':4,'c':6}})
        #should do test by RowOrder, but copy the whole row
        d.RowOrder = d.ColOrder = 'ac'
        self.assertEqual(d.getRowsIf(lt_15), \
            {'a':{'a':1,'b':2,'c':3},'c':{'a':3,'b':6,'c':9}})
        #negate should work
        d.RowOrder = d.ColOrder = None
        self.assertEqual(d.getRowsIf(lt_15, negate=True), \
            {'c':{'a':3,'b':6,'c':9}})

    def test_getCols(self):
        """Dict2D getCols should return object with specified cols only"""
        d = Dict2D(self.square)
        self.assertEqual(d.getCols('bc'), {
            'a':{'b':2, 'c':3},
            'b':{'b':4, 'c':6},
            'c':{'b':6,'c':9},
        })
        #check that it works on ragged matrices
        d = Dict2D(self.top_triangle)
        self.assertEqual(d.getCols('ac'), {
            'a':{'a':1, 'c':3}, 'b':{'c':6}, 'c':{'c':9}
        })
        #check that negate works
        d = Dict2D(self.square)
        self.assertEqual(d.getCols('bc', negate=True), {
            'a':{'a':1}, 'b':{'a':2}, 'c':{'a':3},
            })
    
    def test_getColIndices(self):
        """Dict2D getColIndices should return list of indices of matching cols"""
        d = Dict2D(self.square)
        lt_15 = lambda x: sum(x) < 15
        self.assertEqual(d.getColIndices(lt_15), ['a','b'])
        #check that negate works
        self.assertEqual(d.getColIndices(lt_15, negate=True), ['c'])

    def test_getColsIf(self):
        """Dict2D getColsIf should return new Dict2D with matching cols"""
        d = Dict2D(self.square)
        lt_15 = lambda x: sum(x) < 15
        self.assertEqual(d.getColsIf(lt_15), {
            'a':{'a':1,'b':2},'b':{'a':2,'b':4},'c':{'a':3,'b':6}
        })
        #check that negate works
        self.assertEqual(d.getColsIf(lt_15, negate=True), \
            {'a':{'c':3},'b':{'c':6},'c':{'c':9}})

    def test_getItems(self):
        """Dict2D getItems should return list of relevant items"""
        d = Dict2D(self.square)
        self.assertEqual(d.getItems([('a','a'),('b','c'),('c','a'),('a','a')]),\
            [1,6,3,1])
        #should work on ragged matrices...
        d = Dict2D(self.top_triangle)
        self.assertEqual(d.getItems([('a','c'),('c','c')]), [3,9])
        #...unless absent items are asked for...
        self.assertRaises(KeyError, d.getItems, [('a','a'),('c','a')])
        #...unles self.Pad is True
        d.Pad = True
        self.assertEqual(d.getItems([('a','c'),('c','a')]), [3, None])
        #negate should work -- must specify RowOrder and ColOrder to get
        #results in predictable order
        d.Pad = False
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.getItems([('a','c'),('c','a'),('a','a')], \
            negate=True), [2,4,6,9])

    def test_getItemIndices(self):
        """Dict2D getItemIndices should return indices when f(item) is True"""
        lt_5 = lambda x: x < 5
        d = Dict2D(self.square)
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.getItemIndices(lt_5), \
            [('a','a'),('a','b'),('a','c'),('b','a'),('b','b'),('c','a')])
        self.assertEqual(d.getItemIndices(lt_5, negate=True), \
            [('b','c'),('c','b'),('c','c')])
        d = Dict2D(self.top_triangle)
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.getItemIndices(lt_5), \
            [('a','a'),('a','b'),('a','c'),('b','b')])

    def test_getItemsIf(self):
        """Dict2D getItemsIf should return list of items when f(item) is True"""
        lt_5 = lambda x: x < 5
        d = Dict2D(self.square)
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.getItemsIf(lt_5), [1,2,3,2,4,3])
        self.assertEqual(d.getItemsIf(lt_5, negate=True), [6,6,9])
        d = Dict2D(self.top_triangle)
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.getItemsIf(lt_5), [1,2,3,4])
        self.assertEqual(d.getItemsIf(lt_5, negate=True), [6,9])
        

    def test_toLists(self):
        """Dict2D toLists should convert dict into list of lists"""
        d = Dict2D(self.square)
        d.RowOrder = 'abc'
        d.ColOrder = 'abc'
        self.assertEqual(d.toLists(), [[1,2,3],[2,4,6],[3,6,9]])
        self.assertEqual(d.toLists(headers=True), \
            [['-', 'a', 'b', 'c'],
             ['a', 1, 2, 3],
             ['b', 2, 4, 6],
             ['c', 3, 6, 9],
            ])
        #should raise error if called on sparse matrix...	
        self.assertRaises(Dict2DSparseError, Dict2D(self.sparse).toLists)
        #...unless self.Pad is True
        d = Dict2D(self.sparse)
        d.RowOrder = 'ad'
        d.ColOrder = 'abc'
        d.Pad = True
        d.Default = 'x'
        self.assertEqual(d.toLists(headers=True), \
            [['-','a','b','c'],['a',1,'x',3],['d','x',2,'x']])

        #works without RowOrder or ColOrder
        goal = [[1,2,3],[2,4,6],[3,6,9]]
        
        # headers=False
        d = Dict2D(self.square)
        l = d.toLists()
        for r in l:
            r.sort()
        l.sort()
        self.assertEqual(l,goal)
        
        # headers=True
        d.toLists(headers=True)
        l = d.toLists()
        for r in l:
            r.sort()
        l.sort()
        self.assertEqual(l,goal)

    def test_copy(self):
        """Dict2D copy should copy data and selected attributes"""
        #if it works on sparse matrices, it'll work on dense ones
        s = Dict2D(self.sparse)
        s.Pad = True
        s.RowOrder = 'abc'
        s.ColOrder = 'def'
        s.xxx = 'yyy'   #arbitrary attributes won't be copied
        s2 = s.copy()
        self.assertEqual(s, s2)
        assert s is not s2
        assert not hasattr(s2, 'xxx')
        self.assertEqual(s2.RowOrder, 'abc')
        self.assertEqual(s2.ColOrder, 'def')
        self.assertEqual(s2.Pad, True)
        assert 'Default' not in s2.__dict__
        assert 'RowConstructor' not in s2.__dict__
   
    def test_fill(self):
        """Dict2D fill should fill in specified values"""
        #with no parameters, should just fill in elements that exist
        d = Dict2D(self.sparse)
        d.fill('x')
        self.assertEqual(d, {'a':{'a':'x','c':'x'}, 'd':{'b':'x'}})
        #if cols is set, makes sure all the relevant cols exist in each row
        #doesn't delete extra cols if they are present
        d = Dict2D(self.sparse)
        d.fill('x', cols='bc')
        #note that d[a][a] should not be affected by the fill
        self.assertEqual(d, {'a':{'a':1,'b':'x','c':'x'},\
                             'd':{'b':'x','c':'x'}
        })
        #if rows but not cols is set, should create but not fill rows
        d = Dict2D(self.sparse)
        d.fill('y', rows='ab')
        self.assertEqual(d, {'a':{'a':'y','c':'y'},
                             'b':{},        #new row created
                             'd':{'b':2}  #unaffected since not in rows
                             })
        #if both rows and cols are set, should create and fill rows
        d = Dict2D(self.sparse)
        d.fill('z', rows='abc', cols='abc')
        self.assertEqual(d, {'a':{'a':'z','b':'z','c':'z'},
                             'b':{'a':'z','b':'z','c':'z'},
                             'c':{'a':'z','b':'z','c':'z'},
                             'd':{'b':2}    #unaffected since col skipped
                            })
        #if set_orders is True, should reset RowOrder and ColOrder
        d = Dict2D(self.sparse)
        d.fill('z', rows='abc', cols='xyz', set_orders=True)
        self.assertEqual(d.RowOrder, 'abc')
        self.assertEqual(d.ColOrder, 'xyz')
        d.fill('a', set_orders=True)
        self.assertEqual(d.RowOrder, None)
        self.assertEqual(d.ColOrder, None)

    def test_setDiag(self):
        """Dict2D setDiag should set diagonal to specified value"""
        #should have no effect on empty dict2d
        d = Dict2D(self.empty)
        d.setDiag(0)
        self.assertEqual(d, {})
        #should work on one-element dict
        d = Dict2D(self.single_same)
        d.setDiag(0)
        self.assertEqual(d, {'a':{'a':0}})

        d = Dict2D(self.single_diff)
        d.setDiag(0)
        self.assertEqual(d, {'a':{'a':0,'b':3}})
        #should work on dense dict
        d = Dict2D(self.square)
        d.setDiag(9)
        self.assertEqual(d, {
            'a':{'a':9,'b':2,'c':3},
            'b':{'a':2,'b':9,'c':6},
            'c':{'a':3,'b':6,'c':9},
            })
        #should work on sparse dict, creating cols for rows but not vice versa
        d = Dict2D(self.sparse)
        d.setDiag(-1)
        self.assertEqual(d, {'a':{'a':-1,'c':3},'d':{'b':2,'d':-1}})
    
    def test_scale(self):
        """Dict2D scale should apply f(x) to each d[i][j]"""
        doubler = lambda x: x * 2
        #should have no effect on empty Dict2D
        d = Dict2D(self.empty)
        d.scale(doubler)
        self.assertEqual(d, {})
        #should work on single-element dict
        d = Dict2D(self.single_diff)
        d.scale(doubler)
        self.assertEqual(d, {'a':{'b':6}})
        #should work on dense dict
        d = Dict2D(self.square)
        d.scale(doubler)
        self.assertEqual(d, {
            'a':{'a':2,'b':4,'c':6},
            'b':{'a':4,'b':8,'c':12},
            'c':{'a':6,'b':12,'c':18},
            })
        #should work on sparse dict, not creating any new elements
        d = Dict2D(self.sparse)
        d.scale(doubler)
        self.assertEqual(d, {'a':{'a':2,'c':6},'d':{'b':4}})
   
    def test_transpose(self):
        """Dict2D transpose should work on both dense and sparse matrices, in place"""
        #should do nothing to empty matrix
        d = Dict2D(self.empty)
        d.transpose()
        self.assertEqual(d, {})
        #should do nothing to single-element square matrix
        d = Dict2D(self.single_same)
        d.transpose()
        self.assertEqual(d, {'a':{'a':2}})
        #should reverse single-element non-square matrix
        d = Dict2D(self.single_diff)
        d.transpose()
        self.assertEqual(d, {'b':{'a':3}})
        #should work on sparse matrix
        d = Dict2D(self.sparse)
        d.transpose()
        self.assertEqual(d, {'a':{'a':1}, 'c':{'a':3},'b':{'d':2}})
        #should reverse row and col order
        d = Dict2D(self.dense)
        d.RowOrder = 'ab'
        d.ColOrder = 'abc'
        d.transpose()
        self.assertEqual(d, \
            {'a':{'a':1,'b':2},'b':{'a':2,'b':4},'c':{'a':3,'b':6}})
        self.assertEqual(d.ColOrder, 'ab')
        self.assertEqual(d.RowOrder, 'abc')

    def test_reflect(self):
        """Dict2D reflect should reflect square matrices across diagonal."""
        d = Dict2D(self.top_triangle)
        #should fail if RowOrder and/or ColOrder are unspecified
        self.assertRaises(Dict2DError, d.reflect)
        self.assertRaises(Dict2DError, d.reflect, upper_to_lower)
        d.RowOrder = 'abc'
        self.assertRaises(Dict2DError, d.reflect)
        d.RowOrder = None
        d.ColOrder = 'abc'
        self.assertRaises(Dict2DError, d.reflect)
        #should work if RowOrder and ColOrder are both set
        d.RowOrder = 'abc'
        d.reflect(upper_to_lower)
        self.assertEqual(d, self.square)
        #try it on lower triangle as well -- note that the diagonal won't be
        #set if it's absent.
        d  = Dict2D(self.bottom_triangle)
        d.ColOrder = 'abc'
        d.RowOrder = 'abc'
        d.reflect(lower_to_upper)
        self.assertEqual(d, {
            'a':{'b':2,'c':3},
            'b':{'a':2,'c':6},
            'c':{'a':3,'b':6},
        })
        d = Dict2D({
                'a':{'a':2,'b':4,'c':6},
                'b':{'a':10,'b':20, 'c':30},
                'c':{'a':30, 'b':60, 'c':90},
        })
        d.ColOrder = d.RowOrder = 'abc'
        d.reflect(average)
        self.assertEqual(d, {
            'a':{'a':2,'b':7,'c':18},
            'b':{'a':7,'b':20,'c':45},
            'c':{'a':18,'b':45,'c':90},
        })

    def test_toDelimited(self):
        """Dict2D toDelimited should return delimited string for printing"""
        d = Dict2D(self.square)
        d.RowOrder = d.ColOrder = 'abc'
        self.assertEqual(d.toDelimited(), \
        '-\ta\tb\tc\na\t1\t2\t3\nb\t2\t4\t6\nc\t3\t6\t9')
        self.assertEqual(d.toDelimited(headers=False), \
            '1\t2\t3\n2\t4\t6\n3\t6\t9')
        #set up a custom formatter...
        def my_formatter(x):
            try:
                return '%1.1f' % x
            except:
                return str(x)
        #...and use it
        self.assertEqual(d.toDelimited(headers=True, item_delimiter='x', \
            row_delimiter='y', formatter=my_formatter), \
            '-xaxbxcyax1.0x2.0x3.0ybx2.0x4.0x6.0ycx3.0x6.0x9.0')

        
        
if __name__ == '__main__':
    main()

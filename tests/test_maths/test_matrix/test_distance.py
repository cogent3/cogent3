#!/usr/bin/env python
"""Unit tests for distance matrices.
"""

from cogent.util.unit_test import TestCase, main
from cogent.maths.matrix.distance import DistanceMatrix
from cogent.util.dict2d import largest, Dict2DError, Dict2DSparseError
from cogent.parse.aaindex import AAIndex1Record
from cogent.maths.stats.util import Freqs
from copy import deepcopy

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Production"

class DistanceMatrixTests(TestCase):

    def setUp(self):

        # v : vector
        # m : matrix
        self.default_keys = list('ACDEFGHIKLMNPQRSTVWY')
        # Set up some matrices
        v1 = {'A':1, 'B':2, 'C':3}
        v2 = {'A':4, 'B':5, 'C':6}
        v3 = {'A':7, 'B':8, 'C':9}
        
        self.m1 = {'A':dict(v1),\
                   'B':dict(v2),\
                   'C':dict(v3)}

        v4 = {'A':0, 'B':1, 'C':5}
        v5 = {'A':5, 'B':0, 'C':4, 'X':99}
        v6 = {'A':5, 'B':8, 'C':0}

        self.m2 = {'A':dict(v4),\
                   'B':dict(v5),\
                   'C':dict(v6)}

        self.matrices = [self.m1,self.m2]

        aar_data = dict(zip(self.default_keys, [i*.15 for i in range(20)]))
        # Setup a AAIndex1Record for testing purposes
        self.aar = AAIndex1Record("5", "Some Info",\
                "25", "Greg", "A test",\
                "something", "This is a test, this is only a test",\
                [0.987, 0.783, 1., 0], aar_data)

        # From test_Dict2D, used in tests at end of this file for 
        # inheritance testing
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

    def test_all_init_parameters(self):
        """ All parameters to init are handled correctly """
        # will fail if any paramters are not recognized
        d = DistanceMatrix()
        d = DistanceMatrix(data={})
        d = DistanceMatrix(RowOrder=[])
        d = DistanceMatrix(ColOrder=[])
        d = DistanceMatrix(Pad=True)
        d = DistanceMatrix(Default=42)
        d = DistanceMatrix(data={},RowOrder=[],ColOrder=[],Pad=True,Default=42)

    def test_attribute_init(self):
        """ Proper initialization of all attributes """

        # proper setting to defaults
        d = DistanceMatrix(data={'a':{'a':1}})
        self.assertEqual(d.RowOrder, self.default_keys)
        self.assertEqual(d.ColOrder, self.default_keys)
        self.assertEqual(d.Pad, True)
        self.assertEqual(d.Default, None)
        self.assertEqual(d.RowConstructor, dict)

        # differ from defaults
        d = DistanceMatrix(data={'a':{'b':1}},RowOrder=['a'],\
                ColOrder=['b'],Pad=False,Default=42,RowConstructor=Freqs)
        self.assertEqual(d.RowOrder, ['a'])
        self.assertEqual(d.ColOrder, ['b'])
        self.assertEqual(d.Pad, False)
        self.assertEqual(d.Default, 42)
        self.assertEqual(d.RowConstructor, Freqs)

        # differ from defaults and no data
        d = DistanceMatrix(RowOrder=['a'],\
                ColOrder=['b'],Pad=False,Default=42,RowConstructor=Freqs)
        self.assertEqual(d.RowOrder, ['a'])
        self.assertEqual(d.ColOrder, ['b'])
        self.assertEqual(d.Pad, False)
        self.assertEqual(d.Default, 42)
        self.assertEqual(d.RowConstructor, Freqs)
    
    def test_Order_defaults(self):
        """ RowOrder and ColOrder are set to default as expected """
        for m in self.matrices:
            dm = DistanceMatrix(data=m)
            self.assertEqual(dm.RowOrder, self.default_keys)
            self.assertEqual(dm.ColOrder, self.default_keys)
                               
    def test_Order_parameters(self):
        """ RowOrder and ColOrder are set to paramters as expected """
        row_order = ['a']
        col_order = ['b']
        for m in self.matrices:
            dm = DistanceMatrix(data=m, RowOrder=row_order, ColOrder=col_order)
            self.assertEqual(dm.RowOrder, row_order)
            self.assertEqual(dm.ColOrder, col_order)

    def test_rowKeys(self):
        """ rowKeys functions properly """
        dm = DistanceMatrix(data={'a':{'b':1}})
        goal = self.default_keys + ['a']
        goal.sort()
        actual = dm.rowKeys()
        actual.sort()
        self.assertEqual(actual,goal)

    def test_colKeys(self):
        """ colKeys functions properly """
        dm = DistanceMatrix(data={'a':{'b':1}})
        goal = self.default_keys + ['b']
        goal.sort()
        actual = dm.colKeys()
        actual.sort()
        self.assertEqual(actual,goal)
        
    def test_sharedColKeys(self):
        """ sharedColKeys functions properly """
        # no shared keys b/c a is not in RowOrder and therefore not padded
        dm = DistanceMatrix(data={'a':{'b':1}})
        self.assertEqual(dm.sharedColKeys(),[])

        # shared should be only self.default_keys b/c 'b' not in ColOrder
        dm = DistanceMatrix(data={'a':{'b':1}},\
                RowOrder=self.default_keys + ['a'])
        actual = dm.sharedColKeys()
        actual.sort()
        self.assertEqual(actual, self.default_keys)

        # shared should be self.default_keys + 'b'
        dm = DistanceMatrix(data={'a':{'b':1}},\
                RowOrder=self.default_keys + ['a'],\
                ColOrder=self.default_keys + ['b'])
        actual = dm.sharedColKeys()
        actual.sort()
        self.assertEqual(actual, self.default_keys + ['b'])
    
    def test_default_padding(self):
        """ Default padding functions as expected """
        for m in self.matrices:
            dm = DistanceMatrix(data=m)
            for r in self.default_keys:
                for c in self.default_keys:
                    dm[r][c]
    
    def test_init_data_types(self):
        """ Correct init from varying data types  """
        # No data
        goal = {}.fromkeys(self.default_keys)
        for r in goal:
            goal[r] = {}.fromkeys(self.default_keys)
        dm = DistanceMatrix()
        self.assertEqual(dm,goal)

        # data is dict of dicts
        dm = DistanceMatrix(data={'a':{'b':1}}, Pad=False)
        self.assertEqual(dm,{'a':{'b':1}})

        # data is list of lists
        dm = DistanceMatrix(data=[[1]],RowOrder=['a'],ColOrder=['b'], Pad=False)
        self.assertEqual(dm,{'a':{'b':1}})
        
        # data is in Indices form
        dm = DistanceMatrix(data=[('a','b',1)], Pad=False)
        self.assertEqual(dm,{'a':{'b':1}})

    def test_sparse_init(self):
        """ Init correctly from a sparse dict """
        d = DistanceMatrix(data={'A':{'C':0.}})
        for r in self.default_keys:
            for c in self.default_keys:
                if (r == 'A') and (c == 'C'):
                    self.assertEqual(d[r][c],0.)
                else:
                    self.assertEqual(d[r][c],None)

    def test_dict_integrity(self):
        """ Integrity of key -> value pairs """
        for m in self.matrices:
            dm = DistanceMatrix(data=m)
            self.assertEqual(dm['A']['A'], m['A']['A'])
            self.assertEqual(dm['B']['C'], m['B']['C'])

    def test_attribute_forwarder_integrity(self):
        """ Integrity of attribute forwarding """
        dm = DistanceMatrix(data=self.m2,info=self.aar)
        self.assertEqual(dm.ID, '5')
        self.assertEqual(dm.Correlating, [0.987, 0.783, 1., 0])
        self.assertEqual(dm.Data['C'], 0.15)

    def test_copy(self):
        """ Copy functions as expected"""
        dm = DistanceMatrix(data=self.m2, RowOrder=self.m2.keys(), info=self.aar)
        c = dm.copy()
        self.assertEqual(c['A']['A'],dm['A']['A'])
        self.assertEqual(c.RowOrder,dm.RowOrder)
        self.assertEqual(c.ColOrder,dm.ColOrder)
        self.assertEqual(c.Pad,dm.Pad)
        self.assertEqual(c.Power,dm.Power)

        # Make sure it's a separate object
        c['A']['A'] = 999
        self.assertNotEqual(c['A']['A'],dm['A']['A'])

    def test_attribute_forwarder_integrity_after_copy(self):
        """ Integrity of attribute forwarding following a copy()"""
        dm = DistanceMatrix(data=self.m2, RowOrder=self.m2.keys(), info=self.aar)
        c = dm.copy()
       
        # dm.ID == '5'
        self.assertEqual(c.ID, dm.ID)
        self.assertEqual(c.Correlating, dm.Correlating)
        self.assertEqual(c.Data['R'], dm.Data['R'])

        c.ID = '0'
        self.assertNotEqual(c.ID,dm.ID)
    
    def test_setDiag(self):
        """ setDiag works as expected """
        for m in self.matrices:
            # create a deep copy so we can test against original
            # matrix without it being effected by altering the object 
            # based on it
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=m.keys())
            
            # set diag to 42
            dm.setDiag(42)
            # test that diag is 42
            for k in dm:
                self.assertEqual(dm[k][k],42)

            # test that no diag is unchanged
            self.assertEqual(dm['B']['A'], m['B']['A'])
            self.assertEqual(dm['B']['C'], m['B']['C'])

    def test_scale(self):
        """ Scale correctly applies function to all elements """
        for m in self.matrices:

            # Test square all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=m.keys(), Pad=False)
            dm.scale(lambda x: x**2)
            self.assertEqual(dm['A']['A'],m['A']['A']**2)
            self.assertEqual(dm['B']['A'],m['B']['A']**2)
            self.assertEqual(dm['B']['C'],m['B']['C']**2)

            # Test cube all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=m.keys(), Pad=False)
            dm.scale(lambda x: x**3)
            self.assertEqual(dm['A']['A'],m['A']['A']**3)
            self.assertEqual(dm['B']['A'],m['B']['A']**3)
            self.assertEqual(dm['B']['C'],m['B']['C']**3)

            # Test linearize all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=m.keys(), Pad=False)
            dm.scale(lambda x: 10**-(x/10.0))
            self.assertFloatEqual(dm['A']['A'],10**-(m['A']['A']/10.))
            self.assertFloatEqual(dm['B']['A'],10**-(m['B']['A']/10.))
            self.assertFloatEqual(dm['B']['C'],10**-(m['B']['C']/10.))

    def test_elementPow_valid(self):
        """ elementPow correctly scales all elements and updates self.Power"""
        for m in self.matrices:

            # Test square all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            dm.elementPow(2)
            self.assertEqual(dm.Power, 2)
            self.assertEqual(dm['A']['A'],m['A']['A']**2)
            self.assertEqual(dm['B']['A'],m['B']['A']**2)
            self.assertEqual(dm['B']['C'],m['B']['C']**2)

            # Test cube square root of all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            dm.elementPow(3)
            dm.elementPow(1./2.)
            self.assertEqual(dm.Power, 3./2.)
            self.assertEqual(dm['A']['A'],m['A']['A']**(3./2.))
            self.assertEqual(dm['B']['A'],m['B']['A']**(3./2.))
            self.assertEqual(dm['B']['C'],m['B']['C']**(3./2.))

    def test_elementPow_ignore_invalid(self):
        """ elementPow correctly detects and ignores invalid data"""
        for m in self.matrices:

            # Test square all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            dm['A']['A'] = 'p'
            dm.elementPow(2)
            self.assertEqual(dm.Power, 2.)
            self.assertEqual(dm['A']['A'],'p')

            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            dm['A']['A'] = None
            dm.elementPow(2)
            self.assertEqual(dm.Power, 2.)
            self.assertEqual(dm['A']['A'],None)

    def test_elementPow_error_on_invalid(self):
        """ elementPow correctly raises error on invalid data"""
        for m in self.matrices:

            # Test square all elements
            # explicit tests
            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            dm['A']['A'] = 'p'
            self.assertRaises(TypeError,dm.elementPow,2,ignore_invalid=False)
            
            dm['A']['A'] = None
            self.assertRaises(TypeError,dm.elementPow,2,ignore_invalid=False)
           
    def test_elementPow_invalid_pow(self):
        """ elementPow correctly raises error on invalid power """
        for m in self.matrices:

            n = deepcopy(m)
            dm = DistanceMatrix(data=n, RowOrder=n.keys(),ColOrder=n.keys(),\
                    Pad=False)
            self.assertRaises(TypeError,dm.elementPow,None,ignore_invalid=False)
            self.assertRaises(TypeError,dm.elementPow,'a',ignore_invalid=False)
        
    def test_transpose(self):
        """ transpose functions as expected """
        for m in self.matrices:
            d = DistanceMatrix(data=m)
            t = d.copy()
            t.transpose()
            # Note, this line will fail on a matrix where transpose = original
            self.assertNotEqual(t,d)
            for r in t:
                for c in t[r]:
                    self.assertEqual(t[r][c],d[c][r])
            t.transpose()
            self.assertEqual(t,d)

    def test_reflect(self):
        """ reflect functions as expected """
        for m in self.matrices:
            d = DistanceMatrix(data=m)
            n = d.copy()
            # Only testing one method, all other are tested in superclass, so
            # redundant testing is probably not necessary
            n.reflect(method=largest)
            for r in d.RowOrder:
                for c in d.ColOrder:
                    if d[r][c] > d[c][r]:
                        goal = d[r][c]
                    else:
                        goal = d[c][r]
                    self.assertEqual(n[r][c],goal)
                    self.assertEqual(n[c][r],goal)
        
    ######
    # Following tests copied (and slightly modified) from test_DistanceMatrix and
    # written by Rob Knight. Intended to test inheritance
    #####

    def test_toDelimited(self):
        """DistanceMatrix toDelimited functions as expected"""
        d = DistanceMatrix(self.square,Pad=False)
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


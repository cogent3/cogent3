#!/usr/bin/env python
"""Unit tests for the microarray module, dealing with fake expression data."""
from cogent.util.unit_test import TestCase, main
from cogent.seqsim.microarray import MicroarrayNode
from numpy import ones

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class MicroarrayNodeTests(TestCase):
    """Tests of the MicroarrayNode class"""
    def test_init_empty(self):
        """MicroarrayNode empty init should return new object as expected"""
        m = MicroarrayNode()
        self.assertEqual(m.Length, 0)
        self.assertEqual(m.Array, None)
        self.assertEqual(m.Name, None)
        self.assertEqual(m.Children, [])
        self.assertEqual(m.Parent, None)

    def test_init(self):
        """MicroarrayNode init should return new object w/ correct attributes"""
        m = MicroarrayNode(Name='x')
        self.assertEqual(m.Length, 0)
        self.assertEqual(m.Array, None)
        self.assertEqual(m.Name, 'x')
        self.assertEqual(m.Children, [])
        self.assertEqual(m.Parent, None)
        n = MicroarrayNode(3, 'xyz', Parent=m)
        self.assertEqual(n.Length, 3)
        self.assertEqual(n.Array, 'xyz')
        self.assertEqual(n.Name, None)
        self.assertEqual(n.Children, [])
        assert n.Parent is m

    def test_mutate(self):
        """Microarray mutate should set arrays appropriately"""
        #check that it works as the root
        a = ones(25, 'float64')
        m = MicroarrayNode()
        m.setExpression(a)
        assert m.Array is not a
        self.assertEqual(m.Array, a)
        #check that it works on a single node w/ branchlength set
        m.Length = 1
        m.setExpression(a)
        self.assertNotEqual(m.Array, a)
        assert min(m.Array) > -4
        assert max(m.Array) < 6
        #check that it works for the children
        m.Length = None
        m2, m3, m4 = MicroarrayNode(), MicroarrayNode(), MicroarrayNode()
        m5, m6, m7 = MicroarrayNode(), MicroarrayNode(), MicroarrayNode()
        m8, m9, m10 = MicroarrayNode(), MicroarrayNode(), MicroarrayNode()
        m.Children = [m2,m3, m4]
        m2.Children = [m5]
        m3.Children = [m6,m7,m8]
        m8.Children = [m9,m10]
        m2.Length = 2 # should be ~ 2 sd from 1
        m3.Length = 0 # should test equal to m.Array
        m4.Length = 0.1 # should be ~ 0.1 sd from 1
        m5.Length = 1 # should be ~ 3 sd from 1
        m6.Length = 0.1 # should be in same bounds as m4
        m7.Length = 2 # should be in same bounds as m2
        m8.Length = 1 # should be ~ 1 sd from 1
        m9.Length = 1 # should be in same bounds as m2
        m10.Length = 0 # should test equal to m8
        m.setExpression(a)
        self.assertNotEqual(m.Array, m2.Array)
        self.assertEqual(m.Array, m3.Array)
        self.assertNotEqual(m.Array, m4.Array)
        self.assertNotEqual(m.Array, m5.Array)
        self.assertNotEqual(m.Array, m6.Array)
        self.assertNotEqual(m.Array, m7.Array)
        self.assertNotEqual(m.Array, m8.Array)
        self.assertNotEqual(m.Array, m9.Array)
        self.assertNotEqual(m.Array, m10.Array)
        self.assertNotEqual(m2.Array, m3.Array)
        self.assertNotEqual(m2.Array, m4.Array)
        self.assertNotEqual(m2.Array, m5.Array)
        self.assertNotEqual(m2.Array, m6.Array)
        self.assertNotEqual(m2.Array, m7.Array)
        self.assertNotEqual(m2.Array, m8.Array)
        self.assertNotEqual(m2.Array, m9.Array)
        self.assertNotEqual(m2.Array, m10.Array)
        self.assertNotEqual(m3.Array, m4.Array)
        self.assertNotEqual(m3.Array, m5.Array)
        self.assertNotEqual(m3.Array, m6.Array)
        self.assertNotEqual(m3.Array, m7.Array)
        self.assertNotEqual(m3.Array, m8.Array)
        self.assertNotEqual(m3.Array, m9.Array)
        self.assertNotEqual(m3.Array, m10.Array)
        self.assertNotEqual(m4.Array, m5.Array)
        self.assertNotEqual(m4.Array, m6.Array)
        self.assertNotEqual(m4.Array, m7.Array)
        self.assertNotEqual(m4.Array, m8.Array)
        self.assertNotEqual(m4.Array, m9.Array)
        self.assertNotEqual(m4.Array, m10.Array)
        self.assertNotEqual(m5.Array, m6.Array)
        self.assertNotEqual(m5.Array, m7.Array)
        self.assertNotEqual(m5.Array, m8.Array)
        self.assertNotEqual(m5.Array, m9.Array)
        self.assertNotEqual(m5.Array, m10.Array)
        self.assertNotEqual(m6.Array, m7.Array)
        self.assertNotEqual(m6.Array, m8.Array)
        self.assertNotEqual(m6.Array, m9.Array)
        self.assertNotEqual(m6.Array, m10.Array)
        self.assertNotEqual(m7.Array, m8.Array)
        self.assertNotEqual(m7.Array, m9.Array)
        self.assertNotEqual(m7.Array, m10.Array)
        self.assertNotEqual(m8.Array, m9.Array)
        self.assertEqual(m8.Array, m10.Array)
        self.assertNotEqual(m9.Array, m10.Array)

        #check that amount of change is about right
        #assert 1 > min(m2.Array) > -15
        #assert 1 < max(m2.Array) < 15
        #assert 1 > min(m4.Array) > 0.4
        #assert 1 < max(m4.Array) < 1.6
        # might want stochastic tests here...
        self.assertIsBetween(m2.Array, -11, 13)
        self.assertIsBetween(m4.Array, 0.4, 1.6)
        self.assertIsBetween(m5.Array, -15, 17)
        self.assertIsBetween(m6.Array, 0.4, 1.6)
        self.assertIsBetween(m7.Array, -11, 13)
        self.assertIsBetween(m8.Array, -4, 7)
        self.assertIsBetween(m9.Array, -11, 13)
        
if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""Unit tests for the microarray module, dealing with fake expression data."""
from cogent.util.unit_test import TestCase, main
from cogent.seqsim.microarray import MicroarrayNode
from numpy import ones

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0.1"
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
        assert min(m.Array) > -3
        assert max(m.Array) < 5
        #check that it works for the children
        m.Length = None
        m2, m3, m4 = MicroarrayNode(), MicroarrayNode(), MicroarrayNode()
        m.Children = [m2,m3, m4]
        m2.Length = 2
        m3.Length = 0
        m4.Length = 0.1
        m.setExpression(a)
        self.assertNotEqual(m.Array, m2.Array)
        self.assertEqual(m.Array, m3.Array)
        self.assertNotEqual(m.Array, m4.Array)
        self.assertNotEqual(m2.Array, m3.Array)
        self.assertNotEqual(m2.Array, m4.Array)
        self.assertNotEqual(m3.Array, m4.Array)
        #check that amount of change is about right
        assert 1 > min(m2.Array) > -15
        assert 1 < max(m2.Array) < 15
        assert 1 > min(m4.Array) > 0.4
        assert 1 < max(m4.Array) < 1.6
        


if __name__ == '__main__':
    main()

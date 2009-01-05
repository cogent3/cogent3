#!/usr/bin/env python
"""Unit tests for tree parsers.
"""
from cogent.parse.newick import parse_string , _Tokeniser as DndTokenizer, \
        TreeParseError as RecordError
        
from cogent.core.tree import TreeBuilder

def DndParser(data):
    assert isinstance(data, basestring), data
    constructor = TreeBuilder().createEdge
    return parse_string(data, constructor)
    
from cogent.core.tree import PhyloNode
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

sample = """
(
(
xyz:0.28124,
(
def:0.24498,
mno:0.03627)
:0.17710)
:0.04870,

abc:0.05925,
(
ghi:0.06914,
jkl:0.13776)
:0.09853);
"""

node_data_sample = """
(
(
xyz:0.28124,
(
def:0.24498,
mno:0.03627)
'A':0.17710)
B:0.04870,

abc:0.05925,
(
ghi:0.06914,
jkl:0.13776)
C:0.09853);
"""


empty = '();'
single = '(abc:3);'
double = '(abc:3, def:4);'
onenest = '(abc:3, (def:4, ghi:5):6 );'
nodedata = '(abc:3, (def:4, ghi:5)jkl:6 );'

class DndTokenizerTests(TestCase):
    """Tests of the DndTokenizer factory function."""

    def test_gdata(self):
        """DndTokenizer should work as expected on real data"""
        exp = \
        ['(', '(', 'xyz', ':', '0.28124',',', '(', 'def', ':', '0.24498',\
        ',', 'mno', ':', '0.03627', ')', ':', '0.17710', ')', ':', '0.04870', \
        ',', 'abc', ':', '0.05925', ',', '(', 'ghi', ':', '0.06914', ',', \
        'jkl', ':', '0.13776', ')', ':', '0.09853', ')', ';', None]
        #split it up for debugging on an item-by-item basis
        obs = list(DndTokenizer(sample))
        self.assertEqual(len(obs), len(exp))
        for i, j in zip(obs, exp):
            self.assertEqual(i, j)
        #try it all in one go
        self.assertEqual(list(DndTokenizer(sample)), exp)


class DndParserTests(TestCase):
    """Tests of the DndParser factory function."""
    def test_gempty(self):
        """DndParser should produce an empty PhyloNode on null data"""
        t = DndParser(empty)
        self.assertFalse(t)
        self.assertEqual(str(t), ';')

    def test_gsingle(self):
        """DndParser should produce a single-child PhyloNode on minimal data"""
        t = DndParser(single)
        self.assertEqual(len(t), 1)
        child = t[0]
        self.assertEqual(child.Name, 'abc')
        self.assertEqual(child.Length, 3)
        self.assertEqual(str(t), '(abc:3.0);')

    def test_gdouble(self):
        """DndParser should produce a double-child PhyloNode from data"""
        t = DndParser(double)
        self.assertEqual(len(t), 2)
        self.assertEqual(str(t), '(abc:3.0,def:4.0);')

    def test_gonenest(self):
        """DndParser should work correctly with nested data"""
        t = DndParser(onenest)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  #first child is terminal
        self.assertEqual(len(t[1]), 2)  #second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0):6.0);')
        
    def test_gnodedata(self):
        """DndParser should assign Name to internal nodes correctly"""
        t = DndParser(nodedata)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  #first child is terminal
        self.assertEqual(len(t[1]), 2)  #second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0)jkl:6.0);')
        info_dict = {}
        for node in t.traverse():
            info_dict[node.Name] = node.Length
        self.assertEqual(info_dict['abc'], 3.0)
        self.assertEqual(info_dict['def'], 4.0)
        self.assertEqual(info_dict['ghi'], 5.0)
        self.assertEqual(info_dict['jkl'], 6.0)

    def test_data(self):
        """DndParser should work as expected on real data"""
        t = DndParser(sample)
        self.assertEqual(str(t), '((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);')
        tdata = DndParser(node_data_sample)
        self.assertEqual(str(tdata), "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)C:0.09853);")

    def test_gbad(self):
        """DndParser should fail if parens unbalanced"""
        left = '((abc:3)'
        right = '(abc:3))'
        self.assertRaises(RecordError, DndParser, left)
        self.assertRaises(RecordError, DndParser, right)

class PhyloNodeTests(TestCase):
    """Check that PhyloNode works the way I think"""
    def test_gops(self):
        """Basic PhyloNode operations should work as expected"""
        p = PhyloNode()
        self.assertEqual(str(p), ';')
        p.Name = 'abc'
        self.assertEqual(str(p), 'abc;')
        p.Length = 3
        self.assertEqual(str(p), 'abc:3;')   #don't suppress branch from root
        q = PhyloNode()
        p.append(q)
        self.assertEqual(str(p), '()abc:3;')
        r = PhyloNode()
        q.append(r)
        self.assertEqual(str(p), '(())abc:3;')
        r.Name = 'xyz'
        self.assertEqual(str(p), '((xyz))abc:3;')
        q.Length = 2
        self.assertEqual(str(p), '((xyz):2)abc:3;')

if __name__ == '__main__':
    main()

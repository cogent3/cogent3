#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from cogent.seqsim.tree import RangeNode, balanced_breakpoints, BalancedTree, \
    RandomTree, CombTree, StarTree, LineTree
from cogent.core.usage import DnaPairs
from copy import deepcopy
from operator import mul, or_, add
from numpy import array, average, diag
from numpy.random import random, randint
from cogent.seqsim.usage import Rates

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class treeTests(TestCase):
    """Tests for top-level functions."""
    def test_init(self):
        """Make sure keyword arguments are being passed to baseclass"""
        node = RangeNode(LeafRange=1, Id=2, Name='foo', Length=42)
        self.assertEqual(node.LeafRange, 1)
        self.assertEqual(node.Id, 2)
        self.assertEqual(node.Name, 'foo')
        self.assertEqual(node.Length, 42)

    def test_balanced_breakpoints(self):
        """balanced_breakpoints should produce expected arrays."""
        self.assertRaises(ValueError, balanced_breakpoints, 1)
        self.assertEqual(balanced_breakpoints(2), array([0]))
        self.assertEqual(balanced_breakpoints(4), array([1,0,2]))
        self.assertEqual(balanced_breakpoints(8), \
            array([3,1,5,0,2,4,6]))
        self.assertEqual(balanced_breakpoints(16), \
            array([7,3,11,1,5,9,13,0,2,4,6,8,10,12,14]))
        self.assertEqual(balanced_breakpoints(32), \
            array([15,7,23,3,11,19,27,1,5,9,13,17,21,25,29,\
            0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]))

    def test_BalancedTree(self):
        """BalancedTree should return a balanced tree"""
        b = BalancedTree(4)
        self.assertEqual(len(list(b.traverse())), 4)
        b.assignIds()
        self.assertEqual(str(b), '((0,1)4,(2,3)5)6')

    def test_RandomTree(self):
        """RandomTree should return correct number of nodes
        
        NOTE: all the work is done in breakpoints, which is thoroughly
        tested indepndently. RandomTree just makes permutations.
        """
        d = {}
        for i in range(10):
            r = RandomTree(100)
            self.assertEqual(len(list(r.traverse())), 100)
            r.assignIds()
            #make sure we get different trees each time...
            s = str(r)
            assert s not in d
            d[s] = None

    def test_CombTree(self):
        """CombTree should return correct topology"""
        c = CombTree(4)
        c.assignIds()
        self.assertEqual(str(c), '(0,(1,(2,3)4)5)6')
        c = CombTree(4, deepest_first=False)
        c.assignIds()
        self.assertEqual(str(c), '(((0,1)4,2)5,3)6')

    def test_StarTree(self):
        """StarTree should return correct star topology and # nodes"""
        t = StarTree(5)
        self.assertEqual(len(t.Children), 5)
        for c in t.Children:
            assert c.Parent is t

    def test_LineTree(self):
        """LineTree should return correct number of nodes"""
        t = LineTree(5)
        depth = 1
        curr = t
        while curr.Children:
            self.assertEqual(len(curr.Children), 1)
            depth += 1
            curr = curr.Children[0]
        self.assertEqual(depth, 5)

class RangeTreeTests(TestCase):
    """Tests of the RangeTree class."""
    def setUp(self):
        """Make some standard objects to test."""
        #Notes on sample string:
        #
        #1. trailing zeros are stripped in conversion to/from float, so result
        #   is only exactly the same without them.
        #
        #2. trailing chars (e.g. semicolon) are not recaptured in the output,
        #   so were deleted from original Newick-format string.
        #
        #3. whitespace is stripped, but is handy for formatting, so is stripped
        #   from original string before comparisons.
        self.sample_tree_string = """
    (
    (
    xyz:0.28124,
    (
    def:0.24498,
    mno:0.03627)
    A:0.1771)
    B:0.0487,

    abc:0.05925,
    (
    ghi:0.06914,
    jkl:0.13776)
    C:0.09853)
    """
        self.t = DndParser(self.sample_tree_string, RangeNode)
        self.i = self.t.indexByAttr('Name')
        
        self.sample_string_2 = '((((a,b),c),(d,e)),((f,g),h))'
        self.t2 = DndParser(self.sample_string_2, RangeNode)
        self.i2 = self.t2.indexByAttr('Name')

        self.sample_string_3 = '(((a,b),c),(d,e))'
        self.t3 = DndParser(self.sample_string_3, RangeNode)

    def test_str(self):
        """RangeNode should round-trip Newick string corrrectly."""

        r = RangeNode()
        self.assertEqual(str(r), '()')
        
        #should work for tree with branch lengths set
        t = DndParser(self.sample_tree_string, RangeNode)
        expected = self.sample_tree_string.replace('\n', '')
        expected = expected.replace(' ', '')
        self.assertEqual(str(t), expected)
        #self.assertEqual(t.getNewick(with_distances=True), expected)
        #should also work for tree w/o branch lengths
        t2 = DndParser(self.sample_string_2, RangeNode)
        self.assertEqual(str(t2), self.sample_string_2)

    def test_traverse(self):
        """RangeTree traverse should visit all nodes in correct order"""
        t = self.t
        i = self.i
        #first, check that lengths are correct
        #naked traverse() only does leaves; should be 6.
        self.assertEqual(len(list(t.traverse())), 6)
        #traverse() with self_before should count all nodes.
        self.assertEqual(len(list(t.traverse(self_before=True))), 10)
        #traverse() with self_after should have same count as self_before
        self.assertEqual(len(list(t.traverse(self_after=True))), 10)
        #traverse() with self_before and self_after should visit internal
        #nodes multiple times
        self.assertEqual(len(list(t.traverse(True,True))), 14)
        
        #now, check that items are in correct order
        exp = ['xyz','def','mno','abc','ghi','jkl']
        obs = [i.Name for i in t.traverse()]
        self.assertEqual(obs, exp)

        exp = [None, 'B', 'xyz', 'A', 'def', 'mno', 'abc', 'C', 'ghi', 'jkl']
        obs = [i.Name for i in t.traverse(self_before=True)]
        self.assertEqual(obs, exp)

        exp = ['xyz', 'def', 'mno', 'A', 'B', 'abc', 'ghi', 'jkl', 'C', None]
        obs = [i.Name for i in t.traverse(self_after=True)]
        self.assertEqual(obs, exp)

        exp = [None, 'B', 'xyz', 'A', 'def', 'mno', 'A', 'B', 'abc', 'C', \
            'ghi', 'jkl', 'C', None]
        obs = [i.Name for i in t.traverse(self_before=True, self_after=True)]
        self.assertEqual(obs, exp)

    def test_indexByAttr(self):
        """RangeNode indexByAttr should make index using correct attr"""
        t = self.t
        i = self.i
        #check that we got the right number of elements
        #all elements unique, so should be same as num nodes
        self.assertEqual(len(i), len(list(t.traverse(self_before=True))))
        #check that we got everything
        i_keys = i.keys()
        i_vals = i.values()
        for node in t.traverse(self_before=True):
            assert node.Name in i_keys
            assert node in i_vals
            #can't predict which node will have None as the key
            if node.Name is not None:
                assert i[node.Name] is node
        #check that it works when elements are not unique
        t = self.t3
        for node in t.traverse(self_before=True):
            node.X = 'b'
        for node in t.traverse():
            node.X = 'a'
        result = t.indexByAttr('X', multiple=True)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result['a']), 5)
        self.assertEqual(len(result['b']), 4)
        for n in t.traverse():
            assert n in result['a']
            assert not n in result['b']
        
    def test_indexbyFunc(self):
        """RangeNode indexByFunc should make index from function"""
        t = self.t
        def f(n):
            try:
                return n.Name.isupper()
            except AttributeError:
                return None
            
        i = self.i
        f_i = t.indexByFunc(f)
        self.assertEqual(len(f_i), 3)
        self.assertEqual(f_i[True], [i['B'], i['A'], i['C']])
        self.assertEqual(f_i[False], [i['xyz'], i['def'], i['mno'], \
            i['abc'], i['ghi'], i['jkl']])
        self.assertEqual(f_i[None], [i[None]])

    def test_assignIds(self):
        """RangeNode assignIds should work as expected"""
        t = self.t2
        index = self.i2
        
        t.assignIds()
        #check that ids were set correctly on the leaves
        for i, a in enumerate('abcdefgh'):
            self.assertEqual(index[a].Id, i)
        #check that ranges were set correctly on the leaves
        for i, a in enumerate('abcdefgh'):
            self.assertEqual(index[a].LeafRange, (i, i+1))

        #check that internal ids were set correctly
        obs = [i.Id for i in t.traverse(self_after=True)]
        exp = [0,1,8,2,9,3,4,10,11,5,6,12,7,13,14]
        self.assertEqual(obs, exp)

        #check that internal ranges were set correctly
        obs = [i.LeafRange for i in t.traverse(self_after=True)]
        exp = [(0,1),(1,2),(0,2),(2,3),(0,3),(3,4),(4,5),(3,5),(0,5), \
            (5,6),(6,7),(5,7),(7,8),(5,8),(0,8)]
        self.assertEqual(obs, exp)

    def test_propagateAttr(self):
        """RangeNode propagateAttr should send attr down tree, unless set"""
        t = self.t
        i = self.i
        for n in t.traverse(self_before=True):
            assert not hasattr(n, 'XYZ')
        t.XYZ = 3
        t.propagateAttr('XYZ')
        for n in t.traverse(self_before=True):
            self.assertEqual(n.XYZ, 3)
        #shouldn't overwrite internal nodes by default
        a_children = list(i['A'].traverse(self_before=True))
        i['A'].GHI = 5
        t.GHI = 1
        t.propagateAttr('GHI')
        for n in t.traverse(self_before=True):
            if n in a_children:
                self.assertEqual(n.GHI, 5)
            else:
                self.assertEqual(n.GHI, 1)

        t.GHI = 0
        t.propagateAttr('GHI', overwrite=True)
        for n in t.traverse(self_before=True):
            self.assertEqual(n.GHI, 0)
        
    def test_delAttr(self):
        """RangeNode delAttr should delete attr from self and children"""
        t = self.t2
        for n in t.traverse(self_before=True):
            assert hasattr(n, 'Name')
        t.delAttr('Name')
        for n in t.traverse(self_before=True):
            assert not hasattr(n, 'Name')

    def test_accumulateAttr(self):
        """RangeNode accumulateAttr should accumulate attr in right direction"""
        t = self.t3
        #test towards_leaves (the default)
        f = lambda a, b: b + 1
        for n in t.traverse(self_before=True):
            n.Level = 0
        t.accumulateAttr('Level', f=f)
        levels = [i.Level for i in t.traverse(self_before=True)]
        self.assertEqual(levels, [0,1,2,3,3,2,1,2,2])
        for n in t.traverse(self_before=True):
            n.Level=0
        #test away from leaves
        f = lambda a, b : max(a, b+1)
        for n in t.traverse(self_before=True):
            n.Level=0
        t.accumulateAttr('Level', towards_leaves=False,f=f)
        levels = [i.Level for i in t.traverse(self_before=True)]
        self.assertEqual(levels, [3,2,1,0,0,0,1,0,0])
    
    def test_accumulateChildAttr(self):
        """RangeNode accumulateChildAttr should work as expected"""
        t = self.t2
        i = self.i2
        i['a'].x = 3
        i['b'].x = 4
        i['d'].x = 0
        i['f'].x = 1
        i['g'].x = 1
        i['h'].x = 1

        t.accumulateChildAttr('x', f=mul)
        self.assertEqual([i.__dict__.get('x', None) for i in \
            t.traverse(self_after=True)],
            [3, 4, 12, None, 12, 0, None, 0, 0, 1, 1, 1, 1, 1, 0])
        
        t.accumulateChildAttr('x', f=add)
        self.assertEqual([i.__dict__.get('x', None) for i in \
            t.traverse(self_after=True)],
            [3, 4, 7, None, 7, 0, None, 0, 7, 1, 1, 2, 1, 3, 10])
       
    def test_assignLevelsFromRoot(self):
        """RangeNode assignLevelsFromRoot should match hand-calculated levels"""
        t = self.t3
        t.assignLevelsFromRoot()
        levels = [i.Level for i in t.traverse(self_before=True)]
        self.assertEqual(levels, [0,1,2,3,3,2,1,2,2])
     
    def test_assignLevelsFromLeaves(self):
        """RangeNode assignLevelsFromLeaves should match hand-calculated levels"""
        t = self.t3
        t.assignLevelsFromLeaves()
        levels = [i.Level for i in t.traverse(self_before=True)]
        self.assertEqual(levels, [3,2,1,0,0,0,1,0,0])
        t.assignLevelsFromLeaves(use_min=True)
        levels = [i.Level for i in t.traverse(self_before=True)]
        self.assertEqual(levels, [2,1,1,0,0,0,1,0,0])

    def test_attrToList(self):
        """RangeNode attrToList should return correct list of attr"""
        t = self.t3
        t.assignIds()
        t.assignLevelsFromRoot()
        #make sure nodes are in the order we expect
        self.assertEqual([n.Id for n in t.traverse(self_before=True)],
            [8,6,5,0,1,2,7,3,4])
        #by default, should return list containing all nodes
        obs = t.attrToList('Level')
        self.assertEqual(obs, [3,3,2,2,2,2,1,1,0])
        #should be able to do leaves only if specified
        obs = t.attrToList('Level', leaves_only=True)
        self.assertEqual(obs, [3,3,2,2,2])
        #should be able to specify larger size
        obs=t.attrToList('Level', size=12)
        self.assertEqual(obs, [3,3,2,2,2,2,1,1,0,None,None,None])
        #should be able to set default
        obs=t.attrToList('Level', default='x', size=12)
        self.assertEqual(obs, [3,3,2,2,2,2,1,1,0,'x','x','x'])

    def test_attrFromList(self):
        """RangeNode attrFromList should set values correctly"""
        t = self.t3
        t.assignIds()
        #by default, should set all nodes from array
        t.attrFromList('Level', [3,3,2,2,2,2,1,1,0])
        self.assertEqual([n.Level for n in t.traverse(self_before=True)], \
            [0,1,2,3,3,2,1,2,2])
        #should also work if we choose to set only the leaves (rest should
        #stay at default values)
        t.Level = -1
        t.propagateAttr('Level', overwrite=True)
        t.attrFromList('Level', [3,3,2,2,2,2,1,1,0], leaves_only=True)
        self.assertEqual([n.Level for n in t.traverse(self_before=True)], \
            [-1,-1,-1,3,3,2,-1,2,2])
        
    def test_toBreakpoints(self):
        """RangeNode toBreakpoints should give expected list"""
        t = self.t2
        t.assignIds()
        self.assertEqual(t.toBreakpoints(), [4,2,1,0,3,6,5])

    def test_fromBreakpoints(self):
        """RangeNode fromBreakpoints should have correct topology"""
        breakpoints = [4,2,1,0,3,6,5]
        t = RangeNode.fromBreakpoints(breakpoints)
        #check number of leaves
        self.assertEqual(len(list(t.traverse())), 8)
        self.assertEqual(len(list(t.traverse(self_before=True))), 15)
        #check that leaves were created in right order
        self.assertEqual([i.Id for i in t.traverse()], range(8))
        #check that whole toplogy is right wrt ids...
        nodes = list(t.traverse(self_before=True))
        obs = [i.Id for i in nodes]
        exp = [8, 9, 11, 13, 0, 1, 2, 12, 3, 4, 10, 14, 5, 6, 7]
        self.assertEqual(obs, exp)
        #...and ranges
        obs = [i.LeafRange for i in nodes]
        exp = [(0,8),(0,5),(0,3),(0,2),(0,1),(1,2),(2,3),(3,5),(3,4),(4,5), \
            (5,8),(5,7),(5,6),(6,7),(7,8)]
        self.assertEqual(obs, exp)

    def test_leafLcaDepths(self):
        """RangeNode leafLcaDepths should return expected depths"""
        t = self.t3
        result = t.leafLcaDepths()
        self.assertEqual(result, array([[0,1,2,3,3],
                                        [1,0,2,3,3],
                                        [2,2,0,3,3],
                                        [3,3,3,0,1],
                                        [3,3,3,1,0]]))

    def test_randomNode(self):
        """RandomNode should hit all nodes equally"""
        t = self.t3
        result = {}
        for i in range(100):
            ans = id(t.randomNode())
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 9)
        for node in t.traverse(self_before=True):
            assert id(node) in result

    def test_randomLeaf(self):
        """RandomLeaf should hit all leaf nodes equally"""
        t = self.t3
        result = {}
        for i in range(100):
            ans = id(t.randomLeaf())
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 5)
        for node in t.traverse():
            assert id(node) in result

    def test_randomNodeWithNLeaves(self):
        """RandomNodeWithNLeaves should return node with correct # leaves"""
        t = self.t3
        #check that only the root gets selected with 5 leaves
        result = {}
        for i in range(20):
            ans = id(t.randomNodeWithNLeaves(5))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 1)
        assert id(t) in result
        #check that nothing has 6 or 4 (for this tree) leaves
        self.assertRaises(KeyError, t.randomNodeWithNLeaves, 6)
        self.assertRaises(KeyError, t.randomNodeWithNLeaves, 4)
        #check that it works with fewer than 5 leaves
        #SINGLE LEAF:
        result = {}
        for i in range(40):
            ans = id(t.randomNodeWithNLeaves(1))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 5)
        self.assertEqual(sum(result.values()), 40)
        #TWO LEAVES:
        result = {}
        for i in range(20):
            ans = id(t.randomNodeWithNLeaves(2))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 2)
        self.assertEqual(sum(result.values()), 20)
        #THREE LEAVES:
        result = {}
        for i in range(20):
            ans = id(t.randomNodeWithNLeaves(3))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 1)
        self.assertEqual(sum(result.values()), 20)

    def test_randomNodeAtLevel(self):
        """RangeNode randomNodeAtLevel should return random node at correct level"""
        t = self.t3
        #LEAVES:
        result = {}
        for i in range(40):
            ans = id(t.randomNodeAtLevel(0))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 5)
        self.assertEqual(sum(result.values()), 40)
        
        #BACK ONE LEVEL:
        result = {}
        for i in range(20):
            ans = id(t.randomNodeAtLevel(1))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 2)
        self.assertEqual(sum(result.values()), 20)
         
        #BACK TWO LEVELS:
        result = {}
        for i in range(20):
            ans = id(t.randomNodeAtLevel(2))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 1)
        self.assertEqual(sum(result.values()), 20)
             
        #BACK THREE LEVELS (to root):
        result = {}
        for i in range(20):
            ans = id(t.randomNodeAtLevel(3))
            if ans not in result:
                result[ans] = 0
            result[ans] += 1
        self.assertEqual(len(result), 1)
        self.assertEqual(sum(result.values()), 20)
        self.assertEqual(result.keys()[0], id(t))
         

    def test_outgroupLast(self):
        """RangeNode outgroupLast should reorder nodes to put outgroup last"""
        t = self.t3
        a, b, c, d, e = t.traverse()
        self.assertEqual(t.outgroupLast(c,a,b), (a, b, c))
        self.assertEqual(t.outgroupLast(c,b,a), (b, a, c))
        self.assertEqual(t.outgroupLast(b,d,a), (b, a, d))
        self.assertEqual(t.outgroupLast(c,d,e), (d, e, c))
        self.assertEqual(t.outgroupLast(a,d,e), (d, e, a))
        self.assertEqual(t.outgroupLast(a,d,b), (a, b, d))
        #check that it works if we suppress the cache
        self.assertEqual(t.outgroupLast(c,a,b, False), (a, b, c))
        self.assertEqual(t.outgroupLast(c,b,a, False), (b, a, c))
        self.assertEqual(t.outgroupLast(b,d,a, False), (b, a, d))
        self.assertEqual(t.outgroupLast(c,d,e, False), (d, e, c))
        self.assertEqual(t.outgroupLast(a,d,e, False), (d, e, a))
        self.assertEqual(t.outgroupLast(a,d,b, False), (a, b, d))

    def test_filter(self):
        """RangeNode filter should keep or omit selected nodes."""
        t_orig = self.t2
        
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        to_keep = map(idx.__getitem__, 'abch')
        curr_leaves = list(t.traverse())
        t.filter(to_keep)
        curr_leaves = list(t.traverse())
        for i in to_keep:
            assert i in curr_leaves
        for i in map(idx.__getitem__, 'defg'):
            assert i not in curr_leaves
        #note that it collapses one-child nodes
        self.assertEqual(str(t), '(((a,b),c),h)')

        #test same thing but omitting
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        to_omit = map(idx.__getitem__, 'abch')
        t.filter(to_omit, keep=False)
        curr_leaves = list(t.traverse())
        for i in to_omit:
            assert i not in curr_leaves
        for i in map(idx.__getitem__, 'defg'):
            assert i in curr_leaves
        #note that it collapses one-child nodes
        self.assertEqual(str(t), '((d,e),(f,g))')

        #test that it works with internal nodes
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        to_omit = [idx['a'].Parent.Parent]
        t.filter(to_omit, keep=False)
        self.assertEqual(str(t), '((d,e),((f,g),h))')

        #test that it adds branch lengths
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        for i in t.traverse(self_after=True):
            i.BranchLength = 1
        to_omit = map(idx.__getitem__, 'abdefg')
        t.filter(to_omit, keep=False)
        self.assertEqual(str(t), '(c,h)')

        #test that it got rid of the temporary '_selected' attribute
        for node in t.traverse(self_before=True):
            assert not hasattr(node, '_selected')

        #if nothing valid in to_keep, should return empty tree
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        to_keep = []
        t.filter(to_keep, keep=True)
        curr_leaves = list(t.traverse())
        assert len(curr_leaves), 0

        #if nothing valid in to_keep, should return empty tree
        t = deepcopy(t_orig)
        idx = t.indexByAttr('Name')
        to_keep = list('abcde')     #note: just labels, not nodes
        t.filter(to_keep, keep=True)
        curr_leaves = list(t.traverse())
        assert len(curr_leaves), 0

    def test_addChildren(self):
        """RangeNode add_children should add specified # children to list"""
        t = RangeNode()
        t2 = RangeNode(Parent=t)
        t.addChildren(5)
        self.assertEqual(len(t.Children), 6)
        assert t.Children[0] is t2
        for c in t.Children:
            assert c.Parent is t

class OldPhyloNodeTests(TestCase):
    """Tests of the PhyloNode class -- these are all now methods of RangeNode."""
    def setUp(self):
        """Make a couple of standard trees"""
        self.t1 = DndParser('((a,(b,c)),(d,e))', RangeNode)
        #selt.t1 indices: ((0,(1,2)5)6,(3,4)7)8
    
    def test_makeIdIndex(self):
        """RangeNode makeIdIndex should assign ids to every node"""
        self.t1.makeIdIndex()
        result = self.t1.IdIndex
        nodes = list(self.t1.traverse(self_before=True))
        #check we got an entry for each node
        self.assertEqual(len(result), len(nodes))
        #check the ids are in the result
        for i in nodes:
            assert hasattr(i, 'Id')
            assert i.Id in result
            
    def test_assignQ_single_passed(self):
        """RangeNode assignQ should propagate single Q param down tree"""
        #should work if Q explicitly passed
        t = self.t1
        Q = ['a']
        t.assignQ(Q)
        for node in t.traverse(self_before=True):
            assert node.Q is Q

    def test_assignQ_single_set(self):
        """RangeNode assignQ should propagate single Q if set"""
        t = self.t1
        Q = ['a']
        assert not hasattr(t, 'Q')
        t.Q = Q
        t.assignQ()
        for node in t.traverse(self_before=True):
            assert node.Q is Q

    def test_assignQ_single_overwrite(self):
        """RangeNode assignQ should overwrite root Q if new Q passed"""
        t = self.t1
        Q = ['a']
        Q2 = ['b']
        t.Q = Q
        t.assignQ(Q2)
        for node in t.traverse(self_before=True):
            assert node.Q is Q2
            assert not node.Q is Q

    def test_assignQ_multiple(self):
        """RangeNode assignQ should propagate multiple Qs"""
        t = self.t1
        Q1 = ['a']
        Q2 = ['b']
        Q3 = ['c']
        t.makeIdIndex()
        t.IdIndex[7].Q = Q1
        t.IdIndex[5].Q = Q2
        t.assignQ(Q3) 
        result = [i.Q for i in t.traverse(self_after=True)]
        assert t.Q is Q3
        self.assertEqual(result, [Q3,Q2,Q2,Q2,Q3,Q1,Q1,Q1,Q3])

    def test_assignQ_multiple_overwrite(self):
        """RangeNode assignQ should allow overwrite"""
        t = self.t1
        Q1 = ['a']
        Q2 = ['b']
        Q3 = ['c']
        t.makeIdIndex()
        t.IdIndex[7].Q = Q1
        t.IdIndex[5].Q = Q2
        t.assignQ(Q3, overwrite=True)
        for i in t.traverse(self_after=True):
            assert i.Q is Q3

    def test_assignQ_special(self):
       """RangeNode assignQ should work with special Qs"""
       t = self.t1
       Q1 = 'a'
       Q2 = 'b'
       Q3 = 'c'
       t.makeIdIndex()
       special = {7:Q1, 1:Q2}
       #won't work if no Q at root
       self.assertRaises(ValueError, t.assignQ, special_qs=special)
       t.assignQ(Q3, special_qs=special)
       result = [i.Q for i in t.traverse(self_after=True)]
       self.assertEqual(result, ['c','b','c','c','c','a','a','a','c'])
        

    def test_assignP(self):
        """RangeNode assignP should work when Qs set."""
        t = self.t1
        for i in t.traverse(self_before=True):
            i.Length = random() * 0.5 #range 0 to 0.5
        t.Q = Rates.random(DnaPairs)
        t.assignQ()
        t.assignP()
        t.assignIds()
        for node in t.traverse(self_after=True):
            if node.Parent is not None:
                self.assertFloatEqual(average(1-diag(node.P._data), axis=0), \
                    node.Length)
       
    def test_assignLength(self):
        """RangeNode assignLength should set branch length"""
        t = self.t1
        t.assignLength(0.3)
        for i in t.traverse(self_before=True):
            self.assertEqual(i.Length, 0.3)
                
    def test_evolve(self):
        """RangeNode evolve should work on a starting vector"""
        t = self.t1
        t.Q = Rates.random(DnaPairs)
        t.assignQ()
        t.assignLength(0.1)
        t.assignP()
        start = array([1,0,2,1,0,0,2,1,2,0,1,2,1,0,2,0,0,3,0,2,1,0,3,1,0,2,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3])
        t.evolve(start)
        for i in t.traverse():
            self.assertEqual(len(i.Sequence), len(start))
            self.assertNotEqual(i.Sequence, start)
        #WARNING: Doesn't test base freqs etc. at this point, but those aren't
        #really evolve()'s responsibity (tested as self.P.mutate(seq) once
        #P is set, which we've already demonstrated works.)

    def test_assignPs(self):
        """RangeNode assignPs should assign multiple scaled P matrices"""
        t = self.t1
        for i in t.traverse(self_before=True):
            i.Length = random() * 0.5 #range 0 to 0.5
        t.Q = Rates.random(DnaPairs)
        t.assignQ()
        t.assignPs([1, 0.5, 0.25])
        t.assignIds()
        for node in t.traverse(self_after=True):
            if node.Parent is not None:
                self.assertEqual(len(node.Ps), 3)
                self.assertFloatEqual(average(1-diag(node.Ps[0]._data), axis=0), \
                    node.Length)
                self.assertFloatEqual(average(1-diag(node.Ps[1]._data), axis=0), \
                    0.5*node.Length)
                self.assertFloatEqual(average(1-diag(node.Ps[2]._data), axis=0), \
                    0.25*node.Length)

    def test_evolveSeqs(self):
        """PhlyoNode evolveSeqs should evolve multiple sequences"""
        t = self.t1
        for i in t.traverse(self_before=True):
            i.Length = 0.5
        t.Q = Rates.random(DnaPairs)
        t.assignQ()
        t.assignPs([1, 1, 0.1])
        t.assignIds()
        orig_seqs = [array(i) for i in [randint(0,4,200), randint(0,4,200), \
            randint(0,4,200)]]
        t.evolveSeqs(orig_seqs)
        for node in t.traverse():   #only look at leaves
            if node.Parent is not None:
                self.assertEqual(len(node.Sequences), 3)
                for orig, new in zip(orig_seqs, node.Sequences):
                    self.assertEqual(len(orig), len(new))
                    self.assertNotEqual(orig, new)
                assert sum(orig_seqs[1]!=node.Sequences[1]) > \
                        sum(orig_seqs[2]!=node.Sequences[2])

 


if __name__ == "__main__":
    main()

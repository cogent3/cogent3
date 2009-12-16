#!/usr/bin/env python
"""Tests of classes for dealing with trees and phylogeny.
"""
from copy import copy, deepcopy
from cogent import LoadTree
from cogent.core.tree import TreeNode, PhyloNode, TreeError
from cogent.parse.tree import DndParser
from cogent.maths.stats.test import correlation
from cogent.util.unit_test import TestCase, main
from numpy import array, arange

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


class TreeTests(TestCase):
    """Tests of top-level functions."""

    def test_LoadTree(self):
        """LoadTree should load a tree from a file or a string"""
        #NOTE: This method now sits in cogent.__init__

        t_str = '(a_a:10,(b_b:2,c_c:4):5);'
        #NOTE: Tree quotes these labels because they have underscores in them.
        result_str = "('a_a':10.0,('b_b':2.0,'c_c':4.0):5.0);"
        t = LoadTree(treestring=t_str)
        #t = DndParser(t_str)
        names = [i.Name for i in t.tips()]
        self.assertEqual(names, ['a_a', 'b_b', 'c_c'])
        self.assertEqual(str(t),result_str) 
        self.assertEqual(t.getNewick(with_distances=True), result_str) 
        t_str = '(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);'
        #NOTE: Tree silently converts spaces to underscores (only for output),
        #presumably for Newick compatibility.
        result_str = "(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);"
        t = LoadTree(treestring=t_str, underscore_unmunge=True)
        #t = DndParser(t_str, unescape_name=True)
        names = [i.Name for i in t.tips()]
        self.assertEqual(names, ['a a', 'b b', 'c c'])
        self.assertEqual(str(t),result_str) 
        self.assertEqual(t.getNewick(with_distances=True),result_str) 



def _new_child(old_node, constructor):
    """Returns new_node which has old_node as its parent."""
    new_node = constructor()
    new_node.Parent = old_node
    if old_node is not None:
        if new_node not in old_node.Children:
            old_node.Children.append(new_node)
    return new_node

tree_std = """\
        ((a:1, b:2, c:3)abc:0.1, (d:4, (e:5, f:6)ef:0.2)def:0.3);
"""
tree_std_dist = \
      [[  0. ,   3. ,   4. ,   5.4,   6.6,   7.6],
       [  3. ,   0. ,   5. ,   6.4,   7.6,   8.6],
       [  4. ,   5. ,   0. ,   7.4,   8.6,   9.6],
       [  5.4,   6.4,   7.4,   0. ,   9.2,  10.2],
       [  6.6,   7.6,   8.6,   9.2,   0. ,  11. ],
       [  7.6,   8.6,   9.6,  10.2,  11. ,   0. ]]
tree_std_tips = ['a', 'b', 'c', 'd', 'e', 'f']

tree_one_level = """(a:1, b:2, c:3)abc;"""

tree_two_level = """((a:1, b:2, c:3)abc:0.1, d:0.3)abcd;"""

tree_one_child = """((a:1, b:2, c:3)abc:0.1, (d:0.2)d_:0.3)abcd;"""
tree_one_child_dist = \
      [[ 0. ,  3. ,  4. ,  1.6],
       [ 3. ,  0. ,  5. ,  2.6],
       [ 4. ,  5. ,  0. ,  3.6],
       [ 1.6,  2.6,  3.6,  0. ]]
tree_one_child_tips = ['a', 'b', 'c', 'd']

class TreeNodeTests(TestCase):
    """Tests of the TreeNode class."""

    def setUp(self):
        """Define some standard TreeNode for testing"""
        self.Empty = TreeNode()
        self.Single = TreeNode(Name='a')
        self.Child = TreeNode(Name='b')
        self.OneChild = TreeNode(Name='a', Children=[self.Child])
        self.Multi = TreeNode(Name = 'a', Children='bcd')
        self.Repeated = TreeNode(Name='x', Children='aaa')
        self.BigName = map(TreeNode, '0123456789')
        self.BigParent = TreeNode(Name = 'x', Children = self.BigName)
        self.Comparisons = map(TreeNode, 'aab')
        
        nodes = dict([(x, TreeNode(x)) for x in 'abcdefgh'])
        nodes['a'].append(nodes['b'])
        nodes['b'].append(nodes['c'])
        nodes['c'].append(nodes['d'])
        nodes['c'].append(nodes['e'])
        nodes['c'].append(nodes['f'])
        nodes['f'].append(nodes['g'])
        nodes['a'].append(nodes['h'])
        self.TreeNode = nodes
        self.TreeRoot = nodes['a']

        self.s = '((H,G),(R,M));'
        self.t = DndParser(self.s, TreeNode)
        self.s2 = '(((H,G),R),M);'
        self.t2 = DndParser(self.s2, TreeNode)
        self.s4 = '(((H,G),(O,R)),X);'
        self.t4 = DndParser(self.s4, TreeNode)
   
    def test_init_empty(self):
        """Empty TreeNode should init OK"""
        t = self.Empty
        self.assertEqual(t.Name, None)
        self.assertEqual(t.Parent, None)
        self.assertEqual(len(t), 0)

    def test_init_full(self):
        """TreeNode should init OK with parent, data, and children"""
        t = self.Empty
        u = TreeNode(Parent=t, Name='abc', Children='xyz')
        self.assertEqual(u.Name, 'abc')
        assert u.Parent is t
        assert u in t
        self.assertEqual(u[0].Name, 'x')
        self.assertEqual(u[1].Name, 'y')
        self.assertEqual(u[2].Name, 'z')
        self.assertEqual(len(u), 3)

    def test_str(self):
        """TreeNode str should give Newick-style representation"""
        #note: name suppressed if None
        self.assertEqual(str(self.Empty), ';')
        self.assertEqual(str(self.OneChild), '(b)a;')
        self.assertEqual(str(self.BigParent), '(0,1,2,3,4,5,6,7,8,9)x;')
        self.BigParent[-1].extend('abc')
        self.assertEqual(str(self.BigParent), '(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;')

    def test_getNewick(self):
        """Should return Newick-style representation"""
        self.assertEqual(self.Empty.getNewick(), ';')
        self.assertEqual(self.OneChild.getNewick(), '(b)a;')
        self.assertEqual(self.BigParent.getNewick(), \
                '(0,1,2,3,4,5,6,7,8,9)x;')
        self.BigParent[-1].extend('abc')
        self.assertEqual(self.BigParent.getNewick(), \
                '(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;')

    def test_cmp(self):
        """TreeNode cmp should compare using id"""
        nodes = self.TreeNode
        self.assertEqual(cmp(nodes['a'], nodes['a']), 0)
        self.assertNotEqual(cmp(nodes['b'], nodes['a']), 0)
        self.assertNotEqual(cmp(nodes['a'], nodes['b']), 0)

    def test_compareName(self):
        """Compare names between TreeNodes"""
        nodes = self.TreeNode
        self.assertEqual(nodes['a'].compareName(nodes['a']), 0)
        self.assertEqual(nodes['a'].compareName(nodes['b']), -1)
        self.assertEqual(nodes['b'].compareName(nodes['a']), 1)

    def test_compareByNames(self):
        """Compare names between trees"""
        self.assertTrue(self.t.compareByNames(self.t2))
        self.assertTrue(self.t.compareByNames(self.t))
        self.assertFalse(self.t.compareByNames(self.t4))

    def test_eq(self):
        """TreeNode should compare equal if same id"""
        t, u, v = self.Comparisons
        self.assertEqual(t, t)
        assert t is not u
        self.assertNotEqual(t, u)
        self.assertNotEqual(t, v)
    
        f = TreeNode(1.0)
        g = TreeNode(1)
        self.assertNotEqual(f, g)
        f.Name += 0.1
        self.assertNotEqual(f, g)

        #however, two TreeNodes that have no name should not compare equal
        f = TreeNode()
        g = TreeNode()
        self.assertNotEqual(f,g)

        f = TreeNode(Name='foo')
        g = f.copy()
        self.assertNotEqual(f, g)

    def test_ne(self):
        """TreeNode should compare ne by id or data"""
        t, u, v = self.Comparisons
        self.assertFalse(t != t)
        self.assertTrue(t != u)

        f = TreeNode(Name='foo')
        g = f.copy()
        self.assertTrue(f != g)


    def test_append(self):
        """TreeNode append should add item to end of self"""
        self.OneChild.append(TreeNode('c'))
        self.assertEqual(len(self.OneChild), 2)
        self.assertEqual(self.OneChild[-1].Name, 'c')
        self.OneChild.append(6)
        self.assertEqual(len(self.OneChild), 3)
        self.assertEqual(self.OneChild[-1].Name, 6)
        #check that refs are updated when moved from one tree to another
        empty = TreeNode()
        empty.append(self.OneChild[-1])
        self.assertEqual(len(empty), 1)
        self.assertEqual(empty[0].Name, 6)
        self.assertEqual(empty[0].Parent, empty)
        self.assertEqual(self.OneChild[-1].Name, 'c')

    def test_extend(self):
        """TreeNode extend should add many items to end of self"""
        self.Empty.extend('abcdefgh')
        data = ''.join([i.Name for i in self.Empty])
        self.assertEqual(data, 'abcdefgh')

    def test_insert(self):
        """TreeNode insert should insert item at specified index"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        parent.insert(3, 5)
        self.assertEqual(len(parent), 11)
        self.assertEqual(parent[3].Name, 5)
        self.assertEqual(parent[4].Name, '3')
        parent.insert(-1, 123)
        self.assertEqual(len(parent), 12)
        self.assertEqual(parent[-1].Name, '9')
        self.assertEqual(parent[-2].Name, 123)

    def test_pop(self):
        """TreeNode pop should remove and return child at specified index"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        last = parent.pop()
        assert last is nodes[-1]
        assert last.Parent is None
        self.assertEqual(len(parent), 9)
        assert parent[-1] is nodes[-2]
        first = parent.pop(0)
        assert first is nodes[0]
        assert first.Parent is None
        self.assertEqual(len(parent), 8)
        assert parent[0] is nodes[1]
        second_to_last = parent.pop(-2)
        assert second_to_last is nodes[-3]

    def test_remove(self):
        """TreeNode remove should remove first match by value, not id"""
        nodes = map(TreeNode, 'abc'*3)
        parent = TreeNode(Children=nodes)
        self.assertEqual(len(parent), 9)
        parent.remove('a')
        self.assertEqual(len(parent), 8)
        self.assertEqual(''.join([i.Name for i in parent]), 'bcabcabc')
        new_node = TreeNode('a')
        parent.remove(new_node)
        self.assertEqual(len(parent), 7)
        self.assertEqual(''.join([i.Name for i in parent]), 'bcbcabc')

    def test_getitem(self):
        """TreeNode getitem should return item or slice"""
        r = self.TreeRoot
        n = self.TreeNode
        assert r[0] is n['b']
        items = n['c'][0:1]
        self.assertEqual(len(items), 1)
        assert items[0] is n['d']
        items = n['c'][0:2]
        self.assertEqual(len(items), 2)
        assert items[0] is n['d']
        assert items[1] is n['e']
        items = n['c'][:]
        self.assertEqual(len(items), 3)
        assert items[0] is n['d']
        assert items[-1] is n['f']
    
    def test_slice(self):
        """TreeNode slicing should return list, not TreeNode"""
        nodes = self.TreeNode
        c, d, e, f = nodes['c'],nodes['d'],nodes['e'],nodes['f']
        assert c[:] is not c
        self.assertEqual(c[:], [d,e,f])
        self.assertEqual(c[1:2], [e])
        self.assertEqual(c[0:3:2], [d,f])

    def test_setitem(self):
        """TreeNode setitem should set item or extended slice of nodes"""
        parent, nodes = self.BigParent, self.BigName
        t = TreeNode(1)
        parent[0] = t
        assert parent[0] is t
        assert t.Parent is parent
        assert nodes[0].Parent is None
        
        u = TreeNode(2)
        parent[-2] = u
        assert parent[8] is u
        assert u.Parent is parent
        assert nodes[8].Parent is None
        
        parent[1:6:2] = 'xyz'
        for i in [1,3,5]:
            assert nodes[i].Parent is None
        self.assertEqual(parent[1].Name, 'x')
        self.assertEqual(parent[3].Name, 'y')
        self.assertEqual(parent[5].Name, 'z')
        for i in parent:
            assert i.Parent is parent

    def test_setslice(self):
        """TreeNode setslice should set old-style slice of nodes"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        parent[5:] = []
        self.assertEqual(len(parent), 5)
        for i in range(5, 10):
            assert nodes[i].Parent is None
        parent[1:3] = 'abcd'
        self.assertEqual(len(parent), 7)
        for i in parent:
            assert i.Parent is parent
        data_list = [i.Name for i in parent]
        self.assertEqual(data_list, list('0abcd34'))
        parent[1:3] = parent[2:3]
        data_list = [i.Name for i in parent]
        self.assertEqual(data_list, list('0bcd34'))

    def test_delitem(self):
        """TreeNode __delitem__ should delete item and set parent to None"""
        self.assertEqual(self.Child.Parent, self.OneChild)
        self.assertEqual(len(self.OneChild), 1)
        del self.OneChild[0]
        self.assertEqual(self.OneChild.Parent, None)
        self.assertEqual(len(self.OneChild), 0)

        nodes = self.BigName
        parent = self.BigParent
        self.assertEqual(len(parent), 10)
        for n in nodes:
            assert n.Parent is parent
        del parent[-1]
        self.assertEqual(nodes[-1].Parent, None)
        self.assertEqual(len(parent), 9)
        del parent[1:6:2]
        self.assertEqual(len(parent), 6)
        for i, n in enumerate(nodes):
            if i in [0,2,4,6,7,8]:
                assert n.Parent is parent
            else:
                assert n.Parent is None

    def test_delslice(self):
        """TreeNode __delslice__ should delete items from start to end"""
        parent = self.BigParent
        nodes = self.BigName
        self.assertEqual(len(parent), 10)
        del parent[3:-2]
        self.assertEqual(len(parent), 5)
        for i, n in enumerate(nodes):
            if i in [3,4,5,6,7]:
               assert n.Parent is None
            else:
                assert n.Parent is parent

    def test_iter(self):
        """TreeNode iter should iterate over children"""
        r = self.TreeRoot
        n = self.TreeNode
        items = list(r)
        assert items[0] is n['b']
        assert items[1] is n['h']
        self.assertEqual(len(items), 2)

    def test_len(self):
        """TreeNode len should return number of children"""
        r = self.TreeRoot
        self.assertEqual(len(r), 2)

    def test_copyRecursive(self):
        """TreeNode.copyRecursive() should produce deep copy"""
        t = TreeNode(['t'])
        u = TreeNode(['u'])
        t.append(u)

        c = u.copy()
        assert c is not u
        assert c.Name == u.Name
        assert c.Name is not u.Name
        #note: Name _is_ same object if it's immutable, e.g. a string.
        #deepcopy doesn't copy data for immutable objects.
    
        #need to check that we also copy arbitrary attributes
        t.XYZ = [3]
        c = t.copy()
        assert c is not t
        assert c[0] is not u
        assert c[0].Name is not u.Name
        assert c[0].Name == u.Name
        assert c.XYZ == t.XYZ
        assert c.XYZ is not t.XYZ

        t = self.TreeRoot
        c = t.copy()
        self.assertEqual(str(c), str(t))

    def test_copy(self):
        """TreeNode.copy() should work on deep trees"""
        t = comb_tree(1024) # should break recursion limit on regular copy
        t.Name = 'foo' 
        t.XYZ = [3]
        t2 = t.copy()
        t3 = t.copy()
        t3.Name = 'bar'

        self.assertEqual(len(t.tips()), 1024)
        self.assertEqual(len(t2.tips()), 1024)
        self.assertEqual(len(t3.tips()), 1024)

        self.assertNotSameObj(t, t2)
        self.assertEqual(t.Name, t2.Name)
        self.assertNotEqual(t.Name, t3.Name)
    
        self.assertEqual(t.XYZ, t2.XYZ)
        self.assertNotSameObj(t.XYZ, t2.XYZ)

        self.assertEqual(t.getNewick(), t2.getNewick())

        t_simple = TreeNode(['t'])
        u_simple = TreeNode(['u'])
        t_simple.append(u_simple)

        self.assertEqual(str(t_simple.copy()), str(t_simple.copy()))

    def test_copyTopology(self):
        """TreeNode.copyTopology() should produce deep copy ignoring attrs"""
        t = TreeNode(['t'])
        u = TreeNode(['u'])
        t.append(u)

        c = u.copyTopology()
        assert c is not u
        self.assertEqual(c.Name, u.Name)
        #note: Name _is_ same object if it's immutable, e.g. a string.
        #deepcopy doesn't copy data for immutable objects.
    
        #need to check that we do not also copy arbitrary attributes
        t.XYZ = [3]
        c = t.copyTopology()
        assert c is not t
        assert c[0] is not u
        assert c[0].Name is not u.Name
        assert c[0].Name == u.Name
        assert not hasattr(c, 'XYZ')

        t = self.TreeRoot
        c = t.copy()
        self.assertEqual(str(c), str(t))


    def _test_copy_copy(self):
        """copy.copy should raise TypeError on TreeNode"""
        t = TreeNode('t')
        u = TreeNode('u')
        t.append(u)
        self.assertRaises(TypeError, copy, t)
        self.assertRaises(TypeError, copy, u)

    def test_deepcopy(self):
        """copy.deepcopy should work on TreeNode"""
        t = TreeNode(['t'])
        u = TreeNode(['u'])
        t.append(u)

        c = deepcopy(u)
        assert c is not u
        assert c.Name == u.Name
        assert c.Name is not u.Name
        #note: Name _is_ same object if it's immutable, e.g. a string.
        #deepcopy doesn't copy data for immutable objects.
    
        #need to check that we also copy arbitrary attributes
        t.XYZ = [3]
        c = deepcopy(t)
        assert c is not t
        assert c[0] is not u
        assert c[0].Name is not u.Name
        assert c[0].Name == u.Name
        assert c.XYZ == t.XYZ
        assert c.XYZ is not t.XYZ

        t = self.TreeRoot
        c = deepcopy(t)
        self.assertEqual(str(c), str(t))


    def test_Parent(self):
        """TreeNode Parent should hold correct data and be mutable"""
        #check initial conditions
        self.assertEqual(self.Single.Parent, None)
        #set parent and check parent/child relations
        self.Single.Parent = self.Empty
        assert self.Single.Parent is self.Empty
        self.assertEqual(self.Empty[0], self.Single)
        assert self.Single in self.Empty
        self.assertEqual(len(self.Empty), 1)
        #reset parent and check parent/child relations
        self.Single.Parent = self.OneChild
        assert self.Single.Parent is self.OneChild
        assert self.Single not in self.Empty
        assert self.Single is self.OneChild[-1]

        #following is added to check that we don't screw up when there are
        #nodes with different ids that still compare equal
        for i in self.Repeated:
            assert i.Parent is self.Repeated
        last = self.Repeated[-1]
        last.Parent = self.OneChild
        self.assertEqual(len(self.Repeated),  2)
        for i in self.Repeated:
            assert i.Parent is self.Repeated
        assert last.Parent is self.OneChild

    def test_indexInParent(self):
        """TreeNode indexInParent should hold correct data"""
        first = TreeNode('a')
        second = TreeNode('b')
        third = TreeNode('c')
        fourth = TreeNode('0', Children=[first, second, third])
        self.assertEqual(len(fourth), 3)
        self.assertEqual(first.indexInParent(), 0)
        self.assertEqual(second.indexInParent(), 1)
        self.assertEqual(third.indexInParent(), 2)
        del fourth[0]
        self.assertEqual(second.indexInParent(), 0)
        self.assertEqual(third.indexInParent(), 1)
        self.assertEqual(len(fourth), 2)
        assert first.Parent is None

    def test_isTip(self):
        """TreeNode isTip should return True if node is a tip"""
        tips = 'degh'
        for n in self.TreeNode.values():
            if n.Name in tips:
                self.assertEqual(n.isTip(), True)
            else:
                self.assertEqual(n.isTip(), False)

    def test_isRoot(self):
        """TreeNode isRoot should return True if parent is None"""
        r = 'a'
        for n in self.TreeNode.values():
            if n.Name in r:
                self.assertEqual(n.isRoot(), True)
            else:
                self.assertEqual(n.isRoot(), False)

    def test_traverse(self):
        """TreeNode traverse should iterate over nodes in tree."""
        e = self.Empty
        s = self.Single
        o = self.OneChild
        m = self.Multi
        r = self.TreeRoot

        self.assertEqual([i.Name for i in e.traverse()], [None])
        self.assertEqual([i.Name for i in e.traverse(False, False)], [None])
        self.assertEqual([i.Name for i in e.traverse(True, True)], [None])

        self.assertEqual([i.Name for i in s.traverse()], ['a'])
        self.assertEqual([i.Name for i in s.traverse(True, True)], ['a'])
        self.assertEqual([i.Name for i in s.traverse(True, False)], ['a'])
        self.assertEqual([i.Name for i in s.traverse(False, True)], ['a'])
        self.assertEqual([i.Name for i in s.traverse(False, False)], ['a'])

        self.assertEqual([i.Name for i in o.traverse()], ['a','b'])
        self.assertEqual([i.Name for i in o.traverse(True, True)],['a','b','a'])
        self.assertEqual([i.Name for i in o.traverse(True, False)], ['a', 'b'])
        self.assertEqual([i.Name for i in o.traverse(False, True)], ['b', 'a'])
        self.assertEqual([i.Name for i in o.traverse(False, False)], ['b'])

        self.assertEqual([i.Name for i in m.traverse()], ['a','b','c','d'])
        self.assertEqual([i.Name for i in m.traverse(True, True)],\
            ['a','b','c','d','a'])
        self.assertEqual([i.Name for i in m.traverse(True, False)], \
            ['a', 'b','c','d'])
        self.assertEqual([i.Name for i in m.traverse(False, True)], \
            ['b', 'c', 'd', 'a'])
        self.assertEqual([i.Name for i in m.traverse(False, False)], \
            ['b', 'c', 'd'])

        self.assertEqual([i.Name for i in r.traverse()], \
            ['a','b','c','d', 'e', 'f', 'g', 'h'])
        self.assertEqual([i.Name for i in r.traverse(True, True)],\
            ['a','b','c','d','e','f','g','f','c','b','h','a'])
        self.assertEqual([i.Name for i in r.traverse(True, False)], \
            ['a', 'b','c','d','e','f','g','h'])
        self.assertEqual([i.Name for i in r.traverse(False, True)], \
            ['d','e','g','f','c','b','h','a'])
        self.assertEqual([i.Name for i in r.traverse(False, False)], \
            ['d','e','g','h'])
        self.assertEqual([i.Name for i in r.traverse(True, True, False)],\
            ['b','c','d','e','f','g','f','c','b','h'])
        self.assertEqual([i.Name for i in r.traverse(True, False, False)], \
            ['b','c','d','e','f','g','h'])
        self.assertEqual([i.Name for i in r.traverse(False, True, False)], \
            ['d','e','g','f','c','b','h'])
        self.assertEqual([i.Name for i in r.traverse(False, False, False)], \
            ['d','e','g','h'])

        #this previously failed
        t = DndParser('((a:6,(b:1,c:2):8):12,(d:3,(e:1,f:1):4):10);')
        t0 = t.Children[0]
        list(t0.traverse(self_before=False, self_after=True))
        list(t0.traverse(self_before=True, self_after=True))

    def test_levelorder(self):
        t = DndParser("(((A,B)C,(D,E)F,(G,H)I)J,(K,L)M)N;")
        exp = ['N','J','M','C','F','I','K','L','A','B','D','E','G','H']
        names = [n.Name for n in t.levelorder()]
        self.assertEqual(names,exp)

    def test_ancestors(self):
        """TreeNode ancestors should provide list of ancestors, deepest first"""
        nodes, tree = self.TreeNode, self.TreeRoot
        self.assertEqual(nodes['a'].ancestors(), [])
        self.assertEqual(nodes['b'].ancestors(), [nodes['a']])
        self.assertEqual(nodes['d'].ancestors(), nodes['f'].ancestors())
        self.assertEqual(nodes['g'].ancestors(), \
            [nodes['f'], nodes['c'], nodes['b'], nodes['a']])

    def test_root(self):
        """TreeNode root() should find root of tree"""
        nodes, root = self.TreeNode, self.TreeRoot
        for i in nodes.values():
            assert i.root() is root

    def test_children(self):
        """TreeNode Children should allow getting/setting children"""
        nodes = self.TreeNode
        for n in nodes:
            node = nodes[n]
            self.assertEqual(list(node), node.Children)

        t = TreeNode(Children='abc')
        self.assertEqual(len(t), 3)
        u, v = TreeNode('u'), TreeNode('v')

        #WARNING: If you set Children directly, Parent refs will _not_ update!
        t.Children = [u,v]

        assert t[0] is u
        assert t[1] is v
        self.assertEqual(len(t), 2)

    def test_siblings(self):
        """TreeNode siblings() should return all siblings, not self"""
        self.assertEqual(self.Empty.siblings(), [])
        self.assertEqual(self.Child.siblings(), [])
        self.assertEqual(self.OneChild.siblings(), [])
        
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']

        self.assertEqual(g.siblings(), [])
        self.assertEqual(f.siblings(), [d,e])
        self.assertEqual(e.siblings(), [d,f])
        self.assertEqual(d.siblings(), [e,f])
        self.assertEqual(c.siblings(), [])
        self.assertEqual(b.siblings(), [h])
        self.assertEqual(h.siblings(), [b])
        self.assertEqual(a.siblings(), [])

    def test_tips(self):
        """TreeNode tips should return all terminal descendants"""
        self.assertEqual(self.Empty.tips(), [])
        self.assertEqual(self.Child.tips(), [])
        self.assertEqual(self.OneChild.tips(), [self.Child])
        
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']

        self.assertEqual(g.tips(), [])
        self.assertEqual(f.tips(), [g])
        self.assertEqual(e.tips(), [])
        self.assertEqual(d.tips(), [])
        self.assertEqual(c.tips(), [d,e,g])
        self.assertEqual(b.tips(), [d,e,g])
        self.assertEqual(h.tips(), [])
        self.assertEqual(a.tips(), [d,e,g,h])

    def test_itertips(self):
        """TreeNode itertips should iterate over terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.Name for i in tree.iterTips()], list('degh')),

    def test_nontips(self):
        """TreeNode nontips should return all non-terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.Name for i in tree.nontips()], list('bcf'))

    def test_iterNonTips(self):
        """TreeNode iterNontips should iterate over non-terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.Name for i in tree.iterNontips()], list('bcf'))
        

    def test_tipChildren(self):
        """TreeNode tipChildren should return all terminal children"""
        self.assertEqual(self.Empty.tipChildren(), [])
        self.assertEqual(self.Child.tipChildren(), [])
        self.assertEqual(self.OneChild.tipChildren(), [self.Child])
        
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']

        self.assertEqual(g.tipChildren(), [])
        self.assertEqual(f.tipChildren(), [g])
        self.assertEqual(e.tipChildren(), [])
        self.assertEqual(d.tipChildren(), [])
        self.assertEqual(c.tipChildren(), [d,e])
        self.assertEqual(b.tipChildren(), [])
        self.assertEqual(h.tipChildren(), [])
        self.assertEqual(a.tipChildren(), [h])

    def test_nonTipChildren(self):
        """TreeNode nonTipChildren should return all non-terminal children"""
        self.assertEqual(self.Empty.nonTipChildren(), [])
        self.assertEqual(self.Child.nonTipChildren(), [])
        self.assertEqual(self.OneChild.nonTipChildren(), [])
        
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']

        self.assertEqual(g.nonTipChildren(), [])
        self.assertEqual(f.nonTipChildren(), [])
        self.assertEqual(e.nonTipChildren(), [])
        self.assertEqual(d.nonTipChildren(), [])
        self.assertEqual(c.nonTipChildren(), [f])
        self.assertEqual(b.nonTipChildren(), [c])
        self.assertEqual(h.nonTipChildren(), [])
        self.assertEqual(a.nonTipChildren(), [b])

    def test_childGroups(self):
        """TreeNode childGroups should divide children by grandchild presence"""
        parent = TreeNode(Children='aababbbaaabbbababbb')
        for node in parent:
            if node.Name == 'a':
                node.append('def')
        groups = parent.childGroups()
        self.assertEqual(len(groups), 10)
        exp_group_sizes = [2,1,1,3,3,3,1,1,1,3]
        obs_group_sizes = [len(i) for i in groups]
        self.assertEqual(obs_group_sizes, exp_group_sizes)

        parent = TreeNode(Children='aab')
        for node in parent:
            if node.Name == 'a':
                node.append('def')
        groups = parent.childGroups()
        self.assertEqual(len(groups), 2)
        self.assertEqual([len(i) for i in groups], [2,1])

        parent = TreeNode(Children='aaaaa')
        groups = parent.childGroups()
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), 5)

        parent = TreeNode(Children='aaba')
        for node in parent:
            if node.Name == 'a':
                node.append('def')
        groups = parent.childGroups()
        self.assertEqual(len(groups), 3)
        self.assertEqual([len(i) for i in groups], [2,1,1])
        
    def test_removeNode(self):
        """TreeNode removeNode should delete node by id, not value"""
        parent = self.Repeated
        children = list(self.Repeated)
        self.assertEqual(len(parent), 3)
        self.assertEqual(parent.removeNode(children[1]), True)
        self.assertEqual(len(parent), 2)
        assert children[0].Parent is parent
        assert children[1].Parent is None
        assert children[2].Parent is parent
        self.assertEqual(children[0].compareName(children[1]), 0)
        self.assertEqual(parent.removeNode(children[1]), False)
        self.assertEqual(len(parent), 2)
        self.assertEqual(parent.removeNode(children[0]), True)
        self.assertEqual(len(parent), 1)
    
    def test_lastCommonAncestor(self):
        """TreeNode LastCommonAncestor should provide last common ancestor"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']
       
        self.assertEqual(a.lastCommonAncestor(a), a)
        self.assertEqual(a.lastCommonAncestor(b), a)
        self.assertEqual(a.lastCommonAncestor(g), a)
        self.assertEqual(a.lastCommonAncestor(h), a)

        self.assertEqual(b.lastCommonAncestor(g), b)
        self.assertEqual(b.lastCommonAncestor(d), b)
        self.assertEqual(b.lastCommonAncestor(a), a)
        self.assertEqual(b.lastCommonAncestor(h), a)

        self.assertEqual(d.lastCommonAncestor(f), c)
        self.assertEqual(d.lastCommonAncestor(g), c)
        self.assertEqual(d.lastCommonAncestor(a), a)
        self.assertEqual(d.lastCommonAncestor(h), a)

        self.assertEqual(g.lastCommonAncestor(g), g)
        self.assertEqual(g.lastCommonAncestor(f), f)
        self.assertEqual(g.lastCommonAncestor(e), c)
        self.assertEqual(g.lastCommonAncestor(c), c)
        self.assertEqual(g.lastCommonAncestor(b), b)
        self.assertEqual(g.lastCommonAncestor(a), a)
        self.assertEqual(g.lastCommonAncestor(h), a)

        t = TreeNode('h')
        for i in [a,b,c,d,e,f,g,h]:
            self.assertEqual(i.lastCommonAncestor(t), None)
            self.assertEqual(t.lastCommonAncestor(i), None)

        u = TreeNode('a', Children=[t])

    def test_separation(self):
        """TreeNode separation should return correct number of edges"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']

        self.assertEqual(a.separation(a), 0)
        self.assertEqual(c.separation(c), 0)
        self.assertEqual(a.separation(b), 1)
        self.assertEqual(a.separation(h), 1)
        self.assertEqual(g.separation(h), 5)
        self.assertEqual(f.separation(d), 2)
        self.assertEqual(f.separation(c), 1)
        self.assertEqual(c.separation(f), 1)


    def test_nameUnnamedNodes(self):
        """nameUnnamedNodes assigns an arbitrary value when Name == None"""
        tree, tree_nodes = self.TreeRoot, self.TreeNode
        tree_nodes['b'].Name = 'node2'
        tree_nodes['c'].Name = None
        tree_nodes['f'].Name = None
        tree_nodes['e'].Name = 'node3'
        tree.nameUnnamedNodes()
        self.assertEqual(tree_nodes['c'].Name, 'node1')
        self.assertEqual(tree_nodes['f'].Name, 'node4')

    def test_makeTreeArray(self):
        """makeTreeArray maps nodes to the descendants in them"""
        tree = self.TreeRoot
        result, node_list = tree.makeTreeArray()
        self.assertEqual(result, \
                array([[1,1,1,1], [1,1,1,0], [1,1,1,0],[0,0,1,0]]))
        nodes = [node.Name for node in node_list]
        self.assertEqual(nodes, ['a', 'b', 'c', 'f'])
        #test if works with a dec_list supplied
        dec_list = ['d', 'added', 'e', 'g', 'h']
        result2, node_list = tree.makeTreeArray(dec_list)
        self.assertEqual(result2, \
                array([[1,0,1,1,1], [1,0,1,1,0], [1,0,1,1,0], [0,0,0,1,0]]))
   
    def test_reassignNames(self):
        """reassignNames should rename node names based on dict mapping"""
        t = self.TreeRoot
        mapping = dict([(x, str(i)) for i,x in enumerate('abfg')])
        exp_names = ['0','1','2','3','c','d','e','h']
        t.reassignNames(mapping)
        obs_names = sorted(t.getNodeNames())
        self.assertEqual(obs_names, exp_names)

    def test_reassignNames_specific_nodes(self):
        """reassignNames should rename nodes based on dict mapping"""
        t = self.TreeRoot
        nodes = [self.TreeNode['a'], self.TreeNode['b']]
        mapping = dict([(x, str(i)) for i,x in enumerate('abfg')])
        exp_names = ['0','1','c','d','e','f','g','h']
        t.reassignNames(mapping, nodes)
        obs_names = sorted(t.getNodeNames())
        self.assertEqual(obs_names, exp_names)

    def test_getNodesDict(self):
        """getNodesDict returns a dict keyed by name, value is node"""
        t = self.TreeRoot
        nodes = self.TreeNode
        self.assertEqual(t.getNodesDict(), nodes)

    def test_getNodesDict_nonunique_names(self):
        """getNodesDict raises if non unique names are in tree"""
        t = self.TreeRoot
        t.Children[0].Name = 'same'
        t.Children[0].Children[0].Name = 'same'
        self.assertRaises(TreeError, t.getNodesDict)

    def test_removeDeleted(self):
        """removeDeleted should remove all nodes where is_deleted tests true."""
        tree = DndParser('((a:3,(b:2,(c:1,d:1):1):1):2,(e:3,f:3):2);',
            constructor=TreeNode)
        result_not_deleted = deepcopy(tree)
        tree.removeDeleted(lambda x: x.Name in [])
        self.assertEqual(str(tree),str(result_not_deleted))
        deleted = set(['b','d','e','f'])
        result_tree = DndParser('((a:3,((c:1):1):1):2);',constructor=TreeNode)
        is_deleted = lambda x: x.Name in deleted
        tree.removeDeleted(is_deleted)
        self.assertEqual(str(tree),str(result_tree))
    
    def test_prune(self):
        """prune should reconstruct correct topology of tree."""
        tree = DndParser('((a:3,((c:1):1):1):2);',constructor=TreeNode)
        tree.prune()
        result_tree = DndParser('((a:3,c:1));',constructor=TreeNode)
        self.assertEqual(str(tree),str(result_tree))

        samename_bug = DndParser("((A,B)SAMENAME,((C,D)SAMENAME));")
        samename_bug.prune()
        exp_tree_str = '((A,B)SAMENAME,(C,D)SAMENAME);'
        self.assertEqual(str(samename_bug), exp_tree_str)
        
    def test_getNodeMatchingName(self):
        """TreeNode getNodeMatchingName should return node that matches name"""
        nodes = self.TreeNode
        root = self.TreeRoot
        assert root.getNodeMatchingName('g') is nodes['g']

    def test_subset(self):
        """subset should return set of leaves that descends from node"""
        t = self.t
        self.assertEqual(t.subset(), frozenset('HGRM'))
        c = t.Children[0]
        self.assertEqual(c.subset(), frozenset('HG'))
        leaf = c.Children[1] 
        self.assertEqual(leaf.subset(), frozenset('')) 

    def test_subsets(self):
        """subsets should return all subsets descending from a set"""
        t = self.t 
        self.assertEqual(t.subsets(), frozenset(
            [frozenset('HG'), frozenset('RM')]))

    def test_compareBySubsets(self):
        """compareBySubsets should return the fraction of shared subsets"""
        result = self.t.compareBySubsets(self.t)
        self.assertEqual(result, 0)

        result = self.t2.compareBySubsets(self.t2)
        self.assertEqual(result, 0)

        result = self.t.compareBySubsets(self.t2)
        self.assertEqual(result, 0.5)

        result = self.t.compareBySubsets(self.t4)
        self.assertEqual(result, 1-2./5)

        result = self.t.compareBySubsets(self.t4, exclude_absent_taxa=True)
        self.assertEqual(result, 1-2./3)

        result = self.t.compareBySubsets(self.TreeRoot, exclude_absent_taxa=True)
        self.assertEqual(result, 1)
       
        result = self.t.compareBySubsets(self.TreeRoot)
        self.assertEqual(result, 1)


class PhyloNodeTests(TestCase):
    """Tests of phylogeny-specific methods."""
    def setUp(self):
        """Creates a standard tree"""
        nodes = dict([(x, PhyloNode(x)) for x in 'abcdefgh'])
        nodes['a'].append(nodes['b'])
        nodes['b'].append(nodes['c'])
        nodes['c'].append(nodes['d'])
        nodes['c'].append(nodes['e'])
        nodes['c'].append(nodes['f'])
        nodes['f'].append(nodes['g'])
        nodes['a'].append(nodes['h'])
        self.TreeNode = nodes
        self.TreeRoot = nodes['a']
        nodes['a'].Length = None
        nodes['b'].Length = 0
        nodes['c'].Length = 3
        nodes['d'].Length = 1
        nodes['e'].Length = 4
        nodes['f'].Length = 2
        nodes['g'].Length = 3
        nodes['h'].Length = 2

        self.s = '((H:1,G:1):2,(R:0.5,M:0.7):3);'
        self.t = DndParser(self.s,PhyloNode)
        self.s3 = '(((H:1,G:1,O:1):2,R:3):1,X:4);'
        self.t3 = DndParser(self.s3, PhyloNode)

    def test_init(self):
        """Check PhyloNode constructor"""
        n = PhyloNode('foo', Length=10)
        self.assertEqual(n.Name,'foo')
        self.assertEqual(n.Length, 10)

        n = PhyloNode('bar')
        self.assertEqual(n.Name, 'bar')
        self.assertEqual(n.Length, None)

        n = PhyloNode()
        self.assertEqual(n.Name, None)
        self.assertEqual(n.Length, None)

    def test_totalDescendingBranchLength(self):
        """totalDescendingBranchLength returns total branchlength below self"""
        t = self.TreeRoot
        exp = 15
        obs = t.totalDescendingBranchLength()
        self.assertEqual(obs, exp)

        node_c = self.TreeNode['c']
        exp = 10
        obs = node_c.totalDescendingBranchLength()
        self.assertEqual(obs, exp)

    def test_tipsWithinDistance(self):
        """tipsWithinDistance returns tips that are within distance from self"""
        t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2)L;"
        t = DndParser(t_str, constructor=PhyloNode)
        nodes = t.getNodesDict()
        e_node = nodes['E']

        exp_at_dist_2 = []
        exp_at_dist_3 = ['A','C','D']
        exp_at_dist_4 = ['A','B','C','D','F']

        obs_at_dist_2 = sorted([n.Name for n in e_node.tipsWithinDistance(2)])
        obs_at_dist_3 = sorted([n.Name for n in e_node.tipsWithinDistance(3)])
        obs_at_dist_4 = sorted([n.Name for n in e_node.tipsWithinDistance(4)])

        self.assertEqual(obs_at_dist_2, exp_at_dist_2)
        self.assertEqual(obs_at_dist_3, exp_at_dist_3)
        self.assertEqual(obs_at_dist_4, exp_at_dist_4)

    def test_tipsWithinDistance_nodistances(self):
        """tipsWithinDistance returns tips that are within distance from self"""
        t_str = "(A,B,(C,D)E,(F,((G,H)I)J)K)L;"
        t = DndParser(t_str, constructor=PhyloNode)
        nodes = t.getNodesDict()
        e_node = nodes['E']

        exp = sorted([n.Name for n in t.tips()])
        obs = sorted([n.Name for n in e_node.tipsWithinDistance(0)])
        self.assertEqual(obs, exp)

    def test_distance(self):
        """PhyloNode Distance should report correct distance between nodes"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']
       
        self.assertEqual(a.distance(a), 0)
        self.assertEqual(a.distance(b), 0)
        self.assertEqual(a.distance(c), 3)
        self.assertEqual(a.distance(d), 4)
        self.assertEqual(a.distance(e), 7)
        self.assertEqual(a.distance(f), 5)
        self.assertEqual(a.distance(g), 8)
        self.assertEqual(a.distance(h), 2)
    
        self.assertEqual(b.distance(a), 0)
        self.assertEqual(b.distance(b), 0)
        self.assertEqual(b.distance(c), 3)
        self.assertEqual(b.distance(d), 4)
        self.assertEqual(b.distance(e), 7)
        self.assertEqual(b.distance(f), 5)
        self.assertEqual(b.distance(g), 8)
        self.assertEqual(b.distance(h), 2)
        
        self.assertEqual(c.distance(a), 3)
        self.assertEqual(c.distance(b), 3)
        self.assertEqual(c.distance(c), 0)
        self.assertEqual(c.distance(d), 1)
        self.assertEqual(c.distance(e), 4)
        self.assertEqual(c.distance(f), 2)
        self.assertEqual(c.distance(g), 5)
        self.assertEqual(c.distance(h), 5)
        
        self.assertEqual(d.distance(a), 4)
        self.assertEqual(d.distance(b), 4)
        self.assertEqual(d.distance(c), 1)
        self.assertEqual(d.distance(d), 0)
        self.assertEqual(d.distance(e), 5)
        self.assertEqual(d.distance(f), 3)
        self.assertEqual(d.distance(g), 6)
        self.assertEqual(d.distance(h), 6)
        
        self.assertEqual(e.distance(a), 7)
        self.assertEqual(e.distance(b), 7)
        self.assertEqual(e.distance(c), 4)
        self.assertEqual(e.distance(d), 5)
        self.assertEqual(e.distance(e), 0)
        self.assertEqual(e.distance(f), 6)
        self.assertEqual(e.distance(g), 9)
        self.assertEqual(e.distance(h), 9)
        
        self.assertEqual(f.distance(a), 5)
        self.assertEqual(f.distance(b), 5)
        self.assertEqual(f.distance(c), 2)
        self.assertEqual(f.distance(d), 3)
        self.assertEqual(f.distance(e), 6)
        self.assertEqual(f.distance(f), 0)
        self.assertEqual(f.distance(g), 3)
        self.assertEqual(f.distance(h), 7)
        
        self.assertEqual(g.distance(a), 8)
        self.assertEqual(g.distance(b), 8)
        self.assertEqual(g.distance(c), 5)
        self.assertEqual(g.distance(d), 6)
        self.assertEqual(g.distance(e), 9)
        self.assertEqual(g.distance(f), 3)
        self.assertEqual(g.distance(g), 0)
        self.assertEqual(g.distance(h), 10)

        self.assertEqual(h.distance(a), 2)
        self.assertEqual(h.distance(b), 2)
        self.assertEqual(h.distance(c), 5)
        self.assertEqual(h.distance(d), 6)
        self.assertEqual(h.distance(e), 9)
        self.assertEqual(h.distance(f), 7)
        self.assertEqual(h.distance(g), 10)
        self.assertEqual(h.distance(h), 0)

    def test_compareByTipDistances(self):
        obs = self.t.compareByTipDistances(self.t3)
        #note: common taxa are H, G, R (only)
        m1 = array([[0,2,6.5],[2,0,6.5],[6.5,6.5,0]])
        m2 = array([[0,2,6],[2,0,6],[6,6,0]])
        r = correlation(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1-r)/2)

    def test_compareByTipDistances_sample(self):
        obs = self.t.compareByTipDistances(self.t3, sample=3, shuffle_f=sorted)
        #note: common taxa are H, G, R (only)
        m1 = array([[0,2,6.5],[2,0,6.5],[6.5,6.5,0]])
        m2 = array([[0,2,6],[2,0,6],[6,6,0]])
        r = correlation(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1-r)/2)

        # 4 common taxa, still picking H, G, R
        s = '((H:1,G:1):2,(R:0.5,M:0.7,Q:5):3);'
        t = DndParser(self.s,PhyloNode)
        s3 = '(((H:1,G:1,O:1):2,R:3,Q:10):1,X:4);'
        t3 = DndParser(self.s3, PhyloNode)
        obs = t.compareByTipDistances(t3, sample=3, shuffle_f=sorted)

    def test_tipToTipDistances_endpoints(self):
        """Test getting specifc tip distances  with tipToTipDistances"""
        nodes = [self.t.getNodeMatchingName('H'), 
                 self.t.getNodeMatchingName('G'),
                 self.t.getNodeMatchingName('M')]
        names = ['H','G','M']
        exp = (array([[0,2.0,6.7],[2.0,0,6.7],[6.7,6.7,0.0]]), nodes)
        obs = self.t.tipToTipDistances(endpoints=names)
        self.assertEqual(obs, exp)
        
        obs = self.t.tipToTipDistances(endpoints=nodes)
        self.assertEqual(obs, exp)

    def test_prune(self):
        """prune should reconstruct correct topology and Lengths of tree."""
        tree = DndParser('((a:3,((c:1):1):1):2);',constructor=PhyloNode)
        tree.prune()
        result_tree = DndParser('((a:3.0,c:3.0):2.0);',constructor=PhyloNode)
        self.assertEqual(str(tree),str(result_tree))


    def test_str(self):
        """PhyloNode str should give expected results"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes['a']
        b = nodes['b']
        c = nodes['c']
        d = nodes['d']
        e = nodes['e']
        f = nodes['f']
        g = nodes['g']
        h = nodes['h']
        
        self.assertEqual(str(h), 'h:2;')
        self.assertEqual(str(f), '(g:3)f:2;')
        self.assertEqual(str(a), '(((d:1,e:4,(g:3)f:2)c:3)b:0,h:2)a;')
        #check that None isn't converted any more
        h.Length = None
        c.Length = None   #need to test both leaf and internal node
        self.assertEqual(str(a), '(((d:1,e:4,(g:3)f:2)c)b:0,h)a;')

    def test_maxTipTipDistance(self):
        """maxTipTipDistance returns the max dist between any pair of tips"""
        nodes, tree = self.TreeNode, self.TreeRoot
        max_dist, tip_pair = tree.maxTipTipDistance()
        self.assertEqual(max_dist, 10)
        self.assertEqual(tip_pair, ('h', 'g'))

    def test__find_midpoint_nodes(self):
        """_find_midpoint_nodes should return nodes surrounding the midpoint"""
        nodes, tree = self.TreeNode, self.TreeRoot
        max_dist = 10
        tip_pair = ('g', 'h')
        result = tree._find_midpoint_nodes(max_dist, tip_pair)
        self.assertEqual(result, (nodes['b'], nodes['c']))
        tip_pair = ('h', 'g')
        result = tree._find_midpoint_nodes(max_dist, tip_pair)
        self.assertEqual(result, (nodes['f'], nodes['c']))

    def test_rootAtMidpoint(self):
        """rootAtMidpoint performs midpoint rooting"""
        nodes, tree = self.TreeNode, self.TreeRoot
        #works when the midpoint falls on an existing edge
        tree1 = deepcopy(tree)
        result = tree1.rootAtMidpoint()
        self.assertEqual(str(result), '((d:1,e:4,h:5):0.0,g:5.0);')
        #works when the midpoint falls between two existing edges
        nodes['f'].Length = 1
        nodes['c'].Length = 4
        result = tree.rootAtMidpoint()
        self.assertEqual(str(result), '((d:1,e:4,g:4)c:1.0,h:5.0);')

    def test_rootAtMidpoint2(self):
        """rootAtMidpoint works when midpoint is on both sides of root"""
        #also checks whether it works if the midpoint is adjacent to a tip
        nodes, tree = self.TreeNode, self.TreeRoot
        nodes['h'].Length = 20
        result = tree.rootAtMidpoint()
        self.assertEqual(str(result), '(h:14.0,(d:1,e:4,g:5)c:9.0);')
    
    def test_setTipDistances(self):
        """setTipDistances should correctly set tip distances."""
        tree = DndParser('(((A1:.1,B1:.1):.1,(A2:.1,B2:.1):.1):.3,((A3:.1,B3:.1):.1,(A4:.1,B4:.1):.1):.3);',constructor=PhyloNode)
        
        #expected distances for a post order traversal
        expected_tip_distances = [0,0,0.1,0,0,0.1,0.2,0,0,0.1,0,0,0.1,0.2,0.5]
        #tips should have distance of 0
        tree.setTipDistances()
        for node in tree.tips():
            self.assertEqual(node.TipDistance,0)
        idx = 0
        for node in tree.traverse(self_before=False,self_after=True):
            self.assertEqual(node.TipDistance,expected_tip_distances[idx])
            idx+=1
    
    def test_scaleBranchLengths(self):
        """scaleBranchLengths should correclty scale branch lengths."""
        tree = DndParser('(((A1:.1,B1:.1):.1,(A2:.1,B2:.1):.1):.3,((A3:.1,B3:.1):.1,(A4:.1,B4:.1):.1):.3);',constructor=PhyloNode)
        tree.scaleBranchLengths(max_length=100,ultrametric=True)
        expected_tree = '(((A1:20,B1:20):20,(A2:20,B2:20):20):60,((A3:20,B3:20):20,(A4:20,B4:20):20):60);'
        self.assertEqual(str(tree),expected_tree)


class Test_tip_tip_distances_I(object):
    """Abstract class for testing different implementations of tip_to_tip."""
    
    def setUp(self):
        """Define a few standard trees"""
        constructor = PhyloNode
        self.root_std = DndParser(tree_std, constructor)
        self.root_one_level = DndParser(tree_one_level, constructor)
        self.root_two_level = DndParser(tree_two_level, constructor)
        self.root_one_child = DndParser(tree_one_child, constructor)
    
    def test_one_level(self):
        """tip_to_tip should work for one-level multifurcating tree"""
        matrix, order = self.fun(self.root_one_level)
        self.assertEqual([i.Name for i in order], list('abc'))
        self.assertEqual(matrix, array([[0,3,4],[3,0,5],[4,5,0]]))
    
    def test_two_level(self):
        """tip_to_tip should work for two-level tree"""
        matrix, order = self.fun(self.root_two_level)
        self.assertEqual([i.Name for i in order], list('abcd'))
        self.assertFloatEqual(matrix, \
            array([[0,3,4,1.4],[3,0,5,2.4],[4,5,0,3.4],[1.4,2.4,3.4,0]]))

class Test_tip_tip_distances_array(Test_tip_tip_distances_I, TestCase):
    """Tests for the array implementation of tip_to_tip distances"""
    
    def setUp(self):
        """Specify which method to call."""
        self.fun = lambda x: x.tipToTipDistances()
        super(Test_tip_tip_distances_array, self).setUp()
    
    def test_std(self):
        """tip_to_tip should work for small but complex tree"""
        dist, tips = self.fun(self.root_std)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_std_dist)
        self.assertEqual(tips, tree_std_tips)
    
    def test_one_child(self):
        """tip_to_tip should work for tree with a single child"""
        dist, tips = self.fun(self.root_one_child)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_one_child_dist)
        self.assertEqual(tips, tree_one_child_tips)

# for use with testing iterative copy method
def comb_tree(num_leaves):
    """Returns a comb node_class tree."""
    branch_child = 1

    root = TreeNode()
    curr = root

    for i in range(num_leaves-1):
        curr.Children[:] = [TreeNode(Parent=curr),TreeNode(Parent=curr)]
        curr = curr.Children[branch_child]
    return root

#run if called from command line
if __name__ == '__main__':
    main()

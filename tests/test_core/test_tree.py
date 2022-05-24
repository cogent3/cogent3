#!/usr/bin/env python
"""Tests of classes for dealing with trees and phylogeny.
"""
import json
import os
import pathlib

from copy import copy, deepcopy
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from numpy import array
from numpy.testing import assert_allclose, assert_equal

from cogent3 import load_tree, make_tree, open_
from cogent3.core.tree import PhyloNode, TreeError, TreeNode
from cogent3.maths.stats.test import correlation
from cogent3.parse.tree import DndParser
from cogent3.util.misc import get_object_provenance


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Rob Knight",
    "Catherine Lozupone",
    "Daniel McDonald",
    "Peter Maxwell",
    "Gavin Huttley",
    "Andrew Butterfield",
    "Matthew Wakefield",
    "Justin Kuczynski",
    "Jens Reeder",
    "Jose Carlos Clemente Litran",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


base_path = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(base_path, "data")


class TreeTests(TestCase):
    """Tests of top-level functions."""

    def test_make_tree(self):
        """make_tree should load a tree from a file or a string"""
        # NOTE: This method now sits in cogent3.__init__

        t_str = "(a_a:10,(b_b:2,c_c:4):5);"
        # NOTE: Tree quotes these labels because they have underscores in them.
        result_str = "('a_a':10.0,('b_b':2.0,'c_c':4.0):5.0);"
        t = make_tree(treestring=t_str)
        # t = DndParser(t_str)
        names = [i.name for i in t.tips()]
        self.assertEqual(names, ["a_a", "b_b", "c_c"])
        self.assertEqual(str(t), result_str)
        self.assertEqual(t.get_newick(with_distances=True), result_str)
        t_str = "(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);"
        # NOTE: Tree silently converts spaces to underscores (only for output),
        # presumably for Newick compatibility.
        result_str = "(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);"
        t = make_tree(treestring=t_str, underscore_unmunge=True)
        # t = DndParser(t_str, unescape_name=True)
        names = [i.name for i in t.tips()]
        self.assertEqual(names, ["a a", "b b", "c c"])
        self.assertEqual(str(t), result_str)
        self.assertEqual(t.get_newick(with_distances=True), result_str)


def _new_child(old_node, constructor):
    """Returns new_node which has old_node as its parent."""
    new_node = constructor()
    new_node.parent = old_node
    if old_node is not None:
        if new_node not in old_node.children:
            old_node.children.append(new_node)
    return new_node


tree_std = """\
        ((a:1, b:2, c:3)abc:0.1, (d:4, (e:5, f:6)ef:0.2)def:0.3);
"""
tree_std_dist = [
    [0.0, 3.0, 4.0, 5.4, 6.6, 7.6],
    [3.0, 0.0, 5.0, 6.4, 7.6, 8.6],
    [4.0, 5.0, 0.0, 7.4, 8.6, 9.6],
    [5.4, 6.4, 7.4, 0.0, 9.2, 10.2],
    [6.6, 7.6, 8.6, 9.2, 0.0, 11.0],
    [7.6, 8.6, 9.6, 10.2, 11.0, 0.0],
]
tree_std_tips = ["a", "b", "c", "d", "e", "f"]

tree_one_level = """(a:1, b:2, c:3)abc;"""

tree_two_level = """((a:1, b:2, c:3)abc:0.1, d:0.3)abcd;"""

tree_one_child = """((a:1, b:2, c:3)abc:0.1, (d:0.2)d_:0.3)abcd;"""
tree_one_child_dist = [
    [0.0, 3.0, 4.0, 1.6],
    [3.0, 0.0, 5.0, 2.6],
    [4.0, 5.0, 0.0, 3.6],
    [1.6, 2.6, 3.6, 0.0],
]
tree_one_child_tips = ["a", "b", "c", "d"]


class TreeNodeTests(TestCase):
    """Tests of the TreeNode class."""

    def setUp(self):
        """Define some standard TreeNode for testing"""
        self.Empty = TreeNode()
        self.Single = TreeNode(name="a")
        self.Child = TreeNode(name="b")
        self.OneChild = TreeNode(name="a", children=[self.Child])
        self.Multi = TreeNode(name="a", children="bcd")
        self.Repeated = TreeNode(name="x", children="aaa")
        self.BigName = list(map(TreeNode, "0123456789"))
        self.BigParent = TreeNode(name="x", children=self.BigName)
        self.Comparisons = list(map(TreeNode, "aab"))

        nodes = dict([(x, TreeNode(x)) for x in "abcdefgh"])
        nodes["a"].append(nodes["b"])
        nodes["b"].append(nodes["c"])
        nodes["c"].append(nodes["d"])
        nodes["c"].append(nodes["e"])
        nodes["c"].append(nodes["f"])
        nodes["f"].append(nodes["g"])
        nodes["a"].append(nodes["h"])
        self.TreeNode = nodes
        self.TreeRoot = nodes["a"]

        self.s = "((H,G),(R,M));"
        self.t = DndParser(self.s, TreeNode)
        self.s2 = "(((H,G),R),M);"
        self.t2 = DndParser(self.s2, TreeNode)
        self.s4 = "(((H,G),(O,R)),X);"
        self.t4 = DndParser(self.s4, TreeNode)

    def test_init_empty(self):
        """Empty TreeNode should init OK"""
        t = self.Empty
        self.assertEqual(t.name, None)
        self.assertEqual(t.parent, None)
        self.assertEqual(len(t), 0)

    def test_init_full(self):
        """TreeNode should init OK with parent, data, and children"""
        t = self.Empty
        u = TreeNode(parent=t, name="abc", children="xyz")
        self.assertEqual(u.name, "abc")
        assert u.parent is t
        assert u in t
        self.assertEqual(u[0].name, "x")
        self.assertEqual(u[1].name, "y")
        self.assertEqual(u[2].name, "z")
        self.assertEqual(len(u), 3)

    def test_str(self):
        """TreeNode str should give Newick-style representation"""
        # note: name suppressed if None
        self.assertEqual(str(self.Empty), ";")
        self.assertEqual(str(self.OneChild), "(b)a;")
        self.assertEqual(str(self.BigParent), "(0,1,2,3,4,5,6,7,8,9)x;")
        self.BigParent[-1].extend("abc")
        self.assertEqual(str(self.BigParent), "(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;")

    def test_get_newick(self):
        """Should return Newick-style representation"""
        self.assertEqual(self.Empty.get_newick(), ";")
        self.assertEqual(self.OneChild.get_newick(), "(b)a;")
        self.assertEqual(self.BigParent.get_newick(), "(0,1,2,3,4,5,6,7,8,9)x;")
        self.BigParent[-1].extend("abc")
        self.assertEqual(self.BigParent.get_newick(), "(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;")

    def test_to_dict(self):
        """tree produces dict"""
        tr = make_tree(treestring="(a,b,(c,d)e1)")
        got = tr.to_rich_dict()
        attrs = {"length": None}
        expect = {
            "newick": "(a,b,(c,d)e1)",
            "edge_attributes": {
                "a": attrs,
                "b": attrs,
                "c": attrs,
                "d": attrs,
                "e1": attrs,
                "root": attrs,
            },
            "type": get_object_provenance(tr),
            "version": __version__,
        }
        self.assertEqual(got, expect)

        tr = make_tree(treestring="(a:1,b:1,(c:1,d:1)e1:1)")
        got = tr.to_rich_dict()
        attrs = {"length": 1.0}
        expect = {
            "newick": "(a,b,(c,d)e1)",
            "edge_attributes": {
                "a": attrs,
                "b": attrs,
                "c": attrs,
                "d": attrs,
                "e1": attrs,
                "root": {"length": None},
            },
            "type": get_object_provenance(tr),
            "version": __version__,
        }
        self.assertEqual(got, expect)

    def test_to_json(self):
        """tree produces json string that round trips correctly"""
        tr = make_tree(treestring="(a,b,(c,d)e1)")
        got = json.loads(tr.to_json())
        attrs = {"length": None}
        expect = tr.to_rich_dict()
        self.assertEqual(got, expect)

        tr = make_tree(treestring="(a:1,b:1,(c:1,d:1)e1:1)")
        got = json.loads(tr.to_json())
        attrs = {"length": 1.0}
        expect = tr.to_rich_dict()
        self.assertEqual(got, expect)

    def test_write_to_json(self):
        tree = load_tree(filename=os.path.join(data_path, "brca1_5.tree"))
        with TemporaryDirectory(dir=".") as dirname:
            json_path = os.path.join(dirname, "brca1_5.json")
            tree.write(json_path)
            with open_(json_path) as fn:
                got = json.loads(fn.read())
                self.assertEqual(got["type"], get_object_provenance(PhyloNode))
                self.assertEqual(
                    tree.get_newick(semicolon=False, with_node_names=True),
                    got["newick"],
                )
                self.assertEqual(
                    set(tree.get_node_names()), got["edge_attributes"].keys()
                )

    def test_write_to_txt(self):
        """write a tree to newick"""
        tree = load_tree(filename=os.path.join(data_path, "brca1_5.tree"))
        with TemporaryDirectory(dir=".") as dirname:
            out_path = os.path.join(dirname, "brca1_5.txt")
            tree.write(out_path)
            with open_(out_path) as fn:
                got = fn.read()
                self.assertTrue(got.count("(") == got.count(")") == 3)

    def test_write_to_xml(self):
        """write a tree to xml"""
        tree = load_tree(filename=os.path.join(data_path, "brca1_5.tree"))
        with TemporaryDirectory(dir=".") as dirname:
            out_path = os.path.join(dirname, "brca1_5.xml")
            tree.write(out_path)
            with open_(out_path) as fn:
                got = fn.read()
                self.assertTrue(got.count("<clade>") == got.count("</clade>") > 0)

    def test_multifurcating(self):
        """Coerces nodes to have <= n children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)

        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.multifurcating(2)
        self.assertEqual(obs.get_newick(with_distances=True), exp_str)
        self.assertNotEqual(
            t.get_newick(with_distances=True), obs.get_newick(with_distances=True)
        )

        obs = t.multifurcating(2, 0.5)
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.5)d:4.0,((e:5.0,(f:6.0,g:7.0):0.5)h:8.0,(i:9.0,(j:10.0,k:11.0):0.5)l:12.0):0.5)m:14.0;"
        self.assertEqual(obs.get_newick(with_distances=True), exp_str)

        t_str = "((a,b,c)d,(e,f,g)h,(i,j,k)l)m;"
        exp_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
        t = DndParser(t_str, constructor=TreeNode)
        obs = t.multifurcating(2)
        self.assertEqual(obs.get_newick(with_distances=True), exp_str)
        obs = t.multifurcating(2, eps=10)  # no effect on TreeNode type
        self.assertEqual(obs.get_newick(with_distances=True), exp_str)

        self.assertRaises(TreeError, t.multifurcating, 1)

    def test_multifurcating_nameunnamed(self):
        """Coerces nodes to have <= n children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)

        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.multifurcating(2, name_unnamed=True)

        c0, c1 = obs.children
        self.assertTrue(c0.children[1].name.startswith("AUTO"))
        self.assertTrue(c1.name.startswith("AUTO"))
        self.assertTrue(c1.children[0].children[1].name.startswith("AUTO"))
        self.assertTrue(c1.children[1].children[1].name.startswith("AUTO"))
        self.assertEqual(len(c0.children[1].name), 22)
        self.assertEqual(len(c1.name), 22)
        self.assertEqual(len(c1.children[0].children[1].name), 22)
        self.assertEqual(len(c1.children[1].children[1].name), 22)
        names = [n.name for n in t.nontips()]
        self.assertEqual(len(names), len(set(names)))

    def test_bifurcating(self):
        """Coerces nodes to have <= 2 children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)

        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        t.bifurcating()

    def test_eq(self):
        """TreeNode comparison should compare using id"""
        nodes = self.TreeNode
        self.assertEqual(nodes["a"] == nodes["a"])
        self.assertNotEqual(nodes["b"] == nodes["a"])
        self.assertNotEqual(nodes["a"], nodes["b"])

    def test_compare_name(self):
        """Compare names between TreeNodes"""
        nodes = self.TreeNode
        self.assertTrue(nodes["a"].compare_name(nodes["a"]))
        self.assertFalse(nodes["a"].compare_name(nodes["b"]))
        self.assertFalse(nodes["b"].compare_name(nodes["a"]))

    def test_compare_by_names(self):
        """Compare names between trees"""
        self.assertTrue(self.t.compare_by_names(self.t2))
        self.assertTrue(self.t.compare_by_names(self.t))
        self.assertFalse(self.t.compare_by_names(self.t4))

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
        f.name += 0.1
        self.assertNotEqual(f, g)

        # however, two TreeNodes that have no name should not compare equal
        f = TreeNode()
        g = TreeNode()
        self.assertNotEqual(f, g)

        f = TreeNode(name="foo")
        g = f.copy()
        self.assertNotEqual(f, g)

    def test_ne(self):
        """TreeNode should compare ne by id or data"""
        t, u, v = self.Comparisons
        self.assertFalse(t != t)
        self.assertTrue(t != u)

        f = TreeNode(name="foo")
        g = f.copy()
        self.assertTrue(f != g)

    def test_append(self):
        """TreeNode append should add item to end of self"""
        self.OneChild.append(TreeNode("c"))
        self.assertEqual(len(self.OneChild), 2)
        self.assertEqual(self.OneChild[-1].name, "c")
        self.OneChild.append(6)
        self.assertEqual(len(self.OneChild), 3)
        self.assertEqual(self.OneChild[-1].name, 6)
        # check that refs are updated when moved from one tree to another
        empty = TreeNode()
        empty.append(self.OneChild[-1])
        self.assertEqual(len(empty), 1)
        self.assertEqual(empty[0].name, 6)
        self.assertEqual(empty[0].parent, empty)
        self.assertEqual(self.OneChild[-1].name, "c")

    def test_extend(self):
        """TreeNode extend should add many items to end of self"""
        self.Empty.extend("abcdefgh")
        data = "".join([i.name for i in self.Empty])
        self.assertEqual(data, "abcdefgh")

    def test_insert(self):
        """TreeNode insert should insert item at specified index"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        parent.insert(3, 5)
        self.assertEqual(len(parent), 11)
        self.assertEqual(parent[3].name, 5)
        self.assertEqual(parent[4].name, "3")
        parent.insert(-1, 123)
        self.assertEqual(len(parent), 12)
        self.assertEqual(parent[-1].name, "9")
        self.assertEqual(parent[-2].name, 123)

    def test_pop(self):
        """TreeNode pop should remove and return child at specified index"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        last = parent.pop()
        assert last is nodes[-1]
        assert last.parent is None
        self.assertEqual(len(parent), 9)
        assert parent[-1] is nodes[-2]
        first = parent.pop(0)
        assert first is nodes[0]
        assert first.parent is None
        self.assertEqual(len(parent), 8)
        assert parent[0] is nodes[1]
        second_to_last = parent.pop(-2)
        assert second_to_last is nodes[-3]

    def test_remove(self):
        """TreeNode remove should remove first match by value, not id"""
        nodes = list(map(TreeNode, "abc" * 3))
        parent = TreeNode(children=nodes)
        self.assertEqual(len(parent), 9)
        parent.remove("a")
        self.assertEqual(len(parent), 8)
        self.assertEqual("".join([i.name for i in parent]), "bcabcabc")
        new_node = TreeNode("a")
        parent.remove(new_node)
        self.assertEqual(len(parent), 7)
        self.assertEqual("".join([i.name for i in parent]), "bcbcabc")

    def test_getitem(self):
        """TreeNode getitem should return item or slice"""
        r = self.TreeRoot
        n = self.TreeNode
        assert r[0] is n["b"]
        items = n["c"][0:1]
        self.assertEqual(len(items), 1)
        assert items[0] is n["d"]
        items = n["c"][0:2]
        self.assertEqual(len(items), 2)
        assert items[0] is n["d"]
        assert items[1] is n["e"]
        items = n["c"][:]
        self.assertEqual(len(items), 3)
        assert items[0] is n["d"]
        assert items[-1] is n["f"]

    def test_slice(self):
        """TreeNode slicing should return list, not TreeNode"""
        nodes = self.TreeNode
        c, d, e, f = nodes["c"], nodes["d"], nodes["e"], nodes["f"]
        assert c[:] is not c
        self.assertEqual(c[:], [d, e, f])
        self.assertEqual(c[1:2], [e])
        self.assertEqual(c[0:3:2], [d, f])

    def test_setitem(self):
        """TreeNode setitem should set item or extended slice of nodes"""
        parent, nodes = self.BigParent, self.BigName
        t = TreeNode(1)
        parent[0] = t
        assert parent[0] is t
        assert t.parent is parent
        assert nodes[0].parent is None

        u = TreeNode(2)
        parent[-2] = u
        assert parent[8] is u
        assert u.parent is parent
        assert nodes[8].parent is None

        parent[1:6:2] = "xyz"
        for i in [1, 3, 5]:
            assert nodes[i].parent is None
        self.assertEqual(parent[1].name, "x")
        self.assertEqual(parent[3].name, "y")
        self.assertEqual(parent[5].name, "z")
        for i in parent:
            assert i.parent is parent

    def test_setslice(self):
        """TreeNode setslice should set old-style slice of nodes"""
        parent, nodes = self.BigParent, self.BigName
        self.assertEqual(len(parent), 10)
        parent[5:] = []
        self.assertEqual(len(parent), 5)
        for i in range(5, 10):
            assert nodes[i].parent is None
        parent[1:3] = "abcd"
        self.assertEqual(len(parent), 7)
        for i in parent:
            assert i.parent is parent
        data_list = [i.name for i in parent]
        self.assertEqual(data_list, list("0abcd34"))
        parent[1:3] = parent[2:3]
        data_list = [i.name for i in parent]
        self.assertEqual(data_list, list("0bcd34"))

    def test_delitem(self):
        """TreeNode __delitem__ should delete item and set parent to None"""
        self.assertEqual(self.Child.parent, self.OneChild)
        self.assertEqual(len(self.OneChild), 1)
        del self.OneChild[0]
        self.assertEqual(self.OneChild.parent, None)
        self.assertEqual(len(self.OneChild), 0)

        nodes = self.BigName
        parent = self.BigParent
        self.assertEqual(len(parent), 10)
        for n in nodes:
            assert n.parent is parent
        del parent[-1]
        self.assertEqual(nodes[-1].parent, None)
        self.assertEqual(len(parent), 9)
        del parent[1:6:2]
        self.assertEqual(len(parent), 6)
        for i, n in enumerate(nodes):
            if i in [0, 2, 4, 6, 7, 8]:
                assert n.parent is parent
            else:
                assert n.parent is None

    def test_delslice(self):
        """TreeNode __delslice__ should delete items from start to end"""
        parent = self.BigParent
        nodes = self.BigName
        self.assertEqual(len(parent), 10)
        del parent[3:-2]
        self.assertEqual(len(parent), 5)
        for i, n in enumerate(nodes):
            if i in [3, 4, 5, 6, 7]:
                assert n.parent is None
            else:
                assert n.parent is parent

    def test_iter(self):
        """TreeNode iter should iterate over children"""
        r = self.TreeRoot
        n = self.TreeNode
        items = list(r)
        assert items[0] is n["b"]
        assert items[1] is n["h"]
        self.assertEqual(len(items), 2)

    def test_len(self):
        """TreeNode len should return number of children"""
        r = self.TreeRoot
        self.assertEqual(len(r), 2)

    def test_copy_recursive(self):
        """TreeNode.copy_recursive() should produce deep copy"""
        t = TreeNode(["t"])
        u = TreeNode(["u"])
        t.append(u)

        c = u.copy()
        assert c is not u
        assert c.name == u.name
        assert c.name is not u.name
        # note: name _is_ same object if it's immutable, e.g. a string.
        # deepcopy doesn't copy data for immutable objects.

        # need to check that we also copy arbitrary attributes
        t.XYZ = [3]
        c = t.copy()
        assert c is not t
        assert c[0] is not u
        assert c[0].name is not u.name
        assert c[0].name == u.name
        assert c.XYZ == t.XYZ
        assert c.XYZ is not t.XYZ

        t = self.TreeRoot
        c = t.copy()
        self.assertEqual(str(c), str(t))

    def test_copy(self):
        """TreeNode.copy() should work on deep trees"""
        t = comb_tree(1024)  # should break recursion limit on regular copy
        t.name = "foo"
        t.XYZ = [3]
        t2 = t.copy()
        t3 = t.copy()
        t3.name = "bar"

        self.assertEqual(len(t.tips()), 1024)
        self.assertEqual(len(t2.tips()), 1024)
        self.assertEqual(len(t3.tips()), 1024)

        self.assertIsNot(t, t2)
        self.assertEqual(t.name, t2.name)
        self.assertNotEqual(t.name, t3.name)

        self.assertEqual(t.XYZ, t2.XYZ)
        self.assertIsNot(t.XYZ, t2.XYZ)

        self.assertEqual(t.get_newick(), t2.get_newick())

        t_simple = TreeNode(["t"])
        u_simple = TreeNode(["u"])
        t_simple.append(u_simple)

        self.assertEqual(str(t_simple.copy()), str(t_simple.copy()))

    def test_copy_topology(self):
        """TreeNode.copy_topology() should produce deep copy ignoring attrs"""
        t = TreeNode(["t"])
        u = TreeNode(["u"])
        t.append(u)

        c = u.copy_topology()
        assert c is not u
        self.assertEqual(c.name, u.name)
        # note: name _is_ same object if it's immutable, e.g. a string.
        # deepcopy doesn't copy data for immutable objects.

        # need to check that we do not also copy arbitrary attributes
        t.XYZ = [3]
        c = t.copy_topology()
        assert c is not t
        assert c[0] is not u
        assert c[0].name is not u.name
        assert c[0].name == u.name
        assert not hasattr(c, "XYZ")

        t = self.TreeRoot
        c = t.copy()
        self.assertEqual(str(c), str(t))

    def _test_copy_copy(self):
        """copy.copy should raise TypeError on TreeNode"""
        t = TreeNode("t")
        u = TreeNode("u")
        t.append(u)
        self.assertRaises(TypeError, copy, t)
        self.assertRaises(TypeError, copy, u)

    def test_deepcopy(self):
        """copy.deepcopy should work on TreeNode"""
        t = TreeNode(["t"])
        u = TreeNode(["u"])
        t.append(u)

        c = deepcopy(u)
        assert c is not u
        assert c.name == u.name
        assert c.name is not u.name
        # note: name _is_ same object if it's immutable, e.g. a string.
        # deepcopy doesn't copy data for immutable objects.

        # need to check that we also copy arbitrary attributes
        t.XYZ = [3]
        c = deepcopy(t)
        assert c is not t
        assert c[0] is not u
        assert c[0].name is not u.name
        assert c[0].name == u.name
        assert c.XYZ == t.XYZ
        assert c.XYZ is not t.XYZ

        t = self.TreeRoot
        c = deepcopy(t)
        self.assertEqual(str(c), str(t))

    def test_Parent(self):
        """TreeNode parent should hold correct data and be mutable"""
        # check initial conditions
        self.assertEqual(self.Single.parent, None)
        # set parent and check parent/child relations
        self.Single.parent = self.Empty
        assert self.Single.parent is self.Empty
        self.assertEqual(self.Empty[0], self.Single)
        assert self.Single in self.Empty
        self.assertEqual(len(self.Empty), 1)
        # reset parent and check parent/child relations
        self.Single.parent = self.OneChild
        assert self.Single.parent is self.OneChild
        assert self.Single not in self.Empty
        assert self.Single is self.OneChild[-1]

        # following is added to check that we don't screw up when there are
        # nodes with different ids that still compare equal
        for i in self.Repeated:
            assert i.parent is self.Repeated
        last = self.Repeated[-1]
        last.parent = self.OneChild
        self.assertEqual(len(self.Repeated), 2)
        for i in self.Repeated:
            assert i.parent is self.Repeated
        assert last.parent is self.OneChild

    def test_index_in_parent(self):
        """TreeNode index_in_parent should hold correct data"""
        first = TreeNode("a")
        second = TreeNode("b")
        third = TreeNode("c")
        fourth = TreeNode("0", children=[first, second, third])
        self.assertEqual(len(fourth), 3)
        self.assertEqual(first.index_in_parent(), 0)
        self.assertEqual(second.index_in_parent(), 1)
        self.assertEqual(third.index_in_parent(), 2)
        del fourth[0]
        self.assertEqual(second.index_in_parent(), 0)
        self.assertEqual(third.index_in_parent(), 1)
        self.assertEqual(len(fourth), 2)
        assert first.parent is None

    def test_is_tip(self):
        """TreeNode is_tip should return True if node is a tip"""
        tips = "degh"
        for n in list(self.TreeNode.values()):
            if n.name in tips:
                self.assertEqual(n.is_tip(), True)
            else:
                self.assertEqual(n.is_tip(), False)

    def test_isRoot(self):
        """TreeNode isRoot should return True if parent is None"""
        r = "a"
        for n in list(self.TreeNode.values()):
            if n.name in r:
                self.assertEqual(n.is_root(), True)
            else:
                self.assertEqual(n.is_root(), False)

    def test_traverse(self):
        """TreeNode traverse should iterate over nodes in tree."""
        e = self.Empty
        s = self.Single
        o = self.OneChild
        m = self.Multi
        r = self.TreeRoot

        self.assertEqual([i.name for i in e.traverse()], [None])
        self.assertEqual([i.name for i in e.traverse(False, False)], [None])
        self.assertEqual([i.name for i in e.traverse(True, True)], [None])

        self.assertEqual([i.name for i in s.traverse()], ["a"])
        self.assertEqual([i.name for i in s.traverse(True, True)], ["a"])
        self.assertEqual([i.name for i in s.traverse(True, False)], ["a"])
        self.assertEqual([i.name for i in s.traverse(False, True)], ["a"])
        self.assertEqual([i.name for i in s.traverse(False, False)], ["a"])

        self.assertEqual([i.name for i in o.traverse()], ["a", "b"])
        self.assertEqual([i.name for i in o.traverse(True, True)], ["a", "b", "a"])
        self.assertEqual([i.name for i in o.traverse(True, False)], ["a", "b"])
        self.assertEqual([i.name for i in o.traverse(False, True)], ["b", "a"])
        self.assertEqual([i.name for i in o.traverse(False, False)], ["b"])

        self.assertEqual([i.name for i in m.traverse()], ["a", "b", "c", "d"])
        self.assertEqual(
            [i.name for i in m.traverse(True, True)], ["a", "b", "c", "d", "a"]
        )
        self.assertEqual(
            [i.name for i in m.traverse(True, False)], ["a", "b", "c", "d"]
        )
        self.assertEqual(
            [i.name for i in m.traverse(False, True)], ["b", "c", "d", "a"]
        )
        self.assertEqual([i.name for i in m.traverse(False, False)], ["b", "c", "d"])

        self.assertEqual(
            [i.name for i in r.traverse()], ["a", "b", "c", "d", "e", "f", "g", "h"]
        )
        self.assertEqual(
            [i.name for i in r.traverse(True, True)],
            ["a", "b", "c", "d", "e", "f", "g", "f", "c", "b", "h", "a"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(True, False)],
            ["a", "b", "c", "d", "e", "f", "g", "h"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(False, True)],
            ["d", "e", "g", "f", "c", "b", "h", "a"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(False, False)], ["d", "e", "g", "h"]
        )
        self.assertEqual(
            [i.name for i in r.traverse(True, True, False)],
            ["b", "c", "d", "e", "f", "g", "f", "c", "b", "h"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(True, False, False)],
            ["b", "c", "d", "e", "f", "g", "h"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(False, True, False)],
            ["d", "e", "g", "f", "c", "b", "h"],
        )
        self.assertEqual(
            [i.name for i in r.traverse(False, False, False)], ["d", "e", "g", "h"]
        )

        # this previously failed
        t = DndParser("((a:6,(b:1,c:2):8):12,(d:3,(e:1,f:1):4):10);")
        t0 = t.children[0]
        list(t0.traverse(self_before=False, self_after=True))
        list(t0.traverse(self_before=True, self_after=True))

    def test_levelorder(self):
        t = DndParser("(((A,B)C,(D,E)F,(G,H)I)J,(K,L)M)N;")
        exp = ["N", "J", "M", "C", "F", "I", "K", "L", "A", "B", "D", "E", "G", "H"]
        names = [n.name for n in t.levelorder()]
        self.assertEqual(names, exp)

    def test_ancestors(self):
        """TreeNode ancestors should provide list of ancestors, deepest first"""
        nodes, tree = self.TreeNode, self.TreeRoot
        self.assertEqual(nodes["a"].ancestors(), [])
        self.assertEqual(nodes["b"].ancestors(), [nodes["a"]])
        self.assertEqual(nodes["d"].ancestors(), nodes["f"].ancestors())
        self.assertEqual(
            nodes["g"].ancestors(), [nodes["f"], nodes["c"], nodes["b"], nodes["a"]]
        )

    def test_newick_with_labelled_nodes(self):
        """return newick with internal nodes labelled"""
        treestrings = (
            "(a,b,(c,(d,e)));",
            "(a:0.1,b:0.2,(c:0.3,(d:0.4,e:0.5):0.6):0.7);",
            "(a,b,(c,(d,e)edge.0)edge.1);",
        )
        expect = (
            "(a,b,(c,(d,e)edge.0)edge.1);",
            "(a:0.1,b:0.2,(c:0.3,(d:0.4,e:0.5)edge.0:0.6)edge.1:0.7);",
            "(a,b,(c,(d,e)edge.0)edge.1);",
        )
        for i, treestring in enumerate(treestrings):
            if i < 2:
                continue
            tree = make_tree(treestring=treestring)
            nwk = tree.get_newick(with_node_names=True, with_distances=True)
            self.assertEqual(nwk, expect[i])
            nwk = tree.get_newick_recursive(with_node_names=True, with_distances=True)
            self.assertEqual(nwk, expect[i])

    def test_root(self):
        """TreeNode root() should find root of tree"""
        nodes, root = self.TreeNode, self.TreeRoot
        for i in list(nodes.values()):
            assert i.root() is root

    def test_children(self):
        """TreeNode children should allow getting/setting children"""
        nodes = self.TreeNode
        for n in nodes:
            node = nodes[n]
            self.assertEqual(list(node), node.children)

        t = TreeNode(children="abc")
        self.assertEqual(len(t), 3)
        u, v = TreeNode("u"), TreeNode("v")

        # WARNING: If you set children directly, parent refs will _not_ update!
        t.children = [u, v]

        assert t[0] is u
        assert t[1] is v
        self.assertEqual(len(t), 2)

    def test_siblings(self):
        """TreeNode siblings() should return all siblings, not self"""
        self.assertEqual(self.Empty.siblings(), [])
        self.assertEqual(self.Child.siblings(), [])
        self.assertEqual(self.OneChild.siblings(), [])

        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(g.siblings(), [])
        self.assertEqual(f.siblings(), [d, e])
        self.assertEqual(e.siblings(), [d, f])
        self.assertEqual(d.siblings(), [e, f])
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
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(g.tips(), [])
        self.assertEqual(f.tips(), [g])
        self.assertEqual(e.tips(), [])
        self.assertEqual(d.tips(), [])
        self.assertEqual(c.tips(), [d, e, g])
        self.assertEqual(b.tips(), [d, e, g])
        self.assertEqual(h.tips(), [])
        self.assertEqual(a.tips(), [d, e, g, h])

    def test_itertips(self):
        """TreeNode itertips should iterate over terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.name for i in tree.iter_tips()], list("degh")),

    def test_nontips(self):
        """TreeNode nontips should return all non-terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.name for i in tree.nontips()], list("bcf"))

    def test_iterNonTips(self):
        """TreeNode iter_nontips should iterate over non-terminal descendants"""
        tree = self.TreeRoot
        self.assertEqual([i.name for i in tree.iter_nontips()], list("bcf"))

    def test_tip_children(self):
        """TreeNode tip_children should return all terminal children"""
        self.assertEqual(self.Empty.tip_children(), [])
        self.assertEqual(self.Child.tip_children(), [])
        self.assertEqual(self.OneChild.tip_children(), [self.Child])

        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(g.tip_children(), [])
        self.assertEqual(f.tip_children(), [g])
        self.assertEqual(e.tip_children(), [])
        self.assertEqual(d.tip_children(), [])
        self.assertEqual(c.tip_children(), [d, e])
        self.assertEqual(b.tip_children(), [])
        self.assertEqual(h.tip_children(), [])
        self.assertEqual(a.tip_children(), [h])

    def test_non_tip_children(self):
        """TreeNode non_tip_children should return all non-terminal children"""
        self.assertEqual(self.Empty.non_tip_children(), [])
        self.assertEqual(self.Child.non_tip_children(), [])
        self.assertEqual(self.OneChild.non_tip_children(), [])

        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(g.non_tip_children(), [])
        self.assertEqual(f.non_tip_children(), [])
        self.assertEqual(e.non_tip_children(), [])
        self.assertEqual(d.non_tip_children(), [])
        self.assertEqual(c.non_tip_children(), [f])
        self.assertEqual(b.non_tip_children(), [c])
        self.assertEqual(h.non_tip_children(), [])
        self.assertEqual(a.non_tip_children(), [b])

    def test_child_groups(self):
        """TreeNode child_groups should divide children by grandchild presence"""
        parent = TreeNode(children="aababbbaaabbbababbb")
        for node in parent:
            if node.name == "a":
                node.append("def")
        groups = parent.child_groups()
        self.assertEqual(len(groups), 10)
        exp_group_sizes = [2, 1, 1, 3, 3, 3, 1, 1, 1, 3]
        obs_group_sizes = [len(i) for i in groups]
        self.assertEqual(obs_group_sizes, exp_group_sizes)

        parent = TreeNode(children="aab")
        for node in parent:
            if node.name == "a":
                node.append("def")
        groups = parent.child_groups()
        self.assertEqual(len(groups), 2)
        self.assertEqual([len(i) for i in groups], [2, 1])

        parent = TreeNode(children="aaaaa")
        groups = parent.child_groups()
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), 5)

        parent = TreeNode(children="aaba")
        for node in parent:
            if node.name == "a":
                node.append("def")
        groups = parent.child_groups()
        self.assertEqual(len(groups), 3)
        self.assertEqual([len(i) for i in groups], [2, 1, 1])

    def test_remove_node(self):
        """TreeNode remove_node should delete node by id, not value"""
        parent = self.Repeated
        children = list(self.Repeated)
        self.assertEqual(len(parent), 3)
        self.assertEqual(parent.remove_node(children[1]), True)
        self.assertEqual(len(parent), 2)
        assert children[0].parent is parent
        assert children[1].parent is None
        assert children[2].parent is parent
        self.assertEqual(children[0].compare_name(children[1]), True)
        self.assertEqual(parent.remove_node(children[1]), False)
        self.assertEqual(len(parent), 2)
        self.assertEqual(parent.remove_node(children[0]), True)
        self.assertEqual(len(parent), 1)

    def test_lowest_common_ancestor(self):
        """TreeNode lowest_common_ancestor should return LCA for set of tips"""
        t1 = DndParser("((a,(b,c)d)e,f,(g,h)i)j;")
        t2 = t1.copy()
        t3 = t1.copy()
        t4 = t1.copy()
        input1 = ["a"]  # return self
        input2 = ["a", "b"]  # return e
        input3 = ["b", "c"]  # return d
        input4 = ["a", "h", "g"]  # return j
        exp1 = t1.get_node_matching_name("a")
        exp2 = t2.get_node_matching_name("e")
        exp3 = t3.get_node_matching_name("d")
        exp4 = t4
        obs1 = t1.lowest_common_ancestor(input1)
        obs2 = t2.lowest_common_ancestor(input2)
        obs3 = t3.lowest_common_ancestor(input3)
        obs4 = t4.lowest_common_ancestor(input4)
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)

        # verify multiple calls work
        t_mul = t1.copy()
        exp_1 = t_mul.get_node_matching_name("d")
        exp_2 = t_mul.get_node_matching_name("i")
        obs_1 = t_mul.lowest_common_ancestor(["b", "c"])
        obs_2 = t_mul.lowest_common_ancestor(["g", "h"])
        self.assertEqual(obs_1, exp_1)
        self.assertEqual(obs_2, exp_2)

    def test_lowest_common_ancestor_invalid_tips(self):
        """fail if tips not present"""
        t = DndParser("((a,(b,c)d)e,f,(g,h)i)j;")
        # no tips present in tree should raise exception
        with self.assertRaises(ValueError):
            t.lowest_common_ancestor(["m", "n"])

        # not all tips present in tree should raise exception
        with self.assertRaises(ValueError):
            t.lowest_common_ancestor(["a", "n"])

    def test_last_common_ancestor(self):
        """TreeNode last_common_ancestor should provide last common ancestor"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(a.last_common_ancestor(a), a)
        self.assertEqual(a.last_common_ancestor(b), a)
        self.assertEqual(a.last_common_ancestor(g), a)
        self.assertEqual(a.last_common_ancestor(h), a)

        self.assertEqual(b.last_common_ancestor(g), b)
        self.assertEqual(b.last_common_ancestor(d), b)
        self.assertEqual(b.last_common_ancestor(a), a)
        self.assertEqual(b.last_common_ancestor(h), a)

        self.assertEqual(d.last_common_ancestor(f), c)
        self.assertEqual(d.last_common_ancestor(g), c)
        self.assertEqual(d.last_common_ancestor(a), a)
        self.assertEqual(d.last_common_ancestor(h), a)

        self.assertEqual(g.last_common_ancestor(g), g)
        self.assertEqual(g.last_common_ancestor(f), f)
        self.assertEqual(g.last_common_ancestor(e), c)
        self.assertEqual(g.last_common_ancestor(c), c)
        self.assertEqual(g.last_common_ancestor(b), b)
        self.assertEqual(g.last_common_ancestor(a), a)
        self.assertEqual(g.last_common_ancestor(h), a)

        t = TreeNode("h")
        for i in [a, b, c, d, e, f, g, h]:
            self.assertEqual(i.last_common_ancestor(t), None)
            self.assertEqual(t.last_common_ancestor(i), None)

        u = TreeNode("a", children=[t])

    def test_separation(self):
        """TreeNode separation should return correct number of edges"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

        self.assertEqual(a.separation(a), 0)
        self.assertEqual(c.separation(c), 0)
        self.assertEqual(a.separation(b), 1)
        self.assertEqual(a.separation(h), 1)
        self.assertEqual(g.separation(h), 5)
        self.assertEqual(f.separation(d), 2)
        self.assertEqual(f.separation(c), 1)
        self.assertEqual(c.separation(f), 1)

    def test_name_unnamed_nodes(self):
        """name_unnamed_nodes assigns an arbitrary value when name == None"""
        tree, tree_nodes = self.TreeRoot, self.TreeNode
        tree_nodes["b"].name = "node2"
        tree_nodes["c"].name = None
        tree_nodes["f"].name = None
        tree_nodes["e"].name = "node3"
        tree.name_unnamed_nodes()
        self.assertEqual(tree_nodes["c"].name, "node1")
        self.assertEqual(tree_nodes["f"].name, "node4")

    def test_make_tree_array(self):
        """make_tree_array maps nodes to the descendants in them"""
        tree = self.TreeRoot
        result, node_list = tree.make_tree_array()
        assert_equal(
            result, array([[1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 1, 0], [0, 0, 1, 0]])
        )
        nodes = [node.name for node in node_list]
        self.assertEqual(nodes, ["a", "b", "c", "f"])
        # test if works with a dec_list supplied
        dec_list = ["d", "added", "e", "g", "h"]
        result2, node_list = tree.make_tree_array(dec_list)
        assert_equal(
            result2,
            array([[1, 0, 1, 1, 1], [1, 0, 1, 1, 0], [1, 0, 1, 1, 0], [0, 0, 0, 1, 0]]),
        )

    def test_get_node_names(self):
        """get_node_names works correctly"""
        tree = make_tree(treestring="((a:3,(b:2,(c:1,d:1):1):1):2,(e:3,f:3):2);")
        names = tree.get_node_names(includeself=False, tipsonly=False)
        self.assertTrue(tree.name not in names)
        names = tree.get_node_names(includeself=True, tipsonly=False)
        self.assertTrue(tree.name in names)
        tree.get_node_matching_name("a")

    def test_reassign_names(self):
        """reassign_names should rename node names based on dict mapping"""
        t = self.TreeRoot
        mapping = dict([(x, str(i)) for i, x in enumerate("abfg")])
        exp_names = ["0", "1", "2", "3", "c", "d", "e", "h"]
        t.reassign_names(mapping)
        obs_names = sorted(t.get_node_names())
        self.assertEqual(obs_names, exp_names)

    def test_reassign_names_specific_nodes(self):
        """reassign_names should rename nodes based on dict mapping"""
        t = self.TreeRoot
        nodes = [self.TreeNode["a"], self.TreeNode["b"]]
        mapping = dict([(x, str(i)) for i, x in enumerate("abfg")])
        exp_names = ["0", "1", "c", "d", "e", "f", "g", "h"]
        t.reassign_names(mapping, nodes)
        obs_names = sorted(t.get_node_names())
        self.assertEqual(obs_names, exp_names)

    def test_get_nodes_dict(self):
        """get_nodes_dict returns a dict keyed by name, value is node"""
        t = self.TreeRoot
        nodes = self.TreeNode
        self.assertEqual(t.get_nodes_dict(), nodes)

    def test_get_nodes_dict_nonunique_names(self):
        """get_nodes_dict raises if non unique names are in tree"""
        t = self.TreeRoot
        t.children[0].name = "same"
        t.children[0].children[0].name = "same"
        self.assertRaises(TreeError, t.get_nodes_dict)

    def test_remove_deleted(self):
        """remove_deleted should remove all nodes where is_deleted tests true."""
        tree = DndParser(
            "((a:3,(b:2,(c:1,d:1):1):1):2,(e:3,f:3):2);", constructor=TreeNode
        )
        result_not_deleted = deepcopy(tree)
        tree.remove_deleted(lambda x: x.name in [])
        self.assertEqual(str(tree), str(result_not_deleted))
        deleted = set(["b", "d", "e", "f"])
        result_tree = DndParser("((a:3,((c:1):1):1):2);", constructor=TreeNode)
        is_deleted = lambda x: x.name in deleted
        tree.remove_deleted(is_deleted)
        self.assertEqual(str(tree), str(result_tree))

    def test_prune(self):
        """prune should reconstruct correct topology of tree."""
        tree = DndParser("((a:3,((c:1):1):1):2);", constructor=TreeNode)
        tree.prune()
        result_tree = DndParser("((a:3,c:1));", constructor=TreeNode)
        self.assertEqual(str(tree), str(result_tree))

        samename_bug = DndParser("((A,B)SAMENAME,((C,D)SAMENAME));")
        samename_bug.prune()
        exp_tree_str = "((A,B)SAMENAME,(C,D)SAMENAME);"
        self.assertEqual(str(samename_bug), exp_tree_str)

    def test_get_node_matching_name(self):
        """TreeNode get_node_matching_name should return node that matches name"""
        nodes = self.TreeNode
        root = self.TreeRoot
        assert root.get_node_matching_name("g") is nodes["g"]

    def test_subset(self):
        """subset should return set of leaves that descends from node"""
        t = self.t
        self.assertEqual(t.subset(), frozenset("HGRM"))
        c = t.children[0]
        self.assertEqual(c.subset(), frozenset("HG"))
        leaf = c.children[1]
        self.assertEqual(leaf.subset(), frozenset(""))

    def test_subsets(self):
        """subsets should return all subsets descending from a set"""
        t = self.t
        self.assertEqual(t.subsets(), frozenset([frozenset("HG"), frozenset("RM")]))

    def test_compare_by_subsets(self):
        """compare_by_subsets should return the fraction of shared subsets"""
        result = self.t.compare_by_subsets(self.t)
        self.assertEqual(result, 0)

        result = self.t2.compare_by_subsets(self.t2)
        self.assertEqual(result, 0)

        result = self.t.compare_by_subsets(self.t2)
        self.assertEqual(result, 0.5)

        result = self.t.compare_by_subsets(self.t4)
        self.assertEqual(result, 1 - 2.0 / 5)

        result = self.t.compare_by_subsets(self.t4, exclude_absent_taxa=True)
        self.assertEqual(result, 1 - 2.0 / 3)

        result = self.t.compare_by_subsets(self.TreeRoot, exclude_absent_taxa=True)
        self.assertEqual(result, 1)

        result = self.t.compare_by_subsets(self.TreeRoot)
        self.assertEqual(result, 1)


class PhyloNodeTests(TestCase):
    """Tests of phylogeny-specific methods."""

    def setUp(self):
        """Creates a standard tree"""
        nodes = dict([(x, PhyloNode(x)) for x in "abcdefgh"])
        nodes["a"].append(nodes["b"])
        nodes["b"].append(nodes["c"])
        nodes["c"].append(nodes["d"])
        nodes["c"].append(nodes["e"])
        nodes["c"].append(nodes["f"])
        nodes["f"].append(nodes["g"])
        nodes["a"].append(nodes["h"])
        self.TreeNode = nodes
        self.TreeRoot = nodes["a"]
        nodes["a"].length = None
        nodes["b"].length = 0
        nodes["c"].length = 3
        nodes["d"].length = 1
        nodes["e"].length = 4
        nodes["f"].length = 2
        nodes["g"].length = 3
        nodes["h"].length = 2

        self.s = "((H:1,G:1):2,(R:0.5,M:0.7):3);"
        self.t = DndParser(self.s, PhyloNode)
        self.s3 = "(((H:1,G:1,O:1):2,R:3):1,X:4);"
        self.t3 = DndParser(self.s3, PhyloNode)

    def test_init(self):
        """Check PhyloNode constructor"""
        n = PhyloNode("foo", length=10)
        self.assertEqual(n.name, "foo")
        self.assertEqual(n.length, 10)

        n = PhyloNode("bar")
        self.assertEqual(n.name, "bar")
        self.assertEqual(n.length, None)

        n = PhyloNode()
        self.assertEqual(n.name, None)
        self.assertEqual(n.length, None)

    def test_total_descending_branch_length(self):
        """total_descending_branch_length returns total branchlength below self"""
        t = self.TreeRoot
        exp = 15
        obs = t.total_descending_branch_length()
        self.assertEqual(obs, exp)

        node_c = self.TreeNode["c"]
        exp = 10
        obs = node_c.total_descending_branch_length()
        self.assertEqual(obs, exp)

    def test_total_length(self):
        """total_length returns total branchlength, irrespective of edge"""
        t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2);"
        t = DndParser(t_str, constructor=PhyloNode)
        self.assertEqual(t.total_length(), 21)
        node = t.get_node_matching_name("F")
        length = node.total_length()
        self.assertEqual(length, 21)

    def test_tips_within_distance(self):
        """tips_within_distance returns tips that are within distance from self"""
        t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2)L;"
        t = DndParser(t_str, constructor=PhyloNode)
        nodes = t.get_nodes_dict()
        e_node = nodes["E"]

        exp_at_dist_2 = []
        exp_at_dist_3 = ["A", "C", "D"]
        exp_at_dist_4 = ["A", "B", "C", "D", "F"]

        obs_at_dist_2 = sorted([n.name for n in e_node.tips_within_distance(2)])
        obs_at_dist_3 = sorted([n.name for n in e_node.tips_within_distance(3)])
        obs_at_dist_4 = sorted([n.name for n in e_node.tips_within_distance(4)])

        self.assertEqual(obs_at_dist_2, exp_at_dist_2)
        self.assertEqual(obs_at_dist_3, exp_at_dist_3)
        self.assertEqual(obs_at_dist_4, exp_at_dist_4)

    def test_tips_within_distance_nodistances(self):
        """tips_within_distance returns tips that are within distance from self"""
        t_str = "(A,B,(C,D)E,(F,((G,H)I)J)K)L;"
        t = DndParser(t_str, constructor=PhyloNode)
        nodes = t.get_nodes_dict()
        e_node = nodes["E"]

        exp = sorted([n.name for n in t.tips()])
        obs = sorted([n.name for n in e_node.tips_within_distance(0)])
        self.assertEqual(obs, exp)

    def test_distance(self):
        """PhyloNode Distance should report correct distance between nodes"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        b = nodes["b"]
        c = nodes["c"]
        d = nodes["d"]
        e = nodes["e"]
        f = nodes["f"]
        g = nodes["g"]
        h = nodes["h"]

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

    def test_compare_by_tip_distances(self):
        obs = self.t.compare_by_tip_distances(self.t3)
        # note: common taxa are H, G, R (only)
        m1 = array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = correlation(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1 - r) / 2)

    def test_compare_by_tip_distances_sample(self):
        obs = self.t.compare_by_tip_distances(self.t3, sample=3, shuffle_f=sorted)
        # note: common taxa are H, G, R (only)
        m1 = array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = correlation(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1 - r) / 2)

        # 4 common taxa, still picking H, G, R
        s = "((H:1,G:1):2,(R:0.5,M:0.7,Q:5):3);"
        t = DndParser(self.s, PhyloNode)
        s3 = "(((H:1,G:1,O:1):2,R:3,Q:10):1,X:4);"
        t3 = DndParser(self.s3, PhyloNode)
        obs = t.compare_by_tip_distances(t3, sample=3, shuffle_f=sorted)

    def test_tip_to_tip_distances_endpoints(self):
        """Test getting specifc tip distances  with tip_to_tip_distances"""
        exp_nodes = [
            self.t.get_node_matching_name("H"),
            self.t.get_node_matching_name("G"),
            self.t.get_node_matching_name("M"),
        ]
        names = ["H", "G", "M"]
        exp_dists = array([[0, 2.0, 6.7], [2.0, 0, 6.7], [6.7, 6.7, 0.0]])
        got_dists, got_nodes = self.t.tip_to_tip_distances(endpoints=names)
        assert_equal(got_dists, exp_dists)
        self.assertEqual(got_nodes, exp_nodes)

    def test_prune(self):
        """prune should reconstruct correct topology and Lengths of tree."""
        tree = DndParser("((a:3,((c:1):1):1):2);", constructor=PhyloNode)
        tree.prune()
        result_tree = DndParser("((a:3.0,c:3.0):2.0);", constructor=PhyloNode)
        self.assertEqual(str(tree), str(result_tree))

    def test_str(self):
        """PhyloNode str should give expected results"""
        nodes, tree = self.TreeNode, self.TreeRoot
        a = nodes["a"]
        c = nodes["c"]
        f = nodes["f"]
        h = nodes["h"]

        self.assertEqual(str(h), "h:2;")
        self.assertEqual(str(f), "(g:3)f:2;")
        self.assertEqual(str(a), "(((d:1,e:4,(g:3)f:2)c:3)b:0,h:2)a;")
        # check that None isn't converted any more
        h.length = None
        c.length = None  # need to test both leaf and internal node
        self.assertEqual(str(a), "(((d:1,e:4,(g:3)f:2)c)b:0,h)a;")

    def test_get_max_tip_tip_distance(self):
        """get_max_tip_tip_distance should get max tip distance across tree"""
        nodes, tree = self.TreeNode, self.TreeRoot
        dist, names, node = tree.get_max_tip_tip_distance()
        self.assertEqual(dist, 15.0)  # due to nodes with single descendents!!
        self.assertEqual(sorted(names), ["e", "g"])
        self.assertEqual(node.name, "b")

    def test_set_max_tip_tip_distance(self):
        """set_max_tip_tip_distance sets MaxDistTips across tree"""
        nodes, tree = self.TreeNode, self.TreeRoot
        tree.set_max_tip_tip_distance()
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 10)
        self.assertEqual(sorted([tip_a[1], tip_b[1]]), ["g", "h"])

    def test_max_tip_tip_distance(self):
        """max_tip_tip_distance returns the max dist between any pair of tips"""
        nodes, tree = self.TreeNode, self.TreeRoot
        max_dist, tip_pair = tree.max_tip_tip_distance()
        self.assertEqual(max_dist, 10)
        try:
            self.assertEqual(tip_pair, ("h", "g"))
        except AssertionError:
            self.assertEqual(tip_pair, ("g", "h"))

    def test__find_midpoint_nodes(self):
        """_find_midpoint_nodes should return nodes surrounding the midpoint"""
        nodes, tree = self.TreeNode, self.TreeRoot
        max_dist = 10
        tip_pair = ("g", "h")
        result = tree._find_midpoint_nodes(max_dist, tip_pair)
        self.assertEqual(result, (nodes["b"], nodes["c"]))
        tip_pair = ("h", "g")
        result = tree._find_midpoint_nodes(max_dist, tip_pair)
        self.assertEqual(result, (nodes["f"], nodes["c"]))

    def test_root_at_midpoint(self):
        """root_at_midpoint performs midpoint rooting"""
        nodes, tree = self.TreeNode, self.TreeRoot
        # works when the midpoint falls on an existing edge
        tree1 = deepcopy(tree)
        result = tree1.root_at_midpoint()
        self.assertEqual(result.distance(result.get_node_matching_name("e")), 4)
        self.assertEqual(result.get_distances(), tree1.get_distances())
        # works when the midpoint falls between two existing edges
        nodes["f"].length = 1
        nodes["c"].length = 4
        result = tree.root_at_midpoint()
        self.assertEqual(result.distance(result.get_node_matching_name("e")), 5.0)
        self.assertEqual(result.distance(result.get_node_matching_name("g")), 5.0)
        self.assertEqual(result.distance(result.get_node_matching_name("h")), 5.0)
        self.assertEqual(result.distance(result.get_node_matching_name("d")), 2.0)
        self.assertEqual(result.get_distances(), tree.get_distances())

    def test_root_at_midpoint2(self):
        """root_at_midpoint works when midpoint is on both sides of root"""
        # also checks whether it works if the midpoint is adjacent to a tip
        nodes, tree = self.TreeNode, self.TreeRoot
        nodes["h"].length = 20
        result = tree.root_at_midpoint()
        self.assertEqual(result.distance(result.get_node_matching_name("h")), 14)
        self.assertEqual(result.get_distances(), tree.get_distances())

    def test_root_at_midpoint3(self):
        """midpoint between nodes should behave correctly"""
        tree = DndParser("(a:1,((c:1,d:2.5)n3:1,b:1)n2:1)rt;")
        tmid = tree.root_at_midpoint()
        self.assertEqual(tmid.get_distances(), tree.get_distances())
        tree.get_tip_names()
        [t.name for t in tree.nontips()]
        self.assertTrue(tmid.is_root())
        self.assertEqual(tmid.distance(tmid.get_node_matching_name("d")), 2.75)

    def test_root_at_midpoint4(self):
        """midpoint should be selected correctly when it is an internal node"""
        tree = DndParser("(a:1,((c:1,d:3)n3:1,b:1)n2:1)rt;")
        tmid = tree.root_at_midpoint()
        self.assertEqual(tmid.get_distances(), tree.get_distances())
        tree.get_tip_names()
        [t.name for t in tree.nontips()]
        # for tipname in tipnames:
        #     tmid_tip = tmid.get_node_matching_name(tipname)
        #     orig_tip = tree.get_node_matching_name(tipname)
        #     for nontipname in nontipnames:
        #         tmid_dist=\
        #           tmid.get_node_matching_name(nontipname).distance(tmid_tip)
        #         orig_dist=\
        #           tree.get_node_matching_name(nontipname).distance(orig_tip)
        #         print nontipname, tipname, 'assert'
        # self.assertEqual(tmid_dist, orig_dist)
        self.assertTrue(tmid.is_root())
        self.assertEqual(tmid.distance(tmid.get_node_matching_name("d")), 3)

    def test_root_at_midpoint5(self):
        """midpoint should be selected correctly when on an even 2tip tree"""
        tree = DndParser("""(BLO_1:0.649351,BLO_2:0.649351):0.0;""")
        tmid = tree.root_at_midpoint()
        self.assertEqual(tmid.get_distances(), tree.get_distances())
        tree.get_tip_names()
        [t.name for t in tree.nontips()]

        self.assertTrue(tmid.is_root())
        assert_allclose(tmid.distance(tmid.get_node_matching_name("BLO_2")), 0.649351)
        assert_allclose(tmid.distance(tmid.get_node_matching_name("BLO_1")), 0.649351)
        assert_allclose(tmid[0].distance(tmid[1]), 2.0 * 0.649351)

    def test_set_tip_distances(self):
        """set_tip_distances should correctly set tip distances."""
        tree = DndParser(
            "(((A1:.1,B1:.1):.1,(A2:.1,B2:.1):.1):.3,((A3:.1,B3:.1):.1,(A4:.1,B4:.1):.1):.3);",
            constructor=PhyloNode,
        )

        # expected distances for a post order traversal
        expected_tip_distances = [
            0,
            0,
            0.1,
            0,
            0,
            0.1,
            0.2,
            0,
            0,
            0.1,
            0,
            0,
            0.1,
            0.2,
            0.5,
        ]
        # tips should have distance of 0
        tree.set_tip_distances()
        for node in tree.tips():
            self.assertEqual(node.TipDistance, 0)
        idx = 0
        for node in tree.traverse(self_before=False, self_after=True):
            self.assertEqual(node.TipDistance, expected_tip_distances[idx])
            idx += 1

    def test_scale_branch_lengths(self):
        """scale_branch_lengths should correclty scale branch lengths."""
        tree = DndParser(
            "(((A1:.1,B1:.1):.1,(A2:.1,B2:.1):.1):.3,((A3:.1,B3:.1):.1,(A4:.1,B4:.1):.1):.3);",
            constructor=PhyloNode,
        )
        tree.scale_branch_lengths(max_length=100, ultrametric=True)
        expected_tree = "(((A1:20,B1:20):20,(A2:20,B2:20):20):60,((A3:20,B3:20):20,(A4:20,B4:20):20):60);"
        self.assertEqual(str(tree), expected_tree)

    def test_unrooted(self):
        """unrooted should preserve tips, drop a node"""
        rooted = make_tree(treestring="(B:0.2,(C:0.2,D:0.2)F:0.2)G;")
        unrooted = rooted.unrooted()
        self.assertEqual(
            sorted(rooted.get_tip_names()), sorted(unrooted.get_tip_names())
        )
        self.assertLess(len(unrooted.get_node_names()), len(rooted.get_node_names()))

    def test_get_figure(self):
        """exercising get_figure"""
        t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2)L;"
        t = DndParser(t_str, constructor=PhyloNode)
        _ = t.get_figure(style="square")


class _tip_tip_distances_I:
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
        self.assertEqual([i.name for i in order], list("abc"))
        assert_equal(matrix, array([[0, 3, 4], [3, 0, 5], [4, 5, 0]]))

    def test_two_level(self):
        """tip_to_tip should work for two-level tree"""
        matrix, order = self.fun(self.root_two_level)
        self.assertEqual([i.name for i in order], list("abcd"))
        assert_allclose(
            matrix,
            array([[0, 3, 4, 1.4], [3, 0, 5, 2.4], [4, 5, 0, 3.4], [1.4, 2.4, 3.4, 0]]),
        )


class Test_tip_tip_distances_array(_tip_tip_distances_I, TestCase):
    """Tests for the array implementation of tip_to_tip distances"""

    def setUp(self):
        """Specify which method to call."""
        self.fun = lambda x: x.tip_to_tip_distances()
        super(Test_tip_tip_distances_array, self).setUp()

    def test_std(self):
        """tip_to_tip should work for small but complex tree"""
        dist, tips = self.fun(self.root_std)
        tips = [tip.name for tip in tips]
        assert_equal(dist, tree_std_dist)
        assert_equal(tips, tree_std_tips)

    def test_one_child(self):
        """tip_to_tip should work for tree with a single child"""
        dist, tips = self.fun(self.root_one_child)
        tips = [tip.name for tip in tips]
        assert_equal(dist, tree_one_child_dist)
        assert_equal(tips, tree_one_child_tips)


# for use with testing iterative copy method


def comb_tree(num_leaves):
    """Returns a comb node_class tree."""
    branch_child = 1

    root = TreeNode()
    curr = root

    for i in range(num_leaves - 1):
        curr.children[:] = [TreeNode(parent=curr), TreeNode(parent=curr)]
        curr = curr.children[branch_child]
    return root


# Moved  from test_tree2.py during code sprint on 04/14/10
# Missing tests: edge attributes (name, length, children) only get tested
# in passing by some of these tests.  See also xxx's


class TreeInterfaceForLikelihoodFunction(TestCase):
    default_newick = "((A:1,B:2)ab:3,((C:4,D:5)cd,E:6)cde:7)"

    def _maketree(self, treestring=None):
        if treestring is None:
            treestring = self.default_newick
        return make_tree(treestring=treestring, underscore_unmunge=True)

    def setUp(self):
        self.default_tree = self._maketree()

    def test_get_edge_names(self):
        tree = self._maketree()
        for (a, b, outgroup, result) in [
            ("A", "B", None, ["A", "B"]),
            ("E", "C", None, ["C", "D", "cd", "E"]),
            ("C", "D", "E", ["C", "D"]),
        ]:
            self.assertEqual(tree.get_edge_names(a, b, True, False, outgroup), result)

    def test_parser(self):
        """nasty newick"""
        nasty = "( (A :1.0,'B (b)': 2) [com\nment]pair:3,'longer name''s':4)dash_ed;"
        nice = "((A:1.0,'B (b)':2.0)pair:3.0,'longer name''s':4.0)dash_ed;"
        tree = self._maketree(nasty)
        tidied = tree.get_newick(with_distances=1)
        self.assertEqual(tidied, nice)

    # Likelihood Function Interface

    def test_get_edge_names(self):
        tree = self.default_tree
        clade = tree.get_edge_names("C", "E", stem=0, clade=1)
        clade.sort()
        self.assertEqual(clade, ["C", "D", "E", "cd"])

        all = tree.get_edge_names("C", "E", stem=1, clade=1)
        all.sort()
        self.assertEqual(all, ["C", "D", "E", "cd", "cde"])

        stem = tree.get_edge_names("C", "E", stem=1, clade=0)
        self.assertEqual(stem, ["cde"])

    def test_get_edge_namesUseOutgroup(self):
        t1 = make_tree(treestring="((A,B)ab,(F,(C,D)cd)cdf,E)root;")
        # a, e, ogroup f
        t2 = make_tree(treestring="((E,(A,B)ab)abe,F,(C,D)cd)root;")
        expected = ["A", "B", "E", "ab"]
        for t in [t1, t2]:
            edges = t.get_edge_names(
                "A", "E", stem=False, clade=True, outgroup_name="F"
            )
            edges.sort()
            self.assertEqual(expected, edges)

    def test_get_connecting_node(self):
        tree = self.default_tree
        self.assertEqual(tree.get_connecting_node("A", "B").name, "ab")
        self.assertEqual(tree.get_connecting_node("A", "C").name, "root")

    def test_get_connecting_edges(self):
        """correctly identify connecting edges"""
        tree = make_tree(treestring="(((Human,HowlerMon)a,Mouse)b,NineBande,DogFaced);")
        edges = [e.name for e in tree.get_connecting_edges("Human", "Mouse")]
        self.assertEqual(set(edges), set(["Human", "Mouse", "a"]))

        edges = [e.name for e in tree.get_connecting_edges("b", "Human")]
        self.assertEqual(set(edges), set(["Human", "a", "b"]))

    def test_get_node_matching_name(self):
        tree = self.default_tree
        for (name, expect_tip) in [("A", True), ("ab", False)]:
            edge = tree.get_node_matching_name(name)
            self.assertEqual(edge.name, name)
            self.assertEqual(edge.istip(), expect_tip)

    def test_get_edge_vector(self):
        """correctly return vector of edges from a tree"""
        tree = self.default_tree
        names = [e.name for e in tree.get_edge_vector()]
        self.assertEqual(names, ["A", "B", "ab", "C", "D", "cd", "E", "cde", "root"])

        names = [e.name for e in tree.get_edge_vector(include_root=False)]
        self.assertEqual(names, ["A", "B", "ab", "C", "D", "cd", "E", "cde"])

    def test_get_newick_recursive(self):
        orig = "((A:1.0,B:2.0)ab:3.0,((C:4.0,D:5.0)cd:6.0,E:7.0)cde:8.0)all;"
        unlen = "((A,B)ab,((C,D)cd,E)cde)all;"
        tree = self._maketree(orig)
        self.assertEqual(tree.get_newick_recursive(with_distances=1), orig)
        self.assertEqual(tree.get_newick_recursive(), unlen)

        tree.name = "a'l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a'l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a''l';"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), ugly_name)

        tree.name = "a_l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a_l';"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), ugly_name)

        tree.name = "a l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), ugly_name)

        tree.name = "'a l'"
        quoted_name = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        quoted_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), quoted_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), quoted_name)

    def test_get_newick(self):
        orig = "((A:1.0,B:2.0)ab:3.0,((C:4.0,D:5.0)cd:6.0,E:7.0)cde:8.0)all;"
        unlen = "((A,B)ab,((C,D)cd,E)cde)all;"
        tree = self._maketree(orig)
        self.assertEqual(tree.get_newick(with_distances=1), orig)
        self.assertEqual(tree.get_newick(), unlen)

        tree.name = "a'l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a'l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a''l';"
        self.assertEqual(tree.get_newick(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick(escape_name=False), ugly_name)

        tree.name = "a_l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a_l';"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), ugly_name)

        tree.name = "a l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        self.assertEqual(tree.get_newick_recursive(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.get_newick_recursive(escape_name=False), ugly_name)

        tree.name = "'a l'"
        quoted_name = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        quoted_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        self.assertEqual(tree.get_newick(escape_name=True), quoted_name_esc)
        self.assertEqual(tree.get_newick(escape_name=False), quoted_name)

    def test_XML(self):
        # should add some non-length parameters
        orig = self.default_tree
        xml = orig.get_xml()
        parsed = make_tree(treestring=xml)
        self.assertEqual(str(orig), str(parsed))

    # Magic methods

    def test_str(self):
        """testing (well, exercising at least), __str__"""
        str(self.default_tree)

    def test_repr(self):
        """testing (well, exercising at least), __repr__"""
        repr(self.default_tree)

    def test_eq(self):
        """testing (well, exercising at least), __eq__"""
        # xxx not good enough!
        t1 = self._maketree()
        t2 = self._maketree()
        self.assertTrue(t1 == t1)
        self.assertFalse(t1 == t2)

    def test_balanced(self):
        """balancing an unrooted tree"""
        t = make_tree(treestring="((a,b),((c1,(c2,(c3,(c4,(c5,(c6,c7)))))),(d,e)),f)")
        b = make_tree(treestring="(c1,(c2,(c3,(c4,(c5,(c6,c7))))),((d,e),((a,b),f)))")
        self.assertEqual(str(t.balanced()), str(b))

    def test_params_merge(self):
        t = make_tree(treestring="((((a,b)ab,c)abc),d)")
        for (label, length, beta) in [("a", 1, 20), ("b", 3, 2.0), ("ab", 4, 5.0)]:
            t.get_node_matching_name(label).params = {"length": length, "beta": beta}
        t = t.get_sub_tree(["b", "c", "d"])
        self.assertEqual(
            t.get_node_matching_name("b").params,
            {"length": 7, "beta": float(2 * 3 + 4 * 5) / (3 + 4)},
        )
        self.assertRaises(ValueError, t.get_sub_tree, ["b", "c", "xxx"])
        self.assertEqual(
            str(t.get_sub_tree(["b", "c", "xxx"], ignore_missing=True)), "(b:7,c)root;"
        )

    def test_making_from_list(self):
        tipnames_with_spaces = ["a_b", "a b", "T'lk"]
        tipnames_with_spaces.sort()
        t = make_tree(tip_names=tipnames_with_spaces)
        result = t.get_tip_names()
        result.sort()
        assert result == tipnames_with_spaces

    def test_getset_param_value(self):
        """test getting, setting of param values"""
        t = make_tree(treestring="((((a:.2,b:.3)ab:.1,c:.3)abc:.4),d:.6)")
        self.assertEqual(t.get_param_value("length", "ab"), 0.1, 2)
        t.set_param_value("zz", "ab", 4.321)
        node = t.get_node_matching_name("ab")
        self.assertEqual(4.321, node.params["zz"], 4)


class SmallTreeReshapeTestClass(TestCase):
    def test_rootswaps(self):
        """testing (well, exercising at least), unrooted"""
        new_tree = make_tree(treestring="((a,b),(c,d))")
        new_tree = new_tree.unrooted()
        self.assertTrue(len(new_tree.children) > 2, "not unrooted right")

    def test_reroot(self):
        tree = make_tree(treestring="((a,b),(c,d),e)")
        tree2 = tree.rooted_with_tip("b")
        self.assertEqual(tree2.get_newick(), "(a,b,((c,d),e));")

    def test_same_shape(self):
        """test topology assessment"""
        t1 = make_tree(treestring="(((s1,s5),s3),s2,s4);")
        t2 = make_tree(treestring="((s1,s5),(s2,s4),s3);")
        t3 = make_tree(treestring="((s1,s4),(s2,s5),s3);")
        assert t1.same_topology(t2), (t1, t2)
        assert not t1.same_topology(t3), (t1, t3)
        assert not t2.same_topology(t3), (t2, t3)


# =============================================================================
# these are tests involving tree manipulation methods
# hence, testing them for small and big trees
# the tests are written once for the small tree, the big tree
# tests are performed by inheriting from this class, but over-riding
# the setUp.


class TestTree(TestCase):
    """tests for a single tree-type"""

    def setUp(self):
        self.name = "small tree - "
        self.otu_names = ["NineBande", "Mouse", "HowlerMon", "DogFaced"]
        self.otu_names.sort()
        self.newick = "(((Human,HowlerMon),Mouse),NineBande,DogFaced);"
        self.newick_sorted = "(DogFaced,((HowlerMon,Human),Mouse),NineBande);"
        self.newick_reduced = "((HowlerMon,Mouse),NineBande,DogFaced);"
        self.tree = make_tree(treestring=self.newick)

    def test_sorttree(self):
        """testing (well, exercising at least) treesort"""
        new_tree = self.tree.sorted()
        if hasattr(self, "newick_sorted"):
            self.assertEqual(self.newick_sorted, new_tree.get_newick(with_distances=0))

    def test_getsubtree(self):
        """testing getting a subtree"""
        subtree = self.tree.unrooted().get_sub_tree(self.otu_names)

        new_tree = make_tree(treestring=self.newick_reduced).unrooted()

        # check we get the same names
        self.assertEqual(*[len(t.children) for t in (subtree, new_tree)])
        self.assertEqual(str(subtree), str(new_tree))

    def test_getsubtree_2(self):
        """tree.get_sub_tree() has same pairwise tip dists as tree (len0 node)"""
        t1 = DndParser(
            "((a:1,b:2):4,((c:3, j:17.2):0,(d:1,e:1):2):3)", PhyloNode
        )  # note c,j is len 0 node
        orig_dists = t1.get_distances()
        subtree = t1.get_sub_tree(set(["a", "b", "d", "e", "c"]))
        sub_dists = subtree.get_distances()
        for pair, dist in list(sub_dists.items()):
            self.assertEqual((pair, dist), (pair, orig_dists[pair]))

    def test_getsubtree_3(self):
        """tree.get_sub_tree() has same pairwise tip dists as tree

        (nonzero nodes)
        """
        t1 = DndParser(
            "((a:1,b:2):4,((c:3, j:17):0,(d:1,e:1):2):3)", PhyloNode
        )  # note c,j is len 0 node
        orig_dists = t1.get_distances()
        subtree = t1.get_sub_tree(set(["a", "b", "d", "e", "c"]))
        subtree.get_distances()
        # for pair, dist in sub_dists.items():
        # self.assertEqual((pair,dist), (pair,orig_dists[pair]))
        t2 = DndParser(
            "((a:1,b:2):4,((c:2, j:16):1,(d:1,e:1):2):3)", PhyloNode
        )  # note c,j similar to above
        t2_dists = t2.get_distances()
        # ensure t2 is same as t1, except j->c or c->j
        for pair, dist in list(t2_dists.items()):
            if (pair == ("c", "j")) or (pair == ("j", "c")):
                continue
            self.assertEqual((pair, dist), (pair, orig_dists[pair]))
        sub2 = t2.get_sub_tree(set(["a", "b", "d", "e", "c"]))
        sub2_dists = sub2.get_distances()
        for pair, dist in list(sub2_dists.items()):
            self.assertEqual((pair, dist), (pair, orig_dists[pair]))

    def test_getsubtree_4(self):
        """tree.get_sub_tree() handles keep_root correctly"""
        t1 = DndParser("((a:1,b:2):4,(((c:2)cparent:1, j:17):0,(d:1,e:4):2):3)")
        #           /----4--- /--1-a
        # ---------|          \--2-b
        #          |          /----0--- /-1---cparent---2---c
        #           \---3----|          \--17-j
        #                     \----2--- /--1--d
        #                               \--4--e
        # note c,j is len 0 node

        true_dists = {
            ("a", "b"): 3.0,
            ("a", "c"): 11.0,
            ("a", "d"): 11.0,
            ("a", "e"): 14.0,
            ("a", "j"): 25.0,
            ("b", "a"): 3.0,
            ("b", "c"): 12.0,
            ("b", "d"): 12.0,
            ("b", "e"): 15.0,
            ("b", "j"): 26.0,
            ("c", "a"): 11.0,
            ("c", "b"): 12.0,
            ("c", "d"): 6.0,
            ("c", "e"): 9.0,
            ("c", "j"): 20.0,
            ("d", "a"): 11.0,
            ("d", "b"): 12.0,
            ("d", "c"): 6.0,
            ("d", "e"): 5.0,
            ("d", "j"): 20.0,
            ("e", "a"): 14.0,
            ("e", "b"): 15.0,
            ("e", "c"): 9.0,
            ("e", "d"): 5.0,
            ("e", "j"): 23.0,
            ("j", "a"): 25.0,
            ("j", "b"): 26.0,
            ("j", "c"): 20.0,
            ("j", "d"): 20.0,
            ("j", "e"): 23.0,
        }

        true_root_dists = {"a": 5, "b": 6, "c": 6, "j": 20, "d": 6, "e": 9}

        t1_dists = t1.get_distances()
        subtree = t1.get_sub_tree(set(["d", "e", "c"]))
        sub_dists = subtree.get_distances()
        true_sub_root_dists = {"c": 3, "d": 3, "e": 6}

        sub_sameroot = t1.get_sub_tree(set(["d", "e", "c"]), keep_root=True)
        sub_sameroot_dists = sub_sameroot.get_distances()

        sub_sameroot2 = t1.get_sub_tree(set(["j", "c"]), keep_root=True)
        sub_sameroot_dists2 = sub_sameroot2.get_distances()

        # tip to tip dists should be the same
        for tip_pair in list(sub_dists.keys()):
            self.assertEqual(sub_dists[tip_pair], true_dists[tip_pair])
        for tip_pair in list(t1_dists.keys()):
            self.assertEqual(t1_dists[tip_pair], true_dists[tip_pair])
        for tip_pair in list(sub_sameroot_dists.keys()):
            self.assertEqual(sub_sameroot_dists[tip_pair], true_dists[tip_pair])
        for tip_pair in list(sub_sameroot_dists2.keys()):
            self.assertEqual(sub_sameroot_dists2[tip_pair], true_dists[tip_pair])

        # sameroot should have longer root to tip dists
        for tip in t1.tips():
            assert_allclose(t1.distance(tip), true_root_dists[tip.name])
        for tip in subtree.tips():
            assert_allclose(subtree.distance(tip), true_sub_root_dists[tip.name])
        for tip in sub_sameroot.tips():
            assert_allclose(sub_sameroot.distance(tip), true_root_dists[tip.name])
        for tip in sub_sameroot2.tips():
            assert_allclose(sub_sameroot2.distance(tip), true_root_dists[tip.name])

    def test_getsubtree_5(self):
        """get sub tree correctly uses tips only if specified"""
        names = [
            "Canis familiaris",
            "Mus musculus",
            "Homo sapiens",
            "Ornithorhynchus anatinus",
        ]
        treestring = (
            "((Homo sapiens:2,(Mus_musculus_129S1SvImJ:1,"
            "(Mus_musculus:0.1,Mus_musculus_LPJ:0.2):0.3)Mus_musculus:0.3)"
            ",Canis_familiaris,Ornithorhynchus_anatinus)"
        )
        tree = make_tree(treestring=treestring, underscore_unmunge=True)
        # change internal node names to eliminate ending digits
        for edge in tree.postorder():
            if edge.istip():
                continue
            if "Mus" in edge.name:
                edge.name = "Mus musculus"

        sub1 = tree.get_sub_tree(names, tipsonly=False)
        self.assertTrue(tree.same_topology(sub1))
        expect = make_tree(
            treestring="(Homo_sapiens,Mus_musculus,"
            "(Canis_familiaris,Ornithorhynchus_anatinus))",
            underscore_unmunge=True,
        )
        sub2 = tree.get_sub_tree(names, tipsonly=True)
        self.assertTrue(expect.same_topology(sub2))

    def test_getsubtree_6(self):
        """get sub tree handles non-scalar params"""
        names = ["A", "B", "C"]
        treestring = "((E:2,(F:1,(C:0.1,D:0.2):0.3):0.3):0.3,A:0.2,B:0.2)"
        tree = make_tree(treestring=treestring)
        # change internal node names to eliminate ending digits
        vals = dict(A=42, B=43, C=44)
        for edge in tree.postorder():
            if edge.is_root():
                continue
            edge.params["non-scalar"] = {edge.name: vals.get(edge.name, 99)}

        sub1 = tree.get_sub_tree(names, tipsonly=False)
        self.assertTrue(set(tree.get_tip_names()), set("ABC"))
        # the edge value for "non-scalar" should be same as original tree
        for name in vals:
            node = sub1.get_node_matching_name(name)
            self.assertTrue(node.params["non-scalar"], {name: vals[node.name]})

    def test_load_tree_from_json(self):
        """tests loading a Tree from json file"""
        with TemporaryDirectory(dir=".") as dirname:
            json_path = os.path.join(dirname, "tree.json")
            self.tree.write(json_path)
            got = load_tree(json_path)
            self.assertIsInstance(got, PhyloNode)
            self.assertEqual(
                got.get_newick(),
                # Since json file has rich dict as its content, parameter with_node_names is set True here.
                self.tree.get_newick(with_node_names=True),
            )
            self.assertEqual(got.get_node_names(), self.tree.get_node_names())
            # now try using non json suffix
            json_path = os.path.join(dirname, "tree.txt")
            self.tree.write(json_path, format="json")
            got = load_tree(json_path, format="json")
            self.assertIsInstance(got, PhyloNode)

    def test_load_tree(self):
        """tests loading a newick formatted Tree"""
        with TemporaryDirectory(dir=".") as dirname:
            tree_path = os.path.join(dirname, "tree.tree")
            self.tree.write(tree_path)
            got = load_tree(tree_path)
            self.assertIsInstance(got, PhyloNode)
            self.assertEqual(
                got.get_newick(),
                self.tree.get_newick(),
            )
            self.assertEqual(got.get_node_names(), self.tree.get_node_names())
            # now try specifying path as pathlib.Path
            tree_path = pathlib.Path(tree_path)
            got = load_tree(tree_path)
            self.assertIsInstance(got, PhyloNode)

    def test_ascii(self):
        self.tree.ascii_art()
        # unlabeled internal node
        tr = DndParser("(B:0.2,(C:0.3,D:0.4):0.6)F;")
        obs = tr.ascii_art(show_internal=True, compact=False)
        exp = "          /-B\n-F-------|\n         |          /-C\n          \\--------|\n                    \\-D"
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=True, compact=True)
        exp = "-F------- /-B\n          \\-------- /-C\n                    \\-D"
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=False, compact=False)
        exp = "          /-B\n---------|\n         |          /-C\n          \\--------|\n                    \\-D"
        self.assertEqual(obs, exp)


# the following class repeats the above tests but using a big tree and big
# data-set


class BigTreeSingleTests(TestTree):
    """using the big-tree for single-tree tests"""

    def setUp(self):
        self.name = "big tree - "
        self.otu_names = [
            "Horse",
            "TombBat",
            "Rhino",
            "Pig",
            "AsianElep",
            "SpermWhal",
            "Cat",
            "Gorilla",
            "Orangutan",
            "bandicoot",
            "Hedgehog",
            "Sloth",
            "HairyArma",
            "Manatee",
            "GoldenMol",
            "Pangolin",
        ]
        self.otu_names.sort()
        self.newick = "((((((((FlyingFox,DogFaced),((FreeTaile,LittleBro),(TombBat,RoundEare))),(FalseVamp,LeafNose)),(((Horse,Rhino),(Pangolin,(Cat,Dog))),(Llama,(Pig,(Cow,(Hippo,(SpermWhal,HumpbackW))))))),(Mole,Hedgehog)),(TreeShrew,(FlyingLem,((Jackrabbit,(FlyingSqu,(OldWorld,(Mouse,Rat)))),(Galago,(HowlerMon,(Rhesus,(Orangutan,(Gorilla,(Human,Chimpanzee)))))))))),(((NineBande,HairyArma),(Anteater,Sloth)),(((Dugong,Manatee),((AfricanEl,AsianElep),(RockHyrax,TreeHyrax))),(Aardvark,((GoldenMol,(Madagascar,Tenrec)),(LesserEle,GiantElep)))))),(caenolest,(phascogale,(wombat,bandicoot))));"
        self.newick_reduced = "(((((TombBat,(((Horse,Rhino),(Pangolin,Cat)),(Pig,SpermWhal))),Hedgehog),(Orangutan,Gorilla)),((HairyArma,Sloth),((Manatee,AsianElep),GoldenMol))),bandicoot);"
        self.tree = make_tree(treestring=self.newick)

    def test_get_edge_names(self):
        """testing (well, exercising at least), getedgenames"""
        # Fell over on small tree because "stem descended from root
        # joiner was a tip"
        a, b = self.otu_names[:2]
        self.tree.get_edge_names(a, b, True, False)

    def test_get_tip_names(self):
        """testing (well, exercising at least), get_tip_names"""
        a, b = self.otu_names[:2]
        tips = self.tree.get_tip_names()
        self.assertEqual(len(tips), 55)


# run if called from command line
if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Unit tests for tree parsers.
"""
from unittest import TestCase, main

from cogent3.core.tree import PhyloNode
from cogent3.parse.tree import DndParser, DndTokenizer, RecordError


# from cogent3.parse.newick import parse_string, TreeParseError as RecordError
# def DndParser(data, NodeClass=PhyloNode, unescape_name=True):
#    if not unescape_name:
#        raise NotImplementedError
#    def constructor(children, name, attribs):
#        return NodeClass(children = list(children or []), name=name, params=attribs)
#    return parse_string(data, constructor)

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
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

minimal = "();"
no_names = "((,),(,));"
missing_tip_name = "((a,b),(c,));"

empty = "();"
single = "(abc:3);"
double = "(abc:3, def:4);"
onenest = "(abc:3, (def:4, ghi:5):6 );"
nodedata = "(abc:3, (def:4, ghi:5)jkl:6 );"


class DndTokenizerTests(TestCase):
    """Tests of the DndTokenizer factory function."""

    def test_gdata(self):
        """DndTokenizer should work as expected on real data"""
        exp = [
            "(",
            "(",
            "xyz",
            ":",
            "0.28124",
            ",",
            "(",
            "def",
            ":",
            "0.24498",
            ",",
            "mno",
            ":",
            "0.03627",
            ")",
            ":",
            "0.17710",
            ")",
            ":",
            "0.04870",
            ",",
            "abc",
            ":",
            "0.05925",
            ",",
            "(",
            "ghi",
            ":",
            "0.06914",
            ",",
            "jkl",
            ":",
            "0.13776",
            ")",
            ":",
            "0.09853",
            ")",
            ";",
        ]
        # split it up for debugging on an item-by-item basis
        obs = list(DndTokenizer(sample))
        self.assertEqual(len(obs), len(exp))
        for i, j in zip(obs, exp):
            self.assertEqual(i, j)
        # try it all in one go
        self.assertEqual(list(DndTokenizer(sample)), exp)

    def test_nonames(self):
        """DndTokenizer should work as expected on trees with no names"""
        exp = ["(", "(", ",", ")", ",", "(", ",", ")", ")", ";"]
        obs = list(DndTokenizer(no_names))
        self.assertEqual(obs, exp)

    def test_missing_tip_name(self):
        """DndTokenizer should work as expected on trees with a missing name"""
        exp = ["(", "(", "a", ",", "b", ")", ",", "(", "c", ",", ")", ")", ";"]
        obs = list(DndTokenizer(missing_tip_name))
        self.assertEqual(obs, exp)

    def test_minimal(self):
        """DndTokenizer should work as expected a minimal tree without names"""
        exp = ["(", ")", ";"]
        obs = list(DndTokenizer(minimal))
        self.assertEqual(obs, exp)


class DndParserTests(TestCase):
    """Tests of the DndParser factory function."""

    def test_nonames(self):
        """DndParser should produce the correct tree when there are no names"""
        obs = DndParser(no_names)
        exp = PhyloNode()
        exp.append(PhyloNode())
        exp.append(PhyloNode())
        exp.children[0].append(PhyloNode())
        exp.children[0].append(PhyloNode())
        exp.children[1].append(PhyloNode())
        exp.children[1].append(PhyloNode())
        self.assertEqual(str(obs), str(exp))

    def test_minimal(self):
        """DndParser should produce the correct minimal tree"""
        obs = DndParser(minimal)
        exp = PhyloNode()
        exp.append(PhyloNode())
        self.assertEqual(str(obs), str(exp))

    def test_missing_tip_name(self):
        """DndParser should produce the correct tree when missing a name"""
        obs = DndParser(missing_tip_name)
        exp = PhyloNode()
        exp.append(PhyloNode())
        exp.append(PhyloNode())
        exp.children[0].append(PhyloNode(name="a"))
        exp.children[0].append(PhyloNode(name="b"))
        exp.children[1].append(PhyloNode(name="c"))
        exp.children[1].append(PhyloNode())
        self.assertEqual(str(obs), str(exp))

    def test_gsingle(self):
        """DndParser should produce a single-child PhyloNode on minimal data"""
        t = DndParser(single)
        self.assertEqual(len(t), 1)
        child = t[0]
        self.assertEqual(child.name, "abc")
        self.assertEqual(child.length, 3)
        self.assertEqual(str(t), "(abc:3.0);")

    def test_gdouble(self):
        """DndParser should produce a double-child PhyloNode from data"""
        t = DndParser(double)
        self.assertEqual(len(t), 2)
        self.assertEqual(str(t), "(abc:3.0,def:4.0);")

    def test_gonenest(self):
        """DndParser should work correctly with nested data"""
        t = DndParser(onenest)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  # first child is terminal
        self.assertEqual(len(t[1]), 2)  # second child has two children
        self.assertEqual(str(t), "(abc:3.0,(def:4.0,ghi:5.0):6.0);")

    def test_gnodedata(self):
        """DndParser should assign name to internal nodes correctly"""
        t = DndParser(nodedata)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  # first child is terminal
        self.assertEqual(len(t[1]), 2)  # second child has two children
        self.assertEqual(str(t), "(abc:3.0,(def:4.0,ghi:5.0)jkl:6.0);")
        info_dict = {}
        for node in t.traverse():
            info_dict[node.name] = node.length
        self.assertEqual(info_dict["abc"], 3.0)
        self.assertEqual(info_dict["def"], 4.0)
        self.assertEqual(info_dict["ghi"], 5.0)
        self.assertEqual(info_dict["jkl"], 6.0)

    def test_data(self):
        """DndParser should work as expected on real data"""
        t = DndParser(sample)
        self.assertEqual(
            str(t),
            "((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);",
        )
        tdata = DndParser(node_data_sample, unescape_name=True)
        self.assertEqual(
            str(tdata),
            "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)C:0.09853);",
        )

    def test_gbad(self):
        """DndParser should fail if parens unbalanced"""
        left = "((abc:3)"
        right = "(abc:3))"
        self.assertRaises(RecordError, DndParser, left)
        self.assertRaises(RecordError, DndParser, right)

    def test_DndParser(self):
        """DndParser tests"""
        t_str = "(A_a,(B:1.0,C),'D_e':0.5)E;"
        tree_unesc = DndParser(t_str, PhyloNode, unescape_name=True)
        tree_esc = DndParser(t_str, PhyloNode, unescape_name=False)

        self.assertEqual(tree_unesc.name, "E")
        self.assertEqual(tree_unesc.children[0].name, "A a")
        self.assertEqual(tree_unesc.children[1].children[0].name, "B")
        self.assertEqual(tree_unesc.children[1].children[0].length, 1.0)
        self.assertEqual(tree_unesc.children[1].children[1].name, "C")
        self.assertEqual(tree_unesc.children[2].name, "D_e")
        self.assertEqual(tree_unesc.children[2].length, 0.5)

        self.assertEqual(tree_esc.name, "E")
        self.assertEqual(tree_esc.children[0].name, "A_a")
        self.assertEqual(tree_esc.children[1].children[0].name, "B")
        self.assertEqual(tree_esc.children[1].children[0].length, 1.0)
        self.assertEqual(tree_esc.children[1].children[1].name, "C")
        self.assertEqual(tree_esc.children[2].name, "'D_e'")
        self.assertEqual(tree_esc.children[2].length, 0.5)

        reload_test = tree_esc.get_newick(with_distances=True, escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        self.assertEqual(
            obs.get_newick(with_distances=True),
            tree_esc.get_newick(with_distances=True),
        )
        reload_test = tree_unesc.get_newick(with_distances=True, escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        self.assertEqual(
            obs.get_newick(with_distances=True),
            tree_unesc.get_newick(with_distances=True),
        )


class PhyloNodeTests(TestCase):
    """Check that PhyloNode works the way I think"""

    def test_gops(self):
        """Basic PhyloNode operations should work as expected"""
        p = PhyloNode()
        self.assertEqual(str(p), ";")
        p.name = "abc"
        self.assertEqual(str(p), "abc;")
        p.length = 3
        self.assertEqual(str(p), "abc:3;")  # don't suppress branch from root
        q = PhyloNode()
        p.append(q)
        self.assertEqual(str(p), "()abc:3;")
        r = PhyloNode()
        q.append(r)
        self.assertEqual(str(p), "(())abc:3;")
        r.name = "xyz"
        self.assertEqual(str(p), "((xyz))abc:3;")
        q.length = 2
        self.assertEqual(str(p), "((xyz):2)abc:3;")


if __name__ == "__main__":
    main()

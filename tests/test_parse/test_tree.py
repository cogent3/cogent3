"""Unit tests for tree parsers."""

from unittest import TestCase

import pytest

from cogent3.core.tree import PhyloNode
from cogent3.parse.tree import DndParser, DndTokenizer, RecordError

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
        assert len(obs) == len(exp)
        for i, j in zip(obs, exp, strict=False):
            assert i == j
        # try it all in one go
        assert list(DndTokenizer(sample)) == exp

    def test_nonames(self):
        """DndTokenizer should work as expected on trees with no names"""
        exp = ["(", "(", ",", ")", ",", "(", ",", ")", ")", ";"]
        obs = list(DndTokenizer(no_names))
        assert obs == exp

    def test_missing_tip_name(self):
        """DndTokenizer should work as expected on trees with a missing name"""
        exp = ["(", "(", "a", ",", "b", ")", ",", "(", "c", ",", ")", ")", ";"]
        obs = list(DndTokenizer(missing_tip_name))
        assert obs == exp

    def test_minimal(self):
        """DndTokenizer should work as expected a minimal tree without names"""
        exp = ["(", ")", ";"]
        obs = list(DndTokenizer(minimal))
        assert obs == exp


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
        assert str(obs) == str(exp)

    def test_minimal(self):
        """DndParser should produce the correct minimal tree"""
        obs = DndParser(minimal)
        exp = PhyloNode()
        exp.append(PhyloNode())
        assert str(obs) == str(exp)

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
        assert str(obs) == str(exp)

    def test_gsingle(self):
        """DndParser should produce a single-child PhyloNode on minimal data"""
        t = DndParser(single)
        assert len(t) == 1
        child = t[0]
        assert child.name == "abc"
        assert child.length == 3
        assert str(t) == "(abc:3.0);"

    def test_gdouble(self):
        """DndParser should produce a double-child PhyloNode from data"""
        t = DndParser(double)
        assert len(t) == 2
        assert str(t) == "(abc:3.0,def:4.0);"

    def test_gonenest(self):
        """DndParser should work correctly with nested data"""
        t = DndParser(onenest)
        assert len(t) == 2
        assert len(t[0]) == 0  # first child is terminal
        assert len(t[1]) == 2  # second child has two children
        assert str(t) == "(abc:3.0,(def:4.0,ghi:5.0):6.0);"

    def test_gnodedata(self):
        """DndParser should assign name to internal nodes correctly"""
        t = DndParser(nodedata)
        assert len(t) == 2
        assert len(t[0]) == 0  # first child is terminal
        assert len(t[1]) == 2  # second child has two children
        assert str(t) == "(abc:3.0,(def:4.0,ghi:5.0):6.0);"
        info_dict = {node.name: node.length for node in t.preorder()}
        assert info_dict["abc"] == 3.0
        assert info_dict["def"] == 4.0
        assert info_dict["ghi"] == 5.0
        assert info_dict["jkl"] == 6.0

    def test_data(self):
        """DndParser should work as expected on real data"""
        t = DndParser(sample)
        assert (
            str(t)
            == "((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);"
        )
        tdata = DndParser(node_data_sample, unescape_name=True)
        assert (
            tdata.get_newick(with_distances=True, with_node_names=True)
            == "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)C:0.09853);"
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

        assert tree_unesc.name == "E"
        assert tree_unesc.children[0].name == "A a"
        assert tree_unesc.children[1].children[0].name == "B"
        assert tree_unesc.children[1].children[0].length == 1.0
        assert tree_unesc.children[1].children[1].name == "C"
        assert tree_unesc.children[2].name == "D_e"
        assert tree_unesc.children[2].length == 0.5

        assert tree_esc.name == "E"
        assert tree_esc.children[0].name == "A_a"
        assert tree_esc.children[1].children[0].name == "B"
        assert tree_esc.children[1].children[0].length == 1.0
        assert tree_esc.children[1].children[1].name == "C"
        assert tree_esc.children[2].name == "'D_e'"
        assert tree_esc.children[2].length == 0.5

        reload_test = tree_esc.get_newick(with_distances=True, escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        assert obs.get_newick(with_distances=True) == tree_esc.get_newick(
            with_distances=True,
        )
        reload_test = tree_unesc.get_newick(with_distances=True, escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        assert obs.get_newick(with_distances=True) == tree_unesc.get_newick(
            with_distances=True,
        )


class PhyloNodeTests(TestCase):
    """Check that PhyloNode works the way I think"""

    def test_gops(self):
        """Basic PhyloNode operations should work as expected"""
        p = PhyloNode()
        assert p.get_newick(with_node_names=True) == ";"
        p.name = "abc"
        assert p.get_newick(with_node_names=True, with_root_name=True) == "abc;"
        p.length = 3
        assert (
            p.get_newick(with_node_names=True, with_distances=True, with_root_name=True)
            == "abc:3;"
        )  # don't suppress branch from root
        q = PhyloNode()
        p.append(q)
        assert (
            p.get_newick(with_node_names=True, with_distances=True, with_root_name=True)
            == "()abc:3;"
        )
        r = PhyloNode()
        q.append(r)
        assert (
            p.get_newick(with_node_names=True, with_distances=True, with_root_name=True)
            == "(())abc:3;"
        )
        r.name = "xyz"
        assert (
            p.get_newick(with_node_names=True, with_distances=True, with_root_name=True)
            == "((xyz))abc:3;"
        )
        q.length = 2
        assert (
            p.get_newick(with_node_names=True, with_distances=True, with_root_name=True)
            == "((xyz):2)abc:3;"
        )


def test_make_tree_simple():
    from cogent3 import make_tree

    tree = make_tree(treestring="(a:0.4,b:0.5,c:0.6)")
    assert tree.get_node_matching_name("a").length == 0.4
    assert tree.get_node_matching_name("b").length == 0.5
    assert tree.get_node_matching_name("c").length == 0.6


def test_stores_nhx_data(DATA_DIR):
    """capture extended newick data into <node>.params['other']"""
    with open(DATA_DIR / "nhx.tree") as path:
        got = DndParser(path)
        rodent = got.get_connecting_node("Mouse", "Rat")
        assert set(rodent.params["other"]) == {
            '&sCF="89.13"',
            'sCF/sDF1/sDF2="89.13/5.49/5.38"',
            'sCF_N="107.21"',
            'sCF_N/sDF1_N/sDF2_N="107.21/6.63/6.49"',
            'sDF1="5.49"',
            'sDF1_N="6.63"',
            'sDF2="5.38"',
            'sDF2_N="6.49"',
            'sN="120.33"',
        }


def test_subtree_nhx_data(DATA_DIR):
    """capture extended newick data into <node>.params['other']"""
    with open(DATA_DIR / "nhx.tree") as path:
        tree = DndParser(path)

    subtree = tree.get_sub_tree(["Seal", "Mouse", "Rat", "Platypus", "Opossum"])
    rodent = subtree.get_connecting_node("Mouse", "Rat")
    assert set(rodent.params["other"]) == {
        '&sCF="89.13"',
        'sCF/sDF1/sDF2="89.13/5.49/5.38"',
        'sCF_N="107.21"',
        'sCF_N/sDF1_N/sDF2_N="107.21/6.63/6.49"',
        'sDF1="5.49"',
        'sDF1_N="6.63"',
        'sDF2="5.38"',
        'sDF2_N="6.49"',
        'sN="120.33"',
    }


def test_make_tree_with_nhx(DATA_DIR):
    """capture extended newick data into <node>.params['other']"""
    from cogent3 import make_tree

    nwk = (
        '(t1:0.1,(t2:0.2,t3:0.3):0.4[&sCF="81.73"],'
        '((t4:0.5,t5:0.6):0.7[&sCF="50.68",sDF1="24.42"],'
        't6:0.8):0.9[&sCF="77.1",sDF1="11.53"]);'
    )
    got = make_tree(nwk)
    vert = got.get_connecting_node("t2", "t3")
    assert set(vert.params["other"]) == {'&sCF="81.73"'}
    vert = got.get_connecting_node("t4", "t5")
    assert set(vert.params["other"]) == {'&sCF="50.68"', 'sDF1="24.42"'}
    vert = got.get_connecting_node("t5", "t6")
    assert set(vert.params["other"]) == {'&sCF="77.1"', 'sDF1="11.53"'}


@pytest.mark.parametrize(
    "treestring", ["(a,b)ab;", "(a,b)ab", "(a,b,(c,(d,e)))ab", "(a:3,b)ab"]
)
def test_parse_named_root_nodes(treestring):
    tree = DndParser(treestring, PhyloNode)
    assert tree.name == "ab"
    assert tree.parent is None

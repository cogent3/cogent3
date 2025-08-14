"""Tests of classes for dealing with trees and phylogeny."""

import json
import pathlib
import random
from collections.abc import Callable, Iterable
from copy import deepcopy
from tempfile import TemporaryDirectory
from typing import Any

import pytest
from numpy import array
from numpy.testing import assert_allclose, assert_equal

from cogent3 import get_dataset, load_tree, make_tree, open_
from cogent3._version import __version__
from cogent3.core.tree import PhyloNode, TreeBuilder, TreeError, split_name_and_support
from cogent3.parse.tree import DndParser
from cogent3.util.misc import get_object_provenance


def test_make_tree():
    """make_tree should load a tree from a file or a string"""
    # NOTE: This method now sits in cogent3.__init__

    t_str = "(a_a:10,(b_b:2,c_c:4):5);"
    # Tree quotes these labels because they have underscores in them.
    result_str = "('a_a':10.0,('b_b':2.0,'c_c':4.0):5.0);"
    t = make_tree(treestring=t_str)
    names = [i.name for i in t.tips()]
    assert names == ["a_a", "b_b", "c_c"]
    assert str(t) == result_str
    assert t.get_newick(with_distances=True) == result_str
    t_str = "(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);"
    # NOTE: Tree silently converts spaces to underscores (only for output),
    # presumably for Newick compatibility.
    result_str = "(a_a:10.0,(b_b:2.0,c_c:4.0):5.0);"
    t = make_tree(treestring=t_str, underscore_unmunge=True)
    names = [i.name for i in t.tips()]
    assert names == ["a a", "b b", "c c"]
    assert str(t) == result_str
    assert t.get_newick(with_distances=True) == result_str
    # ensure tip names are converted to strings
    # when creating a tree from a list of integer tip names.
    t = make_tree(tip_names=[1, 2, 3])
    assert t.get_tip_names() == ["1", "2", "3"]


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


@pytest.fixture
def empty_node() -> PhyloNode:
    return PhyloNode("")


def test_init_empty(empty_node: PhyloNode):
    """Empty PhyloNode should init OK"""
    t = empty_node
    assert t.name == ""
    assert t.parent is None
    assert len(t) == 0


def test_init_full(empty_node: PhyloNode):
    """PhyloNode should init OK with parent, data, and children"""
    t = empty_node
    u = PhyloNode("abc", parent=t, children="xyz")
    assert u.name == "abc"
    assert u.parent is t
    assert u in t
    assert u[0].name == "x"
    assert u[1].name == "y"
    assert u[2].name == "z"
    assert len(u) == 3


@pytest.fixture
def child_node():
    return PhyloNode(name="b")


@pytest.fixture
def one_child(child_node: PhyloNode) -> PhyloNode:
    return PhyloNode(name="a", children=[child_node])


@pytest.fixture
def big_name() -> list[PhyloNode]:
    """PhyloNode with many children"""
    return list(map(PhyloNode, "0123456789"))


@pytest.fixture
def big_parent(big_name: list[PhyloNode]) -> PhyloNode:
    """PhyloNode with many children"""
    return PhyloNode(name="x", children=big_name)


def test_str(empty_node: PhyloNode, one_child: PhyloNode, big_parent: PhyloNode):
    """PhyloNode str should give Newick-style representation"""
    # note: name suppressed if None
    assert str(empty_node) == ";"
    assert one_child.get_newick(with_node_names=True, with_root_name=True) == "(b)a;"
    assert (
        big_parent.get_newick(with_node_names=True, with_root_name=True)
        == "(0,1,2,3,4,5,6,7,8,9)x;"
    )
    big_parent[-1].extend("abc")
    assert (
        big_parent.get_newick(with_node_names=True, with_root_name=True)
        == "(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;"
    )


def test_get_newick(empty_node: PhyloNode, one_child: PhyloNode, big_parent: PhyloNode):
    """Should return Newick-style representation"""
    assert empty_node.get_newick() == ";"
    assert one_child.get_newick(with_node_names=True, with_root_name=True) == "(b)a;"
    assert (
        big_parent.get_newick(with_node_names=True, with_root_name=True)
        == "(0,1,2,3,4,5,6,7,8,9)x;"
    )
    big_parent[-1].extend("abc")
    assert (
        big_parent.get_newick(with_node_names=True, with_root_name=True)
        == "(0,1,2,3,4,5,6,7,8,(a,b,c)9)x;"
    )


def test_to_dict():
    """tree produces dict"""
    tr = make_tree(treestring="(a,b,(c:0.3,d)e1/10)")
    for node in tr.preorder():
        _ = node.length
    got = tr.to_rich_dict()
    attrs = {}
    no_length_support = {"length": None, "support": None}
    expect: dict[str, Any] = {
        "newick": "(a,b,(c,d)e1)root",
        "edge_attributes": {
            "a": attrs,
            "b": attrs,
            "c": attrs,
            "d": attrs,
            "e1": attrs,
            "root": attrs,
        },
        "length_and_support": {
            "a": no_length_support,
            "b": no_length_support,
            "c": {"length": 0.3, "support": None},
            "d": no_length_support,
            "e1": {"length": None, "support": 10},
            "root": no_length_support,
        },
        "type": get_object_provenance(tr),
        "version": __version__,
    }
    assert got == expect


def test_to_json():
    """tree produces json string that round trips correctly"""
    tr = make_tree(treestring="(a,b,(c,d)e1)")
    got = json.loads(tr.to_json())
    expect = tr.to_rich_dict()
    assert got == expect

    tr = make_tree(treestring="(a:1,b:1,(c:1,d:1)e1:1)")
    got = json.loads(tr.to_json())
    expect = tr.to_rich_dict()
    assert got == expect


def test_write_to_json(DATA_DIR: pathlib.Path, tmp_path: pathlib.Path):
    tree = load_tree(filename=DATA_DIR / "brca1_5.tree")
    json_path = tmp_path / "brca1_5.json"
    tree.write(json_path)
    with open_(json_path) as fn:
        got = json.loads(fn.read())
        assert got["type"] == get_object_provenance(PhyloNode)
        assert (
            tree.get_newick(semicolon=False, with_node_names=True, with_root_name=True)
            == got["newick"]
        )
        assert set(tree.get_node_names()) == got["edge_attributes"].keys()


def test_write_to_txt(DATA_DIR: pathlib.Path, tmp_path: pathlib.Path):
    """write a tree to newick"""
    tree = load_tree(filename=DATA_DIR / "brca1_5.tree")
    out_path = tmp_path / "brca1_5.txt"
    tree.write(out_path)
    with open_(out_path) as fn:
        got = fn.read()
        assert got.count("(") == got.count(")") == 3


def test_multifurcating():
    """Coerces nodes to have <= n children"""
    t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
    t = DndParser(t_str)

    # can't break up easily... sorry 80char
    exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
    obs = t.multifurcating(2)
    assert (
        obs.get_newick(with_distances=True, with_node_names=True, with_root_name=True)
        == exp_str
    )
    assert t.get_newick(with_distances=True) != obs.get_newick(with_distances=True)

    obs = t.multifurcating(2, 0.5)
    exp_str = "((a:1.0,(b:2.0,c:3.0):0.5)d:4.0,((e:5.0,(f:6.0,g:7.0):0.5)h:8.0,(i:9.0,(j:10.0,k:11.0):0.5)l:12.0):0.5)m:14.0;"
    assert (
        obs.get_newick(with_distances=True, with_node_names=True, with_root_name=True)
        == exp_str
    )

    t_str = "((a,b,c)d,(e,f,g)h,(i,j,k)l)m;"
    exp_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
    t = DndParser(t_str, constructor=PhyloNode)
    obs = t.multifurcating(2)
    assert (
        obs.get_newick(with_distances=True, with_node_names=True, with_root_name=True)
        == exp_str
    )
    obs = t.multifurcating(2, eps=10)  # no effect on PhyloNode type
    assert (
        obs.get_newick(with_distances=True, with_node_names=True, with_root_name=True)
        == exp_str
    )

    with pytest.raises(TreeError):
        # PhyloNode does not support multifurcating with n=1
        t.multifurcating(1)


def test_multifurcating_nameunnamed():
    """Coerces nodes to have <= n children"""
    t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
    t = DndParser(t_str)

    obs = t.multifurcating(2, name_unnamed=True)

    c0, c1 = obs.children
    assert c0.children[1].name.startswith("AUTO")
    assert c1.name.startswith("AUTO")
    assert c1.children[0].children[1].name.startswith("AUTO")
    assert c1.children[1].children[1].name.startswith("AUTO")
    assert len(c0.children[1].name) == 22
    assert len(c1.name) == 22
    assert len(c1.children[0].children[1].name) == 22
    assert len(c1.children[1].children[1].name) == 22
    names = [n.name for n in t.nontips()]
    assert len(names) == len(set(names))


def test_bifurcating():
    """Coerces nodes to have <= 2 children"""
    t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
    t = DndParser(t_str)

    # can't break up easily... sorry 80char
    got = t.bifurcating()
    num_children = {len(c.children) for c in got.nontips()}
    assert num_children == {2}


@pytest.fixture
def tree_nodes() -> dict[str, PhyloNode]:
    nodes = {x: PhyloNode(x) for x in "abcdefgh"}
    nodes["a"].append(nodes["b"])
    nodes["b"].append(nodes["c"])
    nodes["c"].append(nodes["d"])
    nodes["c"].append(nodes["e"])
    nodes["c"].append(nodes["f"])
    nodes["f"].append(nodes["g"])
    nodes["a"].append(nodes["h"])
    return nodes


@pytest.fixture
def tree_root(tree_nodes: dict[str, PhyloNode]) -> PhyloNode:
    return tree_nodes["a"]


def test_eq1(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode comparison should compare using id"""
    nodes = tree_nodes
    assert nodes["a"] == nodes["a"]
    assert nodes["b"] != nodes["a"]
    assert nodes["a"] != nodes["b"]


@pytest.fixture
def comparisons() -> tuple[PhyloNode, PhyloNode, PhyloNode]:
    return PhyloNode("a"), PhyloNode("a"), PhyloNode("b")


def test_eq2(comparisons: tuple[PhyloNode, PhyloNode, PhyloNode]):
    """PhyloNode should compare equal if same id"""
    t, u, v = comparisons
    z = t
    assert t == z
    assert t is not u
    assert t != u
    assert t != v

    f = PhyloNode("a")
    g = PhyloNode("aa")
    assert f != g
    f.name += "a"
    assert f != g

    # however, two PhyloNodes that have no name should not compare equal
    f = PhyloNode("")
    g = PhyloNode("")
    assert f != g

    f = PhyloNode(name="foo")
    g = f.copy()
    assert f != g


def test_compare_name(tree_nodes: dict[str, PhyloNode]):
    """Compare names between PhyloNodes"""
    nodes = tree_nodes
    assert nodes["a"].compare_name(nodes["a"])
    assert not nodes["a"].compare_name(nodes["b"])
    assert not nodes["b"].compare_name(nodes["a"])


@pytest.fixture
def treestring_1() -> str:
    return "((H,G),(R,M));"


@pytest.fixture
def tree_1(treestring_1: str) -> PhyloNode:
    return DndParser(treestring_1, PhyloNode)


@pytest.fixture
def treestring_2() -> str:
    return "(((H,G),R),M);"


@pytest.fixture
def tree_2(treestring_2: str) -> PhyloNode:
    return DndParser(treestring_2, PhyloNode)


@pytest.fixture
def treestring_3() -> str:
    return "(((H,G),(O,R)),X);"


@pytest.fixture
def tree_3(treestring_3: str) -> PhyloNode:
    return DndParser(treestring_3, PhyloNode)


def test_compare_by_names(tree_1: PhyloNode, tree_2: PhyloNode, tree_3: PhyloNode):
    """Compare names between trees"""
    assert tree_1.compare_by_names(tree_2)
    assert tree_1.compare_by_names(tree_1)
    assert not tree_1.compare_by_names(tree_3)


def test_ne(comparisons: tuple[PhyloNode, PhyloNode, PhyloNode]):
    """PhyloNode should compare ne by id or data"""
    t, u, _ = comparisons
    z = t
    assert t == z
    assert t != u

    f = PhyloNode(name="foo")
    g = f.copy()
    assert f != g


def test_append(one_child: PhyloNode):
    """PhyloNode append should add item to end of self"""
    one_child.append(PhyloNode("c"))
    assert len(one_child) == 2
    assert one_child[-1].name == "c"
    one_child.append("6")
    assert len(one_child) == 3
    assert one_child[-1].name == "6"
    # check that refs are updated when moved from one tree to another
    empty = PhyloNode("")
    empty.append(one_child[-1])
    assert len(empty) == 1
    assert empty[0].name == "6"
    assert empty[0].parent == empty
    assert one_child[-1].name == "c"


def test_extend(empty_node: PhyloNode):
    """PhyloNode extend should add many items to end of self"""
    empty_node.extend("abcdefgh")
    data = "".join([i.name for i in empty_node])
    assert data == "abcdefgh"


def test_pop(big_parent: PhyloNode, big_name: list[PhyloNode]):
    """PhyloNode pop should remove and return child at specified index"""
    parent, nodes = big_parent, big_name
    assert len(parent) == 10
    last = parent.pop()
    assert last is nodes[-1]
    assert last.parent is None
    assert len(parent) == 9
    assert parent[-1] is nodes[-2]
    first = parent.pop(0)
    assert first is nodes[0]
    assert first.parent is None
    assert len(parent) == 8
    assert parent[0] is nodes[1]
    second_to_last = parent.pop(-2)
    assert second_to_last is nodes[-3]


def test_remove():
    """PhyloNode remove should remove first match by value, not id"""
    nodes = list(map(PhyloNode, "abc" * 3))
    parent = PhyloNode("", children=nodes)
    assert len(parent) == 9
    parent.remove("a")
    assert len(parent) == 8
    assert "".join([i.name for i in parent]) == "bcabcabc"
    new_node = PhyloNode("a")
    parent.remove(new_node)
    assert len(parent) == 7
    assert "".join([i.name for i in parent]) == "bcbcabc"


def test_insert(big_parent: PhyloNode):
    """PhyloNode insert should insert item at specified index"""
    parent = big_parent
    assert len(parent) == 10
    parent.insert(3, "5")
    assert len(parent) == 11
    assert parent[3].name == "5"
    assert parent[4].name == "3"
    parent.insert(-1, "123")
    assert len(parent) == 12
    assert parent[-1].name == "9"
    assert parent[-2].name == "123"


def test_getitem(tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]):
    """PhyloNode getitem should return item or slice"""
    r = tree_root
    n = tree_nodes
    assert r[0] is n["b"]
    items = n["c"][0:1]
    assert len(items) == 1
    assert items[0] is n["d"]
    items = n["c"][0:2]
    assert len(items) == 2
    assert items[0] is n["d"]
    assert items[1] is n["e"]
    items = n["c"][:]
    assert len(items) == 3
    assert items[0] is n["d"]
    assert items[-1] is n["f"]


def test_slice(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode slicing should return list, not PhyloNode"""
    nodes = tree_nodes
    c, d, e, f = nodes["c"], nodes["d"], nodes["e"], nodes["f"]
    assert c[:] is not c
    assert c[:] == [d, e, f]
    assert c[1:2] == [e]
    assert c[0:3:2] == [d, f]


def test_setitem(big_parent: PhyloNode, big_name: list[PhyloNode]):
    """PhyloNode setitem should set item or extended slice of nodes"""
    parent, nodes = big_parent, big_name
    t = PhyloNode("1")
    parent[0] = t
    assert parent[0] is t
    assert t.parent is parent
    assert nodes[0].parent is None

    u = PhyloNode("2")
    parent[-2] = u
    assert parent[8] is u
    assert u.parent is parent
    assert nodes[8].parent is None

    parent[1:6:2] = "xyz"
    for i in [1, 3, 5]:
        assert nodes[i].parent is None
    assert parent[1].name == "x"
    assert parent[3].name == "y"
    assert parent[5].name == "z"
    for i in parent:
        assert i.parent is parent


def test_setslice(big_parent: PhyloNode, big_name: list[PhyloNode]):
    """PhyloNode setslice should set old-style slice of nodes"""
    parent, nodes = big_parent, big_name
    assert len(parent) == 10
    parent[5:] = []
    assert len(parent) == 5
    for i in range(5, 10):
        assert nodes[i].parent is None
    parent[1:3] = "abcd"
    assert len(parent) == 7
    for i in parent:
        assert i.parent is parent
    data_list = [i.name for i in parent]
    assert data_list == list("0abcd34")
    parent[1:3] = parent[2:3]
    data_list = [i.name for i in parent]
    assert data_list == list("0bcd34")


def test_delitem(
    child_node: PhyloNode,
    one_child: PhyloNode,
    big_name: list[PhyloNode],
    big_parent: PhyloNode,
):
    """PhyloNode __delitem__ should delete item and set parent to None"""
    assert child_node.parent == one_child
    assert len(one_child) == 1
    del one_child[0]
    assert one_child.parent is None
    assert len(one_child) == 0

    nodes = big_name
    parent = big_parent
    assert len(parent) == 10
    for n in nodes:
        assert n.parent is parent
    del parent[-1]
    assert nodes[-1].parent is None
    assert len(parent) == 9
    del parent[1:6:2]
    assert len(parent) == 6
    for i, n in enumerate(nodes):
        if i in [0, 2, 4, 6, 7, 8]:
            assert n.parent is parent
        else:
            assert n.parent is None


def test_delslice(big_name: list[PhyloNode], big_parent: PhyloNode):
    """PhyloNode __delslice__ should delete items from start to end"""
    parent = big_parent
    nodes = big_name
    assert len(parent) == 10
    del parent[3:-2]
    assert len(parent) == 5
    for i, n in enumerate(nodes):
        if i in [3, 4, 5, 6, 7]:
            assert n.parent is None
        else:
            assert n.parent is parent


def test_len(tree_root: PhyloNode):
    """PhyloNode len should return number of children"""
    r = tree_root
    assert len(r) == 2


def test_copy():
    """PhyloNode.copy() should work on deep trees"""
    t = comb_tree(1024)  # should break recursion limit on regular copy
    t.name = "foo"
    t.params["XYZ"] = [3]
    t2 = t.copy()
    t3 = t.copy()
    t3.name = "bar"

    assert len(t.tips()) == 1024
    assert len(t2.tips()) == 1024
    assert len(t3.tips()) == 1024

    assert t is not t2
    assert t.name == t2.name
    assert t.name != t3.name

    assert t.params["XYZ"] == t2.params["XYZ"]
    assert t.params["XYZ"] is not t2.params["XYZ"]

    assert t.get_newick() == t2.get_newick()

    t_simple = PhyloNode("t")
    u_simple = PhyloNode("u")
    t_simple.append(u_simple)

    assert str(t_simple.copy()) == str(t_simple.copy())


def test_iter(tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]):
    """PhyloNode iter should iterate over children"""
    r = tree_root
    n = tree_nodes
    items = list(r)
    assert items[0] is n["b"]
    assert items[1] is n["h"]
    assert len(items) == 2


def test_deepcopy(tree_root: PhyloNode):
    """copy.deepcopy should work on PhyloNode"""
    t = PhyloNode('["t"]')
    u = PhyloNode('["u"]')
    t.append(u)

    c = deepcopy(u)
    assert c is not u
    assert c.name == u.name
    # note: name _is_ same object if it's immutable, e.g. a string.
    # deepcopy doesn't copy data for immutable objects.

    # need to check that we also copy arbitrary attributes
    t.params["XYZ"] = [3]
    c = deepcopy(t)
    assert c is not t
    assert c[0] is not u
    assert c[0].name == u.name
    assert c.params["XYZ"] == t.params["XYZ"]
    assert c.params["XYZ"] is not t.params["XYZ"]

    t = tree_root
    c = deepcopy(t)
    assert str(c) == str(t)


@pytest.fixture
def single_node() -> PhyloNode:
    return PhyloNode(name="a")


@pytest.fixture
def repeated_child() -> PhyloNode:
    return PhyloNode(name="x", children="aaa")


def test_Parent(
    single_node: PhyloNode,
    empty_node: PhyloNode,
    one_child: PhyloNode,
    repeated_child: PhyloNode,
):
    """PhyloNode parent should hold correct data and be mutable"""
    # check initial conditions
    assert single_node.parent is None
    # set parent and check parent/child relations
    single_node.parent = empty_node
    assert single_node.parent is empty_node
    assert empty_node[0] == single_node
    assert single_node in empty_node
    assert len(empty_node) == 1
    # reset parent and check parent/child relations
    single_node.parent = one_child
    assert single_node.parent is one_child
    assert single_node not in empty_node
    assert single_node is one_child[-1]

    # following is added to check that we don't screw up when there are
    # nodes with different ids that still compare equal
    for i in repeated_child:
        assert i.parent is repeated_child
    last = repeated_child[-1]
    last.parent = one_child
    assert len(repeated_child) == 2
    for i in repeated_child:
        assert i.parent is repeated_child
    assert last.parent is one_child


def test_index_in_parent():
    """PhyloNode index_in_parent should hold correct data"""
    first = PhyloNode("a")
    second = PhyloNode("b")
    third = PhyloNode("c")
    fourth = PhyloNode("0", children=[first, second, third])
    assert len(fourth) == 3
    assert first.index_in_parent() == 0
    assert second.index_in_parent() == 1
    assert third.index_in_parent() == 2
    del fourth[0]
    assert second.index_in_parent() == 0
    assert third.index_in_parent() == 1
    assert len(fourth) == 2
    assert first.parent is None


def test_is_tip(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode is_tip should return True if node is a tip"""
    tips = "degh"
    for n in list(tree_nodes.values()):
        if n.name in tips:
            assert n.is_tip() is True
        else:
            assert n.is_tip() is False


def test_isRoot(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode isRoot should return True if parent is None"""
    r = "a"
    for n in list(tree_nodes.values()):
        if n.name in r:
            assert n.is_root() is True
        else:
            assert n.is_root() is False


@pytest.fixture
def multi_child():
    return PhyloNode(name="a", children="bcd")


def test_levelorder():
    t = DndParser("(((A,B)C,(D,E)F,(G,H)I)J,(K,L)M)N;")
    exp = ["N", "J", "M", "C", "F", "I", "K", "L", "A", "B", "D", "E", "G", "H"]
    names = [n.name for n in t.levelorder()]
    assert names == exp


def test_ancestors(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode ancestors should provide list of ancestors, deepest first"""
    nodes = tree_nodes
    assert nodes["a"].ancestors() == []
    assert nodes["b"].ancestors() == [nodes["a"]]
    assert nodes["d"].ancestors() == nodes["f"].ancestors()
    assert nodes["g"].ancestors() == [
        nodes["f"],
        nodes["c"],
        nodes["b"],
        nodes["a"],
    ]


def test_newick_with_labelled_nodes():
    """return newick with internal nodes labelled"""
    treestrings = (
        "(a,b,(c,(d,e)))",
        "(a:0.1,b:0.2,(c:0.3,(d:0.4,e:0.5):0.6):0.7)",
        "(a,b,(c,(d,e)edge.0)edge.1)",
    )
    expect = (
        "(a,b,(c,(d,e)edge.0)edge.1)root",
        "(a:0.1,b:0.2,(c:0.3,(d:0.4,e:0.5)edge.0:0.6)edge.1:0.7)root",
        "(a,b,(c,(d,e)edge.0)edge.1)root",
    )
    for i, treestring in enumerate(treestrings):
        if i < 2:
            continue
        tree = make_tree(treestring=treestring)
        nwk = tree.get_newick(
            with_node_names=True,
            with_distances=True,
            semicolon=False,
            with_root_name=True,
        )
        assert nwk == expect[i]


def test_root(tree_nodes: dict[str, PhyloNode], tree_root: PhyloNode):
    """PhyloNode root() should find root of tree"""
    nodes, root = tree_nodes, tree_root
    for i in list(nodes.values()):
        assert i.get_root() is root


def test_children(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode children should allow getting/setting children"""
    nodes = tree_nodes
    for n in nodes:
        node = nodes[n]
        assert list(node) == node.children

    t = PhyloNode("", children="abc")
    assert len(t) == 3
    u, v = PhyloNode("u"), PhyloNode("v")

    # WARNING: If you set children directly, parent refs will _not_ update!
    t.children = [u, v]

    assert t[0] is u
    assert t[1] is v
    assert len(t) == 2


def test_siblings(
    tree_nodes: dict[str, PhyloNode],
    empty_node: PhyloNode,
    child_node: PhyloNode,
    one_child: PhyloNode,
):
    """PhyloNode siblings() should return all siblings, not self"""
    assert empty_node.siblings() == []
    assert child_node.siblings() == []
    assert one_child.siblings() == []

    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert g.siblings() == []
    assert f.siblings() == [d, e]
    assert e.siblings() == [d, f]
    assert d.siblings() == [e, f]
    assert c.siblings() == []
    assert b.siblings() == [h]
    assert h.siblings() == [b]
    assert a.siblings() == []


def test_tips(
    empty_node: PhyloNode,
    child_node: PhyloNode,
    one_child: PhyloNode,
    tree_nodes: dict[str, PhyloNode],
):
    """PhyloNode tips should return all terminal descendants"""
    assert empty_node.tips() == []
    assert child_node.tips() == []
    assert one_child.tips() == [child_node]

    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert g.tips() == []
    assert f.tips() == [g]
    assert e.tips() == []
    assert d.tips() == []
    assert c.tips() == [d, e, g]
    assert b.tips() == [d, e, g]
    assert h.tips() == []
    assert a.tips() == [d, e, g, h]


def test_itertips(tree_root: PhyloNode):
    """PhyloNode itertips should iterate over terminal descendants"""
    tree = tree_root
    assert [i.name for i in tree.iter_tips()] == list("degh")


def test_nontips(tree_root: PhyloNode):
    """PhyloNode nontips should return all non-terminal descendants"""
    tree = tree_root
    assert [i.name for i in tree.nontips()] == list("bcf")


def test_iterNonTips(tree_root: PhyloNode):
    """PhyloNode iter_nontips should iterate over non-terminal descendants"""
    tree = tree_root
    assert [i.name for i in tree.iter_nontips()] == list("bcf")


def test_tip_children(
    empty_node: PhyloNode,
    child_node: PhyloNode,
    one_child: PhyloNode,
    tree_nodes: dict[str, PhyloNode],
):
    """PhyloNode tip_children should return all terminal children"""
    assert empty_node.tip_children() == []
    assert child_node.tip_children() == []
    assert one_child.tip_children() == [child_node]

    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert g.tip_children() == []
    assert f.tip_children() == [g]
    assert e.tip_children() == []
    assert d.tip_children() == []
    assert c.tip_children() == [d, e]
    assert b.tip_children() == []
    assert h.tip_children() == []
    assert a.tip_children() == [h]


def test_non_tip_children(
    empty_node: PhyloNode,
    child_node: PhyloNode,
    one_child: PhyloNode,
    tree_nodes: dict[str, PhyloNode],
):
    """PhyloNode non_tip_children should return all non-terminal children"""
    assert empty_node.non_tip_children() == []
    assert child_node.non_tip_children() == []
    assert one_child.non_tip_children() == []

    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert g.non_tip_children() == []
    assert f.non_tip_children() == []
    assert e.non_tip_children() == []
    assert d.non_tip_children() == []
    assert c.non_tip_children() == [f]
    assert b.non_tip_children() == [c]
    assert h.non_tip_children() == []
    assert a.non_tip_children() == [b]


def test_remove_node(repeated_child: PhyloNode):
    """PhyloNode remove_node should delete node by id, not value"""
    parent = repeated_child
    children = list(repeated_child)
    assert len(parent) == 3
    assert parent.remove_node(children[1]) is True
    assert len(parent) == 2
    assert children[0].parent is parent
    assert children[1].parent is None
    assert children[2].parent is parent
    assert children[0].compare_name(children[1]) is True
    assert parent.remove_node(children[1]) is False
    assert len(parent) == 2
    assert parent.remove_node(children[0]) is True
    assert len(parent) == 1


def test_lowest_common_ancestor():
    """PhyloNode lowest_common_ancestor should return LCA for set of tips"""
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
    assert obs1 == exp1
    assert obs2 == exp2
    assert obs3 == exp3
    assert obs4 == exp4

    # verify multiple calls work
    t_mul = t1.copy()
    exp_1 = t_mul.get_node_matching_name("d")
    exp_2 = t_mul.get_node_matching_name("i")
    obs_1 = t_mul.lowest_common_ancestor(["b", "c"])
    obs_2 = t_mul.lowest_common_ancestor(["g", "h"])
    assert obs_1 == exp_1
    assert obs_2 == exp_2


def test_lowest_common_ancestor_invalid_tips():
    """fail if tips not present"""
    t = DndParser("((a,(b,c)d)e,f,(g,h)i)j;")
    # no tips present in tree should raise exception
    with pytest.raises(ValueError):
        t.lowest_common_ancestor(["m", "n"])

    # not all tips present in tree should raise exception
    with pytest.raises(ValueError):
        t.lowest_common_ancestor(["a", "n"])


def test_last_common_ancestor(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode last_common_ancestor should provide last common ancestor"""
    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert a.last_common_ancestor(a) == a
    assert a.last_common_ancestor(b) == a
    assert a.last_common_ancestor(g) == a
    assert a.last_common_ancestor(h) == a

    assert b.last_common_ancestor(g) == b
    assert b.last_common_ancestor(d) == b
    assert b.last_common_ancestor(a) == a
    assert b.last_common_ancestor(h) == a

    assert d.last_common_ancestor(f) == c
    assert d.last_common_ancestor(g) == c
    assert d.last_common_ancestor(a) == a
    assert d.last_common_ancestor(h) == a

    assert g.last_common_ancestor(g) == g
    assert g.last_common_ancestor(f) == f
    assert g.last_common_ancestor(e) == c
    assert g.last_common_ancestor(c) == c
    assert g.last_common_ancestor(b) == b
    assert g.last_common_ancestor(a) == a
    assert g.last_common_ancestor(h) == a

    t = PhyloNode("h")
    for i in [a, b, c, d, e, f, g, h]:
        with pytest.raises(TreeError, match="No common ancestor found"):
            i.last_common_ancestor(t)
        with pytest.raises(TreeError, match="No common ancestor found"):
            t.last_common_ancestor(i)

    PhyloNode("a", children=[t])


def test_separation(tree_nodes: dict[str, PhyloNode]):
    """PhyloNode separation should return correct number of edges"""
    nodes = tree_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert a.separation(a) == 0
    assert c.separation(c) == 0
    assert a.separation(b) == 1
    assert a.separation(h) == 1
    assert g.separation(h) == 5
    assert f.separation(d) == 2
    assert f.separation(c) == 1
    assert c.separation(f) == 1


def test_name_unnamed_nodes(tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]):
    """name_unnamed_nodes assigns an arbitrary value when name == None"""
    tree = tree_root
    tree_nodes["b"].name = "node2"
    tree_nodes["c"].name = ""
    tree_nodes["f"].name = ""
    tree_nodes["e"].name = "node3"
    tree.name_unnamed_nodes()
    assert tree_nodes["c"].name == "node1"
    assert tree_nodes["f"].name == "node4"


def test_make_tree_array(tree_root: PhyloNode):
    """make_tree_array maps nodes to the descendants in them"""
    tree = tree_root
    result, node_list = tree.make_tree_array()
    assert_equal(
        result,
        array([[1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 1, 0], [0, 0, 1, 0]]),
    )
    nodes = [node.name for node in node_list]
    assert nodes == ["a", "b", "c", "f"]
    # test if works with a dec_list supplied
    dec_list = ["d", "added", "e", "g", "h"]
    result2, node_list = tree.make_tree_array(dec_list)
    assert_equal(
        result2,
        array([[1, 0, 1, 1, 1], [1, 0, 1, 1, 0], [1, 0, 1, 1, 0], [0, 0, 0, 1, 0]]),
    )


def test_get_node_names():
    """get_node_names works correctly"""
    tree = make_tree(treestring="((a:3,(b:2,(c:1,d:1):1):1):2,(e:3,f:3):2);")
    names = tree.get_node_names(include_self=False, tips_only=False)
    assert tree.name not in names
    names = tree.get_node_names(include_self=True, tips_only=False)
    assert tree.name in names
    tree.get_node_matching_name("a")


def test_reassign_names(tree_root: PhyloNode):
    """reassign_names should rename node names based on dict mapping"""
    t = tree_root
    mapping = {x: str(i) for i, x in enumerate("abfg")}
    exp_names = ["0", "1", "2", "3", "c", "d", "e", "h"]
    t.reassign_names(mapping)
    obs_names = sorted(t.get_node_names())
    assert obs_names == exp_names


def test_reassign_names_specific_nodes(
    tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]
):
    """reassign_names should rename nodes based on dict mapping"""
    t = tree_root
    nodes = [tree_nodes["a"], tree_nodes["b"]]
    mapping = {x: str(i) for i, x in enumerate("abfg")}
    exp_names = ["0", "1", "c", "d", "e", "f", "g", "h"]
    t.reassign_names(mapping, nodes)
    obs_names = sorted(t.get_node_names())
    assert obs_names == exp_names


def test_get_nodes_dict(tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]):
    """get_nodes_dict returns a dict keyed by name, value is node"""
    t = tree_root
    nodes = tree_nodes
    assert t.get_nodes_dict() == nodes


def test_get_nodes_dict_nonunique_names(tree_root: PhyloNode):
    """get_nodes_dict raises if non unique names are in tree"""
    t = tree_root
    t.children[0].name = "same"
    t.children[0].children[0].name = "same"
    with pytest.raises(TreeError):
        t.get_nodes_dict()


def test_remove_deleted():
    """remove_deleted should remove all nodes where is_deleted tests true."""
    tree = DndParser(
        "((a:3,(b:2,(c:1,d:1):1):1):2,(e:3,f:3):2);",
        constructor=PhyloNode,
    )
    result_not_deleted = deepcopy(tree)
    tree.remove_deleted(lambda x: x.name in [])
    assert str(tree) == str(result_not_deleted)
    deleted = {"b", "d", "e", "f"}
    result_tree = DndParser("((a:3,((c:1):1):1):2);", constructor=PhyloNode)

    def is_deleted(x: PhyloNode):
        return x.name in deleted

    tree.remove_deleted(is_deleted)
    assert str(tree) == str(result_tree)


def test_prune():
    """prune should reconstruct correct topology of tree."""
    tree = DndParser("((a:3,((c:1):1):1):2);", constructor=PhyloNode)
    tree.prune()
    result_tree = DndParser("(a:3,c:3);", constructor=PhyloNode)
    assert str(tree) == str(result_tree)

    samename_bug = DndParser("((A,B)SAMENAME,((C,D)SAMENAME));")
    samename_bug.prune()
    exp_tree_str = "((A,B)SAMENAME,(C,D)SAMENAME);"
    assert samename_bug.get_newick(with_node_names=True) == exp_tree_str


def test_prune_2():
    tree = DndParser("(((d:1.0,e:4.0)x:2.0,c:3.0)y:3.0)root;")
    tree.prune()
    assert len(tree.children) == 2


def test_get_node_matching_name(tree_root: PhyloNode, tree_nodes: dict[str, PhyloNode]):
    """PhyloNode get_node_matching_name should return node that matches name"""
    nodes = tree_nodes
    root = tree_root
    assert root.get_node_matching_name("g") is nodes["g"]


def test_subset(tree_1: PhyloNode):
    """subset should return set of leaves that descends from node"""
    t = tree_1
    assert t.subset() == frozenset("HGRM")
    c = t.children[0]
    assert c.subset() == frozenset("HG")
    leaf = c.children[1]
    assert leaf.subset() == frozenset()


def test_subsets(tree_1: PhyloNode):
    """subsets should return all subsets descending from a set"""
    t = tree_1
    assert t.subsets() == frozenset([frozenset("HG"), frozenset("RM")])


def test_compare_by_subsets(
    tree_1: PhyloNode, tree_2: PhyloNode, tree_3: PhyloNode, tree_root: PhyloNode
):
    """compare_by_subsets should return the fraction of shared subsets"""
    result = tree_1.compare_by_subsets(tree_1)
    assert result == 0

    result = tree_2.compare_by_subsets(tree_2)
    assert result == 0

    result = tree_1.compare_by_subsets(tree_2)
    assert result == 0.5

    result = tree_1.compare_by_subsets(tree_3)
    assert result == 1 - 2.0 / 5

    result = tree_1.compare_by_subsets(tree_3, exclude_absent_taxa=True)
    assert result == 1 - 2.0 / 3

    result = tree_1.compare_by_subsets(tree_root, exclude_absent_taxa=True)
    assert result == 1

    result = tree_1.compare_by_subsets(tree_root)
    assert result == 1


def test_treenode_comparison_with_none_name(
    empty_node: PhyloNode, single_node: PhyloNode
):
    assert empty_node < single_node
    assert single_node > empty_node
    assert single_node > PhyloNode(name="")
    assert PhyloNode(name="") < single_node
    assert PhyloNode(name="test") > empty_node
    assert empty_node < PhyloNode(name="test")


@pytest.fixture
def phylo_nodes() -> dict[str, PhyloNode]:
    nodes = {x: PhyloNode(x) for x in "abcdefgh"}
    nodes["a"].append(nodes["b"])
    nodes["b"].append(nodes["c"])
    nodes["c"].append(nodes["d"])
    nodes["c"].append(nodes["e"])
    nodes["c"].append(nodes["f"])
    nodes["f"].append(nodes["g"])
    nodes["a"].append(nodes["h"])
    nodes["a"].length = None
    nodes["b"].length = 0
    nodes["c"].length = 3
    nodes["d"].length = 1
    nodes["e"].length = 4
    nodes["f"].length = 2
    nodes["g"].length = 3
    nodes["h"].length = 2
    return nodes


@pytest.fixture
def phylo_root(phylo_nodes: dict[str, PhyloNode]) -> PhyloNode:
    """Returns the root of the phylogenetic tree."""
    tree = phylo_nodes["a"]
    tree.prune()
    return tree


@pytest.fixture
def phylostring_1() -> str:
    """Returns a sample phylogenetic string."""
    return "((H:1,G:1):2,(R:0.5,M:0.7):3);"


@pytest.fixture
def phylo_1(phylostring_1: str) -> PhyloNode:
    """Returns a parsed phylogenetic tree from the string."""
    return DndParser(phylostring_1, PhyloNode)


@pytest.fixture
def phylostring_2() -> str:
    return "(((H:1,G:1,O:1):2,R:3):1,X:4);"


@pytest.fixture
def phylo_2(phylostring_2: str) -> PhyloNode:
    """Returns a parsed phylogenetic tree from the second string."""
    return DndParser(phylostring_2, PhyloNode)


def test_init():
    """Check PhyloNode constructor"""
    n = PhyloNode("foo", length=10)
    assert n.name == "foo"
    assert n.length == 10

    n = PhyloNode("bar")
    assert n.name == "bar"
    assert n.length is None

    n = PhyloNode("")
    assert n.name == ""
    assert n.length is None


def test_total_descending_branch_length(
    phylo_root: PhyloNode, phylo_nodes: dict[str, PhyloNode]
):
    """total_descending_branch_length returns total branchlength below self"""
    t = phylo_root
    exp = 15
    obs = t.total_descending_branch_length()
    assert obs == exp

    node_c = phylo_nodes["c"]
    exp = 10
    obs = node_c.total_descending_branch_length()
    assert obs == exp


def test_total_length():
    """total_length returns total branchlength, irrespective of edge"""
    t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2);"
    t = DndParser(t_str, constructor=PhyloNode)
    assert t.total_length() == 21
    node = t.get_node_matching_name("F")
    length = node.total_length()
    assert length == 21


def test_tips_within_distance():
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

    assert obs_at_dist_2 == exp_at_dist_2
    assert obs_at_dist_3 == exp_at_dist_3
    assert obs_at_dist_4 == exp_at_dist_4


def test_tips_within_distance_nodistances():
    """tips_within_distance returns tips that are within distance from self"""
    t_str = "(A,B,(C,D)E,(F,((G,H)I)J)K)L;"
    t = DndParser(t_str, constructor=PhyloNode)
    nodes = t.get_nodes_dict()
    e_node = nodes["E"]

    exp = sorted([n.name for n in t.tips()])
    obs = sorted([n.name for n in e_node.tips_within_distance(0)])
    assert obs == exp


def test_distance(phylo_nodes: dict[str, PhyloNode]):
    """PhyloNode distance should report correct distance between nodes"""
    nodes = phylo_nodes
    a = nodes["a"]
    b = nodes["b"]
    c = nodes["c"]
    d = nodes["d"]
    e = nodes["e"]
    f = nodes["f"]
    g = nodes["g"]
    h = nodes["h"]

    assert a.distance(a) == 0
    assert a.distance(b) == 0
    assert a.distance(c) == 3
    assert a.distance(d) == 4
    assert a.distance(e) == 7
    assert a.distance(f) == 5
    assert a.distance(g) == 8
    assert a.distance(h) == 2

    assert b.distance(a) == 0
    assert b.distance(b) == 0
    assert b.distance(c) == 3
    assert b.distance(d) == 4
    assert b.distance(e) == 7
    assert b.distance(f) == 5
    assert b.distance(g) == 8
    assert b.distance(h) == 2

    assert c.distance(a) == 3
    assert c.distance(b) == 3
    assert c.distance(c) == 0
    assert c.distance(d) == 1
    assert c.distance(e) == 4
    assert c.distance(f) == 2
    assert c.distance(g) == 5
    assert c.distance(h) == 5

    assert d.distance(a) == 4
    assert d.distance(b) == 4
    assert d.distance(c) == 1
    assert d.distance(d) == 0
    assert d.distance(e) == 5
    assert d.distance(f) == 3
    assert d.distance(g) == 6
    assert d.distance(h) == 6

    assert e.distance(a) == 7
    assert e.distance(b) == 7
    assert e.distance(c) == 4
    assert e.distance(d) == 5
    assert e.distance(e) == 0
    assert e.distance(f) == 6
    assert e.distance(g) == 9
    assert e.distance(h) == 9

    assert f.distance(a) == 5
    assert f.distance(b) == 5
    assert f.distance(c) == 2
    assert f.distance(d) == 3
    assert f.distance(e) == 6
    assert f.distance(f) == 0
    assert f.distance(g) == 3
    assert f.distance(h) == 7

    assert g.distance(a) == 8
    assert g.distance(b) == 8
    assert g.distance(c) == 5
    assert g.distance(d) == 6
    assert g.distance(e) == 9
    assert g.distance(f) == 3
    assert g.distance(g) == 0
    assert g.distance(h) == 10

    assert h.distance(a) == 2
    assert h.distance(b) == 2
    assert h.distance(c) == 5
    assert h.distance(d) == 6
    assert h.distance(e) == 9
    assert h.distance(f) == 7
    assert h.distance(g) == 10
    assert h.distance(h) == 0


def test_tip_to_tip_distances_endpoints(phylo_1: PhyloNode):
    """Test getting specifc tip distances  with tip_to_tip_distances"""
    names = ["H", "G", "M"]
    exp_dists = array([[0.0, 2.0, 6.7], [2.0, 0.0, 6.7], [6.7, 6.7, 0.0]])
    got_dists = phylo_1.tip_to_tip_distances()
    assert_equal(got_dists.take_dists(names).array, exp_dists)


def test_phylo_prune():
    """prune should reconstruct correct topology and Lengths of tree."""
    tree = DndParser("((a:3,((c:1):1):1):2);", constructor=PhyloNode)
    tree.prune()
    result_tree = DndParser("(a:3.0,c:3.0);", constructor=PhyloNode)
    assert str(tree) == str(result_tree)


def test_phylo_str(phylo_nodes: dict[str, PhyloNode]):
    """PhyloNode str should give expected results"""
    nodes = phylo_nodes
    a = nodes["a"]
    c = nodes["c"]
    f = nodes["f"]
    h = nodes["h"]

    assert str(f) == "(g:3):2"
    assert str(a) == "(((d:1,e:4,(g:3):2):3):0,h:2);"
    # check that None isn't converted any more
    h.length = None
    c.length = None  # need to test both leaf and internal node
    assert str(a) == "(((d:1,e:4,(g:3):2)):0,h);"


def test_get_max_tip_tip_distance(phylo_root: PhyloNode):
    """get_max_tip_tip_distance should get max tip distance across tree"""
    dist, names, node = phylo_root.get_max_tip_tip_distance()
    assert dist == 10.0
    assert set(names) == {"g", "h"}
    assert node.is_root()


def test_max_tip_tip_distance(phylo_root: PhyloNode):
    """max_tip_tip_distance returns the max dist between any pair of tips"""
    tree = phylo_root
    max_dist, tip_pair = tree.max_tip_tip_distance()
    assert max_dist == 10
    try:
        assert tip_pair == ("h", "g")
    except AssertionError:
        assert tip_pair == ("g", "h")


def test_tip_to_root_distances_phlylo():
    treestring = "(a:1,b:4,((c:3,d:2,e:1):5):2);"
    tree = make_tree(treestring)
    expect = {"a": 1, "b": 4, "c": 10, "d": 9, "e": 8}
    got = tree.tip_to_root_distances()
    assert got == expect
    got = tree.tip_to_root_distances(
        names=[
            "c",
            "a",
        ]
    )
    assert got == {"c": 10, "a": 1}


def test_tip_to_root_distances_tree():
    treestring = "(a:1,b:4,((c:3,d:2,e:2):5):2);"
    tree = DndParser(treestring, constructor=PhyloNode)
    expect = {"a": 1, "b": 1, "c": 3, "d": 3, "e": 3}
    got = tree.tip_to_root_distances(node_length=True)
    assert got == expect
    got = tree.tip_to_root_distances(
        names=[
            "c",
            "a",
        ],
        node_length=True,
    )
    assert got == {"c": 3, "a": 1}


def test_tip_to_root_distances_missing_names():
    treestring = "(a:1,b:4,((c:3,d:2,e:2):5):2);"
    tree = DndParser(treestring, constructor=PhyloNode)
    with pytest.raises(TreeError):
        tree.tip_to_root_distances(names=["x"])


@pytest.mark.parametrize(
    ("treestring", "check"),
    [
        ("(h:20,(d:1,(e:4,g:5):1):3.0);", any),
        ("(h:23,d:1,(e:4,g:5):1)a;", all),
        ("(e:4,g:5,(d:1,h:23):1)a;", all),
    ],
)
def test_root_at_midpoint2(treestring: str, check: Callable[[Iterable[object]], bool]):
    tree = make_tree(treestring)
    dists = array(
        [
            [0.0, 24.0, 28.0, 29.0],
            [24.0, 0.0, 6.0, 7.0],
            [28.0, 6.0, 0.0, 9.0],
            [29.0, 7.0, 9.0, 0.0],
        ]
    )
    names = ["h", "d", "e", "g"]

    result = tree.root_at_midpoint()
    assert check("h" in n.name for n in result.children)
    assert len(result.children) == 2
    got_dists = result.tip_to_tip_distances().take_dists(names)
    assert_equal(got_dists.array, dists)


def test_root_at_midpoint_already_rooted():
    tree = make_tree("(h:19,(d:1,(e:4,g:5):1):3.0);")
    result = tree.root_at_midpoint()
    assert any("h" in n.name for n in result.children)
    # already midpoint rooted
    assert result is tree
    node = result.get_node_matching_name("h")
    assert node.length == 14


def test_root_at_midpoint3():
    """midpoint between nodes should behave correctly"""
    tree = DndParser("(a:1,((c:1,d:2.5)mp:1,b:1)n2:1)root;")
    names = tree.get_tip_names()
    orig_dists = tree.tip_to_tip_distances().take_dists(names)
    tmid = tree.root_at_midpoint()
    assert tmid.is_root()
    # check the root children names are correct
    assert all("mp" in n.name for n in tmid.children)
    mp_dists = tmid.tip_to_tip_distances().take_dists(names)
    assert_equal(mp_dists.array, orig_dists.array)


def test_root_at_midpoint4():
    """midpoint should be selected correctly when it is an internal node"""
    # midpoint node is mp
    tree = DndParser("(a:1,((c:0.8,d:2.8)mp:1.2,b:1)n2:1);")
    tmid = tree.root_at_midpoint()
    assert tmid.is_root()
    # check the root children names are correct
    assert all("mp" in n.name for n in tmid.children)
    # check the dists are unchanged
    names = tree.get_tip_names()
    orig_dists = tree.tip_to_tip_distances().take_dists(names)
    mp_dists = tmid.tip_to_tip_distances().take_dists(names)
    assert_equal(
        mp_dists.array,
        orig_dists.array,
    )


def test_root_at_midpoint5():
    """midpoint should be selected correctly when on an even 2tip tree"""
    tree = DndParser("""(A:0.649351,B:0.649351):0.0;""")
    tmid = tree.root_at_midpoint()
    assert tmid is tree


def test_root_at_midpoint6():
    """already at midpoint"""
    tree = DndParser("(a:10,((c:1,d:3)n3:1,b:1)n2:1)root;")
    tmid = tree.root_at_midpoint()
    assert tmid is tree


def test_root_at_midpoint7():
    """midpoint is a tip"""
    tree = make_tree("(a:1,((c:1,d:7):1,b:1):1)root;")
    tmid = tree.root_at_midpoint()
    assert all("d" in n.name for n in tmid.children)
    # now check the root children lengths are correctly split
    mp = 5.0
    t2r = tmid.tip_to_root_distances()
    assert t2r["d"] == mp
    assert t2r["a"] == mp


def test_unrooted():
    """unrooted should preserve tips, drop a node"""
    rooted = make_tree(treestring="(B:0.2,(C:0.2,D:0.2)F:0.2)G;")
    unrooted = rooted.unrooted()
    assert sorted(rooted.get_tip_names()) == sorted(unrooted.get_tip_names())
    assert len(unrooted.get_node_names()) < len(rooted.get_node_names())


def test_get_figure():
    """exercising get_figure"""
    t_str = "(A:1,B:2,(C:3,D:3)E:2,(F,((G:1,H:2)I:2)J:3)K:2)L;"
    t = DndParser(t_str, constructor=PhyloNode)
    _ = t.get_figure(style="square")


@pytest.fixture
def root_std() -> PhyloNode:
    """Fixture to create a standard tree for testing."""
    return DndParser(tree_std, PhyloNode)


@pytest.fixture
def root_one_level() -> PhyloNode:
    """Fixture to create a one-level tree for testing."""
    return DndParser(tree_one_level, PhyloNode)


@pytest.fixture
def root_two_level() -> PhyloNode:
    """Fixture to create a two-level tree for testing."""
    return DndParser(tree_two_level, PhyloNode)


@pytest.fixture
def root_one_child() -> PhyloNode:
    """Fixture to create a tree with a single child for testing."""
    return DndParser(tree_one_child, PhyloNode)


def test_one_level(root_one_level: PhyloNode):
    """tip_to_tip should work for one-level multifurcating tree"""
    dmat = root_one_level.tip_to_tip_distances()
    assert_equal(dmat.array, array([[0, 3, 4], [3, 0, 5], [4, 5, 0]]))


def test_two_level(root_two_level: PhyloNode):
    """tip_to_tip should work for two-level tree"""
    dmat = root_two_level.tip_to_tip_distances()
    assert_allclose(
        dmat.array,
        array(
            [
                [0.0, 3.0, 4.0, 1.4],
                [3.0, 0.0, 5.0, 2.4],
                [4.0, 5.0, 0.0, 3.4],
                [1.4, 2.4, 3.4, 0.0],
            ]
        ),
    )


def test_std(root_std: PhyloNode):
    """tip_to_tip should work for small but complex tree"""
    dmat = root_std.tip_to_tip_distances()
    assert_equal(dmat.array, tree_std_dist)


def test_one_child(root_one_child: PhyloNode):
    """tip_to_tip should work for tree with a single child"""
    dmat = root_one_child.tip_to_tip_distances()
    assert_equal(dmat.array, tree_one_child_dist)


# for use with testing iterative copy method


def comb_tree(num_leaves: int) -> PhyloNode:
    """Returns a comb node_class tree."""
    branch_child = 1

    root = PhyloNode("")
    curr = root

    for _i in range(num_leaves - 1):
        curr.children[:] = [PhyloNode("", parent=curr), PhyloNode("", parent=curr)]
        curr = curr.children[branch_child]
    return root


# Moved  from test_tree2.py during code sprint on 04/14/10
# Missing tests: edge attributes (name, length, children) only get tested
# in passing by some of these tests.  See also xxx's

_default_newick = "((A:1,B:2)ab:3,((C:4,D:5)cd,E:6)cde:7)"


def _maketree(treestring: str = _default_newick) -> PhyloNode:
    return make_tree(treestring=treestring, underscore_unmunge=True)


@pytest.mark.parametrize(
    ("a", "b", "outgroup", "expect"),
    [
        ("A", "B", None, {"A", "B"}),
        ("E", "C", None, {"C", "D", "cd", "E"}),
        ("C", "D", "E", {"C", "D"}),
    ],
)
def test_get_edge_names(a: str, b: str, outgroup: str | None, expect: set[str]):
    tree = _maketree()
    assert (
        set(tree.get_edge_names(a, b, clade=True, stem=False, outgroup_name=outgroup))
        == expect
    )


def test_get_edge_names_2():
    tree = _maketree()
    clade = tree.get_edge_names("C", "E", stem=False, clade=True)
    clade.sort()
    assert clade == ["C", "D", "E", "cd"]

    names = tree.get_edge_names("C", "E", stem=True, clade=True)
    names.sort()
    assert names == ["C", "D", "E", "cd", "cde"]

    stem = tree.get_edge_names("C", "E", stem=True, clade=False)
    assert stem == ["cde"]


def test_get_edge_names_3():
    tree = make_tree(treestring="(a,b,(c,(d,e)e1)e2)")
    tip1name = "d"
    tip2name = "c"
    clade = False
    stem = True
    outgroup_name = "e"
    names = tree.get_edge_names(
        tip1name,
        tip2name,
        clade=clade,
        stem=stem,
        outgroup_name=outgroup_name,
    )
    assert names == ["e"]


@pytest.mark.parametrize(
    "treestring", ["((A,B)ab,(F,(C,D)cd)cdf,E)root;", "((E,(A,B)ab)abe,F,(C,D)cd)root;"]
)
def test_get_edge_names_use_outgroup(treestring: str):
    t = make_tree(treestring=treestring)
    # a, e, ogroup f
    expected = {"A", "B", "E", "ab"}
    edges = t.get_edge_names(
        "A",
        "E",
        stem=False,
        clade=True,
        outgroup_name="F",
    )
    assert set(edges) == expected


def test_get_connecting_node():
    tree = _maketree()
    assert tree.get_connecting_node("A", "B").name == "ab"
    assert tree.get_connecting_node("A", "C").name == "root"


def test_get_connecting_edges():
    """correctly identify connecting edges"""
    tree = make_tree(treestring="(((Human,HowlerMon)a,Mouse)b,NineBande,DogFaced);")
    edges = [e.name for e in tree.get_connecting_edges("Human", "Mouse")]
    assert set(edges) == {"Human", "Mouse", "a"}

    edges = [e.name for e in tree.get_connecting_edges("b", "Human")]
    assert set(edges) == {"Human", "a", "b"}


@pytest.mark.parametrize(("name", "expect_tip"), [("A", True), ("ab", False)])
def test_get_node_matching_name_2(name: str, expect_tip: bool):
    tree = _maketree()
    edge = tree.get_node_matching_name(name)
    assert edge.name == name
    assert edge.istip() == expect_tip


def test_get_edge_vector():
    """correctly return vector of edges from a tree"""
    tree = _maketree()
    names = [e.name for e in tree.get_edge_vector()]
    assert names == ["A", "B", "ab", "C", "D", "cd", "E", "cde", "root"]

    names = [e.name for e in tree.get_edge_vector(include_root=False)]
    assert names == ["A", "B", "ab", "C", "D", "cd", "E", "cde"]


def test_get_newick_2():
    orig = "((A:1.0,B:2.0)ab:3.0,((C:4.0,D:5.0)cd:6.0,E:7.0)cde:8.0)all;"
    unlen = "((A,B)ab,((C,D)cd,E)cde)all;"
    tree = _maketree(orig)
    assert (
        tree.get_newick(with_distances=True, with_node_names=True, with_root_name=True)
        == orig
    )
    assert tree.get_newick(with_node_names=True, with_root_name=True) == unlen

    tree.name = "a'l"
    ugly_name = "((A,B)ab,((C,D)cd,E)cde)a'l;"
    ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a''l';"
    assert (
        tree.get_newick(escape_name=True, with_node_names=True, with_root_name=True)
        == ugly_name_esc
    )
    assert (
        tree.get_newick(escape_name=False, with_node_names=True, with_root_name=True)
        == ugly_name
    )

    tree.name = "'a l'"
    quoted_name = "((A,B)ab,((C,D)cd,E)cde)'a l';"
    quoted_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a l';"
    assert (
        tree.get_newick(escape_name=True, with_node_names=True, with_root_name=True)
        == quoted_name_esc
    )
    assert (
        tree.get_newick(escape_name=False, with_node_names=True, with_root_name=True)
        == quoted_name
    )


def test_str_1():
    """testing (well, exercising at least), __str__"""
    got = str(_maketree())
    assert got.count("(") == got.count(")")


def test_repr():
    """testing (well, exercising at least), __repr__"""
    got = repr(_maketree())
    assert isinstance(got, str)
    assert got.count("(") == got.count(")")


def test_eq():
    """testing (well, exercising at least), __eq__"""
    t1 = _maketree()
    t1_again = t1

    t2 = _maketree()
    assert t1 == t1_again
    assert t1 != t2


def test_balanced():
    """balancing an unrooted tree"""
    t = make_tree(treestring="((a,b),((c1,(c2,(c3,(c4,(c5,(c6,c7)))))),(d,e)),f)")
    b = make_tree(treestring="(c1,(c2,(c3,(c4,(c5,(c6,c7))))),((d,e),((a,b),f)))")
    assert str(t.balanced()) == str(b)


def test_params_merge():
    t = make_tree(treestring="((((a:1,b:3)ab:4,c)abc),d)")
    params: list[tuple[str, float]] = [("a", 20), ("b", 2.0), ("ab", 5.0)]
    for label, beta in params:
        t.get_node_matching_name(label).params |= {"beta": beta}
    t = t.get_sub_tree(["b", "c", "d"])
    # previous implementation on merge was
    # dividing the sum of parameters across nodes by the lengths
    # we no longer try and support this
    assert t.get_node_matching_name("b").params == {
        "beta": 2.0,
    }
    assert t.get_node_matching_name("b").length == 7

    ex = make_tree("(b:7,c);")
    st = t.get_sub_tree(["b", "c", "xxx"], ignore_missing=True)
    assert ex.same_topology(st)
    with pytest.raises(ValueError):
        # should raise ValueError if a tip is not found
        t.get_sub_tree(["b", "c", "xxx"])


def test_making_from_list():
    tipnames_with_spaces = {"a_b", "a b", "T'lk"}
    t = make_tree(tip_names=list(tipnames_with_spaces))
    result = t.get_tip_names()
    assert set(result) == tipnames_with_spaces


def test_getset_param_value():
    """test getting, setting of param values"""
    t = make_tree(treestring="((((a:.2,b:.3)ab:.1,c:.3)abc:.4),d:.6)")
    t.set_param_value("zz", "ab", 4.321)
    assert t.get_param_value("zz", "ab") == 4.321
    node = t.get_node_matching_name("ab")
    assert node.params["zz"] == 4.321


def test_rootswaps():
    """testing (well, exercising at least), unrooted"""
    new_tree = make_tree(treestring="((a,b),(c,d))")
    new_tree = new_tree.unrooted()
    assert len(new_tree.children) > 2, "not unrooted right"


def test_reroot():
    tree = make_tree(treestring="((a,b),(c,d),e)")
    tree2 = tree.rooted_with_tip("b")
    assert tree2.get_newick() == "(a,b,((c,d),e));"


def test_same_shape():
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


@pytest.fixture
def otu_names_small() -> list[str]:
    return sorted(["NineBande", "Mouse", "HowlerMon", "DogFaced"])


@pytest.fixture
def newick_small() -> str:
    return "(((Human,HowlerMon),Mouse),NineBande,DogFaced);"


@pytest.fixture
def newick_small_reduced() -> str:
    return "((HowlerMon,Mouse),NineBande,DogFaced);"


@pytest.fixture
def newick_small_sorted() -> str:
    return "(DogFaced,((HowlerMon,Human),Mouse),NineBande);"


@pytest.fixture
def tree_small(newick_small: str):
    return make_tree(treestring=newick_small)


@pytest.fixture
def otu_names_big() -> list[str]:
    return sorted(
        [
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
    )


@pytest.fixture
def newick_big_reduced() -> str:
    return "(((((TombBat,(((Horse,Rhino),(Pangolin,Cat)),(Pig,SpermWhal))),Hedgehog),(Orangutan,Gorilla)),((HairyArma,Sloth),((Manatee,AsianElep),GoldenMol))),bandicoot);"


@pytest.fixture
def tree_big() -> PhyloNode:
    return make_tree(
        "((((((((FlyingFox,DogFaced),((FreeTaile,LittleBro),(TombBat,RoundEare))),(FalseVamp,LeafNose)),(((Horse,Rhino),(Pangolin,(Cat,Dog))),(Llama,(Pig,(Cow,(Hippo,(SpermWhal,HumpbackW))))))),(Mole,Hedgehog)),(TreeShrew,(FlyingLem,((Jackrabbit,(FlyingSqu,(OldWorld,(Mouse,Rat)))),(Galago,(HowlerMon,(Rhesus,(Orangutan,(Gorilla,(Human,Chimpanzee)))))))))),(((NineBande,HairyArma),(Anteater,Sloth)),(((Dugong,Manatee),((AfricanEl,AsianElep),(RockHyrax,TreeHyrax))),(Aardvark,((GoldenMol,(Madagascar,Tenrec)),(LesserEle,GiantElep)))))),(caenolest,(phascogale,(wombat,bandicoot))));"
    )


def test_sorttree(tree_small: PhyloNode, newick_small_sorted: str):
    """testing (well, exercising at least) treesort"""
    new_tree = tree_small.sorted()
    assert newick_small_sorted == new_tree.get_newick(with_distances=False)


@pytest.mark.parametrize(
    ("tree_fxt", "otu_fxt", "nwk_fxt"),
    [
        ("tree_small", "otu_names_small", "newick_small_reduced"),
        ("tree_big", "otu_names_big", "newick_big_reduced"),
    ],
)
def test_getsubtree(
    request: pytest.FixtureRequest,
    tree_fxt: str,
    otu_fxt: str,
    nwk_fxt: str,
):
    """testing getting a subtree"""
    tree = request.getfixturevalue(tree_fxt)
    otu_names = request.getfixturevalue(otu_fxt)
    newick_reduced = request.getfixturevalue(nwk_fxt)
    tree = tree.unrooted()
    subtree = tree.get_sub_tree(otu_names)

    exp_tree = make_tree(treestring=newick_reduced).unrooted()

    # names match
    assert set(subtree.get_tip_names()) == set(exp_tree.get_tip_names())
    # topology matches
    assert exp_tree.same_topology(subtree)
    # newick correctly formed means we reconstruct the same topology
    from_nwk = make_tree(subtree.get_newick())
    assert exp_tree.same_topology(from_nwk)


def test_getsubtree_2():
    """tree.get_sub_tree() has same pairwise tip dists as tree (len0 node)"""
    t1 = DndParser(
        "((a:1,b:2):4,((c:3, j:17.2):0,(d:1,e:1):2):3)",
        PhyloNode,
    )  # note c,j is len 0 node
    names = "a", "b", "d", "e", "c"
    orig = t1.tip_to_tip_distances().take_dists(names)
    subtree = t1.get_sub_tree(names)
    subset = subtree.tip_to_tip_distances().take_dists(names)
    assert_equal(orig.array, subset.array)


def test_ascii():
    # unlabeled internal node
    tr = DndParser("(B:0.2,(C:0.3,D:0.4):0.6)F;")
    obs = tr.ascii_art(show_internal=True, compact=False)
    exp = "          /-B\n-F-------|\n         |          /-C\n          \\--------|\n                    \\-D"
    assert obs == exp
    obs = tr.ascii_art(show_internal=True, compact=True)
    exp = "-F------- /-B\n          \\-------- /-C\n                    \\-D"
    assert obs == exp
    obs = tr.ascii_art(show_internal=False, compact=False)
    exp = "          /-B\n---------|\n         |          /-C\n          \\--------|\n                    \\-D"
    assert obs == exp


@pytest.mark.parametrize("tree_fxt", ["tree_small", "tree_big"])
def test_load_tree(
    tmp_path: pathlib.Path, tree_fxt: str, request: pytest.FixtureRequest
):
    """tests loading a newick formatted Tree"""
    tree_path = tmp_path / "tree.tree"
    tree = request.getfixturevalue(tree_fxt)
    tree.write(tree_path)
    got = load_tree(tree_path)
    assert isinstance(got, PhyloNode)
    assert got.get_newick() == tree.get_newick()
    assert set(got.get_node_names()) == set(tree.get_node_names())


@pytest.mark.parametrize("tree_fxt", ["tree_small", "tree_big"])
def test_load_tree_from_json(
    tmp_path: pathlib.Path, tree_fxt: str, request: pytest.FixtureRequest
):
    """tests loading a Tree from json file"""
    json_path = tmp_path / "tree.json"
    tree = request.getfixturevalue(tree_fxt)
    tree.write(json_path)
    got = load_tree(json_path)
    assert isinstance(got, PhyloNode)
    assert got.get_newick() == tree.get_newick()
    assert set(got.get_node_names()) == set(tree.get_node_names())
    # now try using non json suffix
    json_path = tmp_path / "tree.txt"
    tree.write(json_path, format_name="json")
    got = load_tree(json_path, format_name="json")
    assert isinstance(got, PhyloNode)


def test_getsubtree_5():
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

    sub1 = tree.get_sub_tree(names, tips_only=False)
    assert tree.same_topology(sub1)
    expect = make_tree(
        treestring="(Homo_sapiens,Mus_musculus,"
        "(Canis_familiaris,Ornithorhynchus_anatinus))",
        underscore_unmunge=True,
    )
    sub2 = tree.get_sub_tree(names, tips_only=True)
    assert expect.same_topology(sub2)


def test_getsubtree_6():
    """get sub tree handles non-scalar params"""
    names = ["A", "B", "C"]
    treestring = "((E:2,(F:1,(C:0.1,D:0.2):0.3):0.3):0.3,A:0.2,B:0.2)"
    tree = make_tree(treestring=treestring)
    # change internal node names to eliminate ending digits
    vals = {"A": 42, "B": 43, "C": 44}
    for edge in tree.postorder(include_self=False):
        edge.params["non-scalar"] = {edge.name: vals.get(edge.name, 99)}

    sub1 = tree.get_sub_tree(names, tips_only=False)
    assert set(tree.get_tip_names()), set("ABC")
    # the edge value for "non-scalar" should be same as original tree
    assert all(sub1.get_node_matching_name(name).params["non-scalar"] for name in names)


def test_get_subtree_7():
    t = make_tree(treestring="((((a:1,b:3)ab:4,c)abc),d)")
    st = t.get_sub_tree(["b", "c"])
    assert st.is_root()


@pytest.mark.parametrize("treestring", ["(a,b,(c,(d,e)))", "((a,b),(c,(d,e)))"])
def test_get_subtree_8(treestring: str):
    # an unrooted tree returns an unrooted subtree
    tree = make_tree(treestring=treestring)
    st = tree.get_sub_tree(["c", "d", "e"], as_rooted=False)
    assert len(st.children) == len(tree.children)


def test_get_subtree_as_rooted():
    # a rooted tree with subtree being naturally multifurcating
    tree = make_tree("((a,b),((c,d),(e,f),(g,h)));")

    got = tree.get_sub_tree(["d", "e", "h"], as_rooted=True)
    assert len(got.children) == 3
    got = tree.get_sub_tree(["d", "e", "f"], as_rooted=False)
    assert len(got.children) == len(tree.children)


def test_get_edge_names_big(tree_big: PhyloNode, otu_names_big: list[str]):
    """testing (well, exercising at least), getedgenames"""
    # Fell over on small tree because "stem descended from root
    # joiner was a tip"
    a, b = otu_names_big[:2]
    got = tree_big.get_edge_names(a, b, clade=True, stem=False)
    assert {type(n) for n in got} == {str}


def test_get_tip_names(tree_big: PhyloNode):
    """testing (well, exercising at least), get_tip_names"""
    tips = tree_big.get_tip_names()
    assert len(tips) == 55


@pytest.mark.parametrize("num_tips", [0, 3])
def test_get_distances_endpoints(num_tips: int):
    names = "G001280565", "G000009905", "G000183545", "G000184705"
    nwk = "({}:0.4,({}:0.25,({}:0.03,{}:0.03):0.23):0.03)".format(*names)
    tree = make_tree(nwk)
    dists = tree.tip_to_tip_distances(names=names[:num_tips] or None)
    num = 4 if num_tips == 0 else num_tips
    assert dists.shape == (num, num)


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("lin_rajan_moret", 3),
        ("lrm", 3),
        ("matching", 3),
        (None, 3),
        ("unrooted_robinson_foulds", 4),
        ("urf", 4),
        ("rf", 4),
    ],
)
def test_tree_distance_unrooted(method: str | None, expected: int):
    tree1 = make_tree(treestring="(a,b,(c,(d,e)));")
    tree2 = make_tree(treestring="((a,c),(b,d),e);")

    assert tree1.tree_distance(tree2, method=method) == expected


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("matching_cluster", 10),
        ("mc", 10),
        ("matching", 10),
        (None, 10),
        ("rooted_robinson_foulds", 6),
        ("rrf", 6),
        ("rf", 6),
    ],
)
def test_tree_distance_rooted(method: str | None, expected: int):
    tree1 = make_tree(treestring="(a,(b,(c,(d,e))));")
    tree2 = make_tree(treestring="(e,(d,(c,(b,a))));")

    assert tree1.tree_distance(tree2, method=method) == expected


@pytest.mark.parametrize(
    "bad_method",
    [
        "matching_cluster ",
        "m",
        "mcc",
        "rr",
        "rff",
        "ur",
        "match",
        " matching",
        "robinson_foulds",
        "None",
    ],
)
def test_tree_distance_method_does_not_exist(bad_method: str):
    unrooted1 = make_tree(treestring="(a,b,(c,(d,e)));")
    unrooted2 = make_tree(treestring="((a,c),(b,d),e);")

    with pytest.raises(ValueError):
        unrooted1.tree_distance(unrooted2, method=bad_method)

    rooted1 = make_tree(treestring="(a,(b,(c,(d,e))));")
    rooted2 = make_tree(treestring="(e,(d,(c,(b,a))));")

    with pytest.raises(ValueError):
        rooted1.tree_distance(rooted2, method=bad_method)


@pytest.mark.parametrize(
    "method",
    [
        "matching_cluster",
        "lin_rajan_moret",
        "rooted_robinson_foulds",
        "unrooted_robinson_foulds",
        "mc",
        "lrm",
        "rrf",
        "urf",
        "matching",
        "rf",
        None,
    ],
)
def test_tree_distance_incompatible_trees(method: str):
    unrooted = make_tree(treestring="(a,b,(c,(d,e)));")
    rooted = make_tree(treestring="(a,(b,(c,(d,e))));")

    with pytest.raises(ValueError):
        rooted.tree_distance(unrooted, method=method)

    with pytest.raises(ValueError):
        unrooted.tree_distance(rooted, method=method)


def test_lrm_method():
    # this test just exercises the method, the tests on the underlying
    # function are in test_tree_distance.py
    a = make_tree(treestring="(1,(((2,3),4),(5,((6,(7,(8,9))),(10,11)))),12);")
    b = make_tree(treestring="(1,((((2,3),4),5),((6,7),((8,9),(10,11)))),12);")
    distance = a.lin_rajan_moret(b)
    assert distance == 8


@pytest.mark.parametrize(
    "nwk",
    ["(((g,b)gb,(c,d)cd),(e,f),a)", "(a,b,c)", ("(a,(b,(c,d)cd))")],
)
def test_child_parent_map(nwk: str):
    tree = make_tree(nwk)
    child_2_parent = tree.child_parent_map()
    all_edges = tree.get_edge_vector(include_root=False)
    assert len(child_2_parent) == len(all_edges)
    assert child_2_parent["a"] == "root"
    edge = random.choice(all_edges)
    assert child_2_parent[edge.name] == edge.parent.name


def test_parser():
    """nasty newick"""
    nasty = "( (A :1.0,'B (b)': 2) [com\nment]pair:3,'longer name''s':4)dash_ed;"
    nice = "((A:1.0,'B (b)':2.0)pair:3.0,'longer name''s':4.0)dash_ed;"
    tree = make_tree(treestring=nasty, underscore_unmunge=True)
    tidied = tree.get_newick(
        with_distances=True, with_node_names=True, with_root_name=True
    )
    assert tidied == nice
    assert tree.get_node_matching_name("pair").params["other"] == ["com\nment"]


def test_load_tree_bad_encoding():
    """
    charset_normalizer rather than chardet as a dependency
    incorrectly detected the encoding of a file as UTF-16LE.
    """
    newick = "(a,b);"

    with TemporaryDirectory(dir=".") as dirname:
        tree_path = pathlib.Path(dirname) / "tree.tree"
        with tree_path.open("wb") as f:
            f.write(newick.encode("ascii"))

        assert load_tree(tree_path).get_newick() == newick


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        (None, (None, None)),
        ("", (None, None)),
        ("edge.98/24", ("edge.98", 24.0)),
        ("edge.98", ("edge.98", None)),
        ("24", (None, 24.0)),
    ],
)
def test_split_name_and_support(
    name: str | None, expected: tuple[str | None, float | None]
):
    assert split_name_and_support(name) == expected


@pytest.mark.parametrize("invalid", ["edge.98/invalid", "edge.98/23/invalid"])
def test_split_name_and_support_invalid_support(invalid: str):
    with pytest.raises(ValueError):
        split_name_and_support(invalid)


def test_phylonode_support():
    tip_names = [str(i) for i in range(1, 13)]
    tree = make_tree(
        treestring="(1,(((2,3)53,4)/53,(5,((6,(7,(8,9))def/25),(10,11)abc))),12);",
    )
    assert tree.get_tip_names() == tip_names

    node_names = set(tree.get_node_names())
    # check that no node name is an empty string or None
    assert all(node_names)

    # parent of 4 is node with only support value
    just_support = tree.get_node_matching_name("4").parent
    assert just_support.support == 53.0

    # parent of 2 has the same support as parent of 4
    same_support = tree.get_node_matching_name("2").parent
    assert same_support.support == just_support.support

    # parent of 10 has a node name only
    just_name = tree.get_node_matching_name("10").parent
    assert just_name.name == "abc"
    assert just_name.support is None

    # the node with name "def/25" correctly resoloved into node
    # name "def" and support 25.0
    name_and_support = tree.get_node_matching_name("def")
    assert name_and_support.name == "def"  # bit redundant given selection process
    assert name_and_support.support == 25.0


def test_source_attr():
    t_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
    t = DndParser(t_str, constructor=PhyloNode)
    t.source = "test_source"
    assert t.source == "test_source"
    assert "source" in t.params
    assert "source" not in t.get_node_matching_name("a").params
    t.source = None
    assert t.source is None
    assert "source" not in t.params
    t.source = ""
    assert t.source is None


def test_dataset_tree_source():
    tr = get_dataset("mammal-tree")
    assert tr.source.endswith("murphy.tree")


def test_rooted_is_rooted_diff_root_name():
    treestring = "(c,b)abc"
    tree = make_tree(treestring=treestring)
    rooted = tree.rooted("b")
    assert rooted.same_topology(tree)


def test_rooted_is_rooted():
    treestring = "((Human,Chimpanzee),Gorilla)"
    tree = make_tree(treestring=treestring)
    rooted = tree.rooted("root")
    assert rooted.same_topology(tree)


def test_rooted_with_root():
    treestring = "(Human,Chimpanzee,Gorilla)"
    tree = make_tree(treestring=treestring)

    assert tree.same_shape(tree.rooted("root"))


def test_rooted_from_nonroot():
    treestring = "(Human,Chimpanzee,Gorilla)"
    tree = make_tree(treestring=treestring)
    child = tree.children[0]
    with pytest.raises(TreeError):
        # cannot root if node is not the root
        child.rooted("Gorilla")


def test_rooted_with_tip():
    treestring = "((((Rhesus:0.02,HowlerMon:0.0487):0.0115,Orangutan:0.008)abcd:0.003,Human:0.006):0.00013,Gorilla:0.0025,Chimpanzee:0.003)"
    tree = make_tree(treestring=treestring)
    tip_name = "HowlerMon"
    rooted = tree.rooted(tip_name)
    assert len(rooted.children) == 2
    assert all(child.length == 0.0487 / 2 for child in rooted.children)
    assert {c.name for c in rooted.children} == {tip_name, f"{tip_name}-root"}


def tested_rooted_newick():
    tree = make_tree(treestring="(A:0.006,B:0.0025,C:0.003)")
    nt = tree.rooted("A")
    got = str(nt.get_newick(with_distances=False))
    assert got in {"(A,(B,C));", "((B,C),A);", "((C,B),A);"}


def test_rooted_with_internal():
    treestring = "((((Rhesus:0.02,HowlerMon:0.0487):0.0115,Orangutan:0.008)abcd:0.003,Human:0.006):0.00013,Gorilla:0.0025,Chimpanzee:0.003)"
    tree = DndParser(treestring, constructor=PhyloNode)
    tree.name = "root"  # have to force this given construction via DndParser
    node_name = "abcd"
    got = tree.rooted(node_name)
    ts1 = "(Human,(Gorilla,Chimpanzee))"
    ts2 = "((Rhesus,HowlerMon),Orangutan)"
    child1, child2 = got.children
    child1.parent = None
    child1.name = "root"
    child2.parent = None
    child2.name = "root"
    ts1, ts2 = (ts1, ts2) if "Human" in child1.get_tip_names() else (ts2, ts1)
    expect1 = make_tree(treestring=ts1)
    expect2 = make_tree(treestring=ts2)
    assert child1.same_topology(expect1)
    assert child2.same_topology(expect2)


def test_unrooted_deepcopy():
    treestring = "(a:1,((c:1,d:2.5)n3:1,b:1)n2:1)rt;"
    tree = make_tree(treestring=treestring)
    node = tree.get_node_matching_name("n3")
    dc = node.unrooted_deepcopy()
    assert all(len(n.children) > 1 for n in dc.iter_nontips())


def test_unrooted_deepcopy_2():
    tree = make_tree(treestring="(a,b,(d,e)de);")
    tip = tree.get_node_matching_name("e")
    unco = tip.unrooted_deepcopy()
    # the root should have one child whose name is 'e'
    assert len(unco.children) == 1
    assert unco.children[0].name == "e"


def test_load_old_tree_json(DATA_DIR: pathlib.Path):
    with pytest.warns(UserWarning, match="Outdated tree json"):
        tree_old = load_tree(DATA_DIR / "tree_json_old.json")
    tree_new = load_tree(DATA_DIR / "tree_json_new.json")

    for old, new in zip(tree_old.preorder(), tree_new.preorder(), strict=True):
        assert old.name == new.name
        if old.length is None:
            assert new.length is None
        else:
            assert pytest.approx(old.length) == new.length

        if old.support is None:
            assert new.length is None
        else:
            assert pytest.approx(old.support) == new.support


def test_support_carries_through_sorted():
    tree = make_tree("(a,((c,d)70,b)30)")
    tree = tree.sorted()
    assert tree[1].support == 30
    assert tree[1][1].support == 70


def test_conflicting_support():
    name_with_support = "innerName/50"

    with pytest.raises(ValueError):
        TreeBuilder(PhyloNode).create_edge(
            [PhyloNode("a"), PhyloNode("b")], name_with_support, {}, None, 30
        )

import pathlib

import pytest

from cogent3 import make_tree
from cogent3.phylo.tree_distance import lin_rajan_moret


def test_different_trees():
    a = make_tree(treestring="(1,(((2,3),4),(5,((6,(7,(8,9))),(10,11)))),12);")
    b = make_tree(treestring="(1,((((2,3),4),5),((6,7),((8,9),(10,11)))),12);")
    distance = lin_rajan_moret(a, b)
    assert distance == 8


def test_dist_matches_original_implementation_small_tree(DATA_DIR):
    path, expect = DATA_DIR / "match_dist_100.tree", 1402
    data = path.read_text().splitlines()
    t1 = make_tree(data[0])
    t2 = make_tree(data[1])
    distance = lin_rajan_moret(t1, t2)
    # Assert that this is same value as original C implementation
    assert distance == expect


@pytest.mark.slow
def test_dist_matches_original_implementation_bigtree(DATA_DIR):
    path, expect = DATA_DIR / "match_dist_2000.tree", 17037
    data = pathlib.Path(path).read_text().splitlines()
    t1 = make_tree(data[0])
    t2 = make_tree(data[1])
    distance = lin_rajan_moret(t1, t2)
    # Assert that this is same value as original C implementation
    assert distance == expect


def test_same_tree():
    a = make_tree(treestring="(a, b, (c, d));")
    b = make_tree(treestring="(c, d, (a, b));")
    distance = lin_rajan_moret(a, b)
    assert distance == 0


@pytest.mark.parametrize(
    "b",
    (
        "(e, d, (a, b));",
        "(a, d, (a, b))",
        "(e, d, (a, b));",
        "((e, d), (a, b));",
        "((e, d), (a, b));",
    ),
)
def test_fails_different_tips(b):
    t1 = make_tree(treestring="(a, b, (c, d));")
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        lin_rajan_moret(t1, t2)


@pytest.mark.parametrize(
    "a,b",
    (
        ("((a, b), (c, d));", "(e, d, (a, b));"),
        ("(a, b, (c, d));", "((e, d), (a, b));"),
        ("((a, b), (c, d));", "((e, d), (a, b));"),
    ),
)
def test_fails_unrooted(a, b):
    t1 = make_tree(treestring=a)
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        lin_rajan_moret(t1, t2)

import pathlib

import pytest

from cogent3 import make_tree
from cogent3.phylo.tree_distance import (
    lin_rajan_moret,
    matching_cluster_distance,
    rooted_robinson_foulds,
    unrooted_robinson_foulds,
)


ROOTED_DISTANCE_MEASURES = matching_cluster_distance, rooted_robinson_foulds
UNROOTED_DISTANCE_MEASURES = lin_rajan_moret, unrooted_robinson_foulds


# unrooted rf distance
def test_unrooted_rf_different_trees():
    a = make_tree(treestring="(a, b, (c, (d, e)));")

    # There are two splits not in common for the two trees
    # abc-de in tree a, and ad-dce in tree b.
    b = make_tree(treestring="(a, b, (d, (c, e)));")
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 2
    distance = unrooted_robinson_foulds(b, a)
    assert distance == 2

    # Identical unrooted tree as before, but a different orientation.
    b = make_tree(treestring="((a, b), d, (c, e));")
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 2
    distance = unrooted_robinson_foulds(b, a)
    assert distance == 2

    # Identical unrooted tree as before, but a different orientation.
    b = make_tree(treestring="(((a, b), d), c, e);")
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 2
    distance = unrooted_robinson_foulds(b, a)
    assert distance == 2

    a = make_tree(treestring="(a, (b, (c, d)), (e, (f, (g, h))));")
    b = make_tree(treestring="(a, (c, (b, d)), (h, (f, (g, e))));")

    # There are six splits not in common.
    # For tree a they are cd-abefgh, abcdef-gh, and abcde-fgh
    # For tree b they are bd-acefgh, abcdhf-ge, and abcdh-fge
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 6
    distance = unrooted_robinson_foulds(b, a)
    assert distance == 6


def test_unrooted_rf_same_tree():
    # Topology of the unrooted trees are identical
    # though they are stored in different orientations.
    a = make_tree(treestring="(a, b, (c, d));")
    b = make_tree(treestring="(c, d, (a, b));")

    distance = unrooted_robinson_foulds(a, b)
    assert distance == 0

    distance = unrooted_robinson_foulds(b, a)
    assert distance == 0


def test_unrooted_rf_distance_multifurcation():
    a = make_tree(treestring="(a, b, (c, (d, e)));")
    b = make_tree(treestring="(a, b, (c, d, e));")

    # There is one split not in common.
    # Tree a: abc-de
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 1

    distance = unrooted_robinson_foulds(b, a)
    assert distance == 1

    a = make_tree(treestring="(a, b, (((c1, c2), (d1, (d2, d3))), (e1, (e2, e3))));")
    b = make_tree(treestring="(a, b, ((c1, c2), (d2, (d1, d3)), (e1, e2, e3)));")

    # There are four splits not in common
    # Tree a: e1-others, d2,d3-others, c's and d's - others
    # Tree b: d1,d3-others
    distance = unrooted_robinson_foulds(a, b)
    assert distance == 4
    distance = unrooted_robinson_foulds(b, a)
    assert distance == 4


# lin_rajan_moret distance
def test_lrm_different_trees():
    a = make_tree(treestring="(1,(((2,3),4),(5,((6,(7,(8,9))),(10,11)))),12);")
    b = make_tree(treestring="(1,((((2,3),4),5),((6,7),((8,9),(10,11)))),12);")
    distance = lin_rajan_moret(a, b)
    assert distance == 8


def test_lrm_matches_original_implementation_small_tree(DATA_DIR):
    path, expect = DATA_DIR / "match_dist_100.tree", 1402
    data = path.read_text().splitlines()
    t1 = make_tree(data[0])
    t2 = make_tree(data[1])
    distance = lin_rajan_moret(t1, t2)
    # Assert that this is same value as original C implementation
    assert distance == expect


@pytest.mark.slow
def test_lrm_matches_original_implementation_big_tree(DATA_DIR):
    path, expect = DATA_DIR / "match_dist_2000.tree", 17037
    data = pathlib.Path(path).read_text().splitlines()
    t1 = make_tree(data[0])
    t2 = make_tree(data[1])
    distance = lin_rajan_moret(t1, t2)
    # Assert that this is same value as original C implementation
    assert distance == expect


def test_lrm_same_tree():
    # Topology of the unrooted trees are identical
    # though they are stored in different orientations.
    a = make_tree(treestring="(a, b, (c, d));")
    b = make_tree(treestring="(c, d, (a, b));")

    distance = lin_rajan_moret(a, b)
    assert distance == 0

    distance = lin_rajan_moret(b, a)
    assert distance == 0


# rooted rf distance
def test_rooted_rf_multifurcation():
    a = make_tree(treestring="((a, b), (c, d));")
    b = make_tree(treestring="((a, c, b), d);")

    # There are three clades not in common
    # Tree a: ab, cd
    # Tree b: acb
    distance = rooted_robinson_foulds(a, b)
    assert distance == 3

    distance = rooted_robinson_foulds(b, a)
    assert distance == 3


def test_rooted_rf_different_trees():
    a = make_tree(
        treestring="((((a1,a2),(b1,b2,b3),x),(c1,c2,c3),(d1,d2)),(e1,(e2,e3)));"
    )
    b = make_tree(
        treestring="((((a1,a2),(b1,b2,b3)),x,(c1,c2,c3),(d1,d2)),(e1,(e2,e3)));"
    )

    # Moving the single x up one level creates two clades that are not in common.
    # Tree a: ABx
    # Tree b: AB
    distance = rooted_robinson_foulds(a, b)
    assert distance == 2

    distance = rooted_robinson_foulds(b, a)
    assert distance == 2

    a = make_tree(treestring="((a,b),c);")
    b = make_tree(treestring="((a,c),b);")

    # There are two clades not in common
    # Tree a: ab
    # Tree b: ac
    distance = rooted_robinson_foulds(a, b)
    assert distance == 2

    distance = rooted_robinson_foulds(b, a)
    assert distance == 2

    a = make_tree(treestring="((a,(b1,(b2,b3))),((c1,c2),(c3,c4)));")
    b = make_tree(treestring="((a,((c1,c2),(c3,c4))),(b1,(b2,b3)));")

    # There are two clades not in common
    # Tree a: aB
    # Tree b: aC
    distance = rooted_robinson_foulds(a, b)
    assert distance == 2

    distance = rooted_robinson_foulds(b, a)
    assert distance == 2

    a = make_tree(treestring="(a,(b,(c,(d,(e,f)))));")
    b = make_tree(treestring="(a,(f,(c,(d,(b,e)))));")

    # There are six clades not in common
    # Tree a: ef, def, cdef
    # Tree b: be, bde, bcde
    distance = rooted_robinson_foulds(a, b)
    assert distance == 6

    distance = rooted_robinson_foulds(b, a)
    assert distance == 6


def test_rooted_rf_same_tree():
    a = make_tree(treestring="((a, b), (c, d));")
    b = make_tree(treestring="((c, d), (a, b));")

    # The distance between trees of the same topology is 0
    distance = rooted_robinson_foulds(a, b)
    assert distance == 0


# matching cluster distance
def test_matching_cluster_distance_multifurcation():
    # The basic example from the matching cluster distance paper
    # see reference in matching_cluster_distance
    a = make_tree(treestring="((a, b), (c, d));")
    b = make_tree(treestring="((a, c, b), d);")

    distance = matching_cluster_distance(a, b)
    assert distance == 3

    distance = matching_cluster_distance(b, a)
    assert distance == 3


def test_matching_cluster_distance_properties():
    # Test properties of the matching cluster distance paper shown as exampels
    # see reference in matching_cluster_distance
    a = make_tree(
        treestring="((((a1,a2),(b1,b2,b3),x),(c1,c2,c3),(d1,d2)),(e1,(e2,e3)));"
    )
    b = make_tree(
        treestring="((((a1,a2),(b1,b2,b3)),x,(c1,c2,c3),(d1,d2)),(e1,(e2,e3)));"
    )

    distance = matching_cluster_distance(a, b)
    assert distance == 1

    distance = matching_cluster_distance(b, a)
    assert distance == 1

    a = make_tree(treestring="((a,b),c);")
    b = make_tree(treestring="((a,c),b);")

    distance = matching_cluster_distance(a, b)
    assert distance == len(a.tips()) - 1

    distance = matching_cluster_distance(b, a)
    assert distance == len(a.tips()) - 1

    a = make_tree(treestring="((a,(b1,(b2,b3))),((c1,c2),(c3,c4)));")
    b = make_tree(treestring="((a,((c1,c2),(c3,c4))),(b1,(b2,b3)));")

    distance = matching_cluster_distance(a, b)
    assert distance == len(a.tips()) - 1

    distance = matching_cluster_distance(b, a)
    assert distance == len(a.tips()) - 1


def test_matching_cluster_distance_same_tree():
    a = make_tree(treestring="((a, b), (c, d));")
    b = make_tree(treestring="((c, d), (a, b));")

    # The distance between trees of the same topology is 0
    distance = matching_cluster_distance(a, b)
    assert distance == 0


# different tips
@pytest.mark.parametrize(
    "b",
    (
        "(e, d, (a, b));",
        "(a, d, (a, b));",
        "(a, b, (c, (d, e)));",
    ),
)
def test_unrooted_rf_fails_different_tips(b):
    t1 = make_tree(treestring="(a, b, (c, d));")
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        unrooted_robinson_foulds(t1, t2)
    with pytest.raises(ValueError):
        unrooted_robinson_foulds(t2, t1)


@pytest.mark.parametrize(
    "b",
    (
        "(e, d, (a, b));",
        "(a, d, (a, b));",
        "(a, b, (c, (d, e)));",
    ),
)
def test_lrm_fails_different_tips(b):
    t1 = make_tree(treestring="(a, b, (c, d));")
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        lin_rajan_moret(t1, t2)
    with pytest.raises(ValueError):
        lin_rajan_moret(t2, t1)


@pytest.mark.parametrize(
    "b",
    (
        "((e, d), (a, b));",
        "((a, d), (a, b));",
        "(a, (b, (c, (d, e))));",
    ),
)
def test_rooted_rf_fails_different_tips(b):
    t1 = make_tree(treestring="((a, b), (c, d));")
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        rooted_robinson_foulds(t1, t2)
    with pytest.raises(ValueError):
        rooted_robinson_foulds(t2, t1)


@pytest.mark.parametrize(
    "b",
    (
        "((e, d), (a, b));",
        "((a, d), (a, b));",
        "(a, (b, (c, (d, e))));",
    ),
)
def test_matching_cluster_fails_different_tips(b):
    t1 = make_tree(treestring="((a, b), (c, d));")
    t2 = make_tree(treestring=b)
    with pytest.raises(ValueError):
        matching_cluster_distance(t1, t2)
    with pytest.raises(ValueError):
        matching_cluster_distance(t2, t1)


# unrooted distance measures applied to examples with rooted trees
@pytest.mark.parametrize(
    "a,b",
    (
        ("((a, b), (c, d));", "(d, c, (a, b));"),
        ("(a, b, (c, d));", "((d, c), (a, b));"),
        ("((a, b), (c, d));", "((a, c), (b, d));"),
    ),
)
def test_unrooted_dist_fails_rooted_trees(a, b):
    t1 = make_tree(treestring=a)
    t2 = make_tree(treestring=b)

    for unrooted_distance_measure in UNROOTED_DISTANCE_MEASURES:
        with pytest.raises(ValueError):
            unrooted_distance_measure(t1, t2)


# rooted distance measures applied to examples with unrooted trees
@pytest.mark.parametrize(
    "a,b",
    (
        ("((a, b), (c, d));", "(d, c, (a, b));"),
        ("(a, b, (c, d));", "((d, c), (a, b));"),
        ("(a, b, (c, d));", "(d, c, (a, b));"),
    ),
)
def test_rooted_dist_fails_unrooted_trees(a, b):
    t1 = make_tree(treestring=a)
    t2 = make_tree(treestring=b)

    for rooted_distance_measure in ROOTED_DISTANCE_MEASURES:
        with pytest.raises(ValueError):
            rooted_distance_measure(t1, t2)

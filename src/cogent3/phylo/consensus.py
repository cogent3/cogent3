#! /usr/bin/env python
"""This module implements methods for generating consensus trees from a list of trees"""
import warnings

from collections import defaultdict
from itertools import product

from cogent3 import make_tree
from cogent3.core.tree import TreeBuilder
from cogent3.util.misc import extend_docstring_from


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"


def majority_rule(trees, strict=False):
    """Determines the consensus tree from a list of rooted trees using the
     majority rules method of Margush and McMorris 1981

    Parameters
    ----------
    trees
        A list of cogent3.evolve.tree objects
    strict
        A boolean flag for strict majority rule tree
        construction when true only nodes occurring >50% will be used
        when false the highest scoring < 50% node will be used if
        there is more than one node with the same score this will be
        arbitrarily chosen on sort order

    Returns:
        a list of cogent3.evolve.tree objects
    """
    trees = [(1, tree) for tree in trees]
    return weighted_majority_rule(trees, strict, "count", method="rooted")


def weighted_majority_rule(
    weighted_trees, strict=False, attr="support", method="unrooted"
):
    """Calculate a greedy consensus tree in the sense of Bryant (2003), if
    weights are taken as counts. Branch lengths calculated as per Holland
    (2006).

    Parameters
    ----------
    weighted_trees : list
        A reverse ordered list of (weight, tree) tuples.
    strict : bool
        Discard splits or clusters with consensus weight <= 0.5.
    attr : str
        Edge parameter in which to store consensus weight.
    method : str
        'unrooted' or 'rooted': treat the trees as if they were such.

    Returns
    -------
    A list of consensus trees. List length will always be one if method is
    'unrooted'.

    Citations
    ---------
    Bryant, D. (2003). A classification of consensus methods for phylogenetics.
    DIMACS series in discrete mathematics and theoretical computer science,
    61:163-184.

    Holland, B. R., Jermiin, L. S., Moulton, V., and Investigators,
    S. T.-N. Y. (2006). Proceedings of the SMBE Tri-National Young
    Investigators' Workshop 2005. Improved consensus network techniques for
    genome-scale phylogeny. Mol Biol Evol, 23(5), 848-855.
    doi:10.1093/molbev/msj061
    """
    if method == "rooted":
        return weighted_rooted_majority_rule(weighted_trees, strict, attr)
    elif method == "unrooted":
        return weighted_unrooted_majority_rule(weighted_trees, strict, attr)
    else:
        raise ValueError('method must be "rooted" or "unrooted"')


@extend_docstring_from(weighted_majority_rule)
def weighted_rooted_majority_rule(weighted_trees, strict=False, attr="support"):
    cladecounts = {}
    edgelengths = {}
    total = 0
    for (weight, tree) in weighted_trees:
        total += weight
        edges = tree.get_edge_vector()
        for edge in edges:
            tips = edge.get_tip_names(includeself=True)
            tips = frozenset(tips)
            if tips not in cladecounts:
                cladecounts[tips] = 0
            cladecounts[tips] += weight
            length = edge.length and edge.length * weight
            if edgelengths.get(tips, None):
                edgelengths[tips] += length
            else:
                edgelengths[tips] = length
    cladecounts = [(count, clade) for (clade, count) in list(cladecounts.items())]
    cladecounts.sort()
    cladecounts.reverse()

    if strict:
        # Remove any with support < 50%
        for index, (count, clade) in enumerate(cladecounts):
            if count <= 0.5 * total:
                cladecounts = cladecounts[:index]
                break

    # Remove conflicts
    accepted_clades = set()
    counts = {}
    for (count, clade) in cladecounts:
        for accepted_clade in accepted_clades:
            if clade.intersection(accepted_clade) and not (
                clade.issubset(accepted_clade) or clade.issuperset(accepted_clade)
            ):
                break
        else:
            accepted_clades.add(clade)
            counts[clade] = count
            weighted_length = edgelengths[clade]
            edgelengths[clade] = weighted_length and weighted_length / count

    nodes = {}
    queue = []
    tree_build = TreeBuilder().create_edge
    for clade in accepted_clades:
        if len(clade) == 1:
            tip_name = next(iter(clade))
            params = {"length": edgelengths[clade], attr: counts[clade]}
            nodes[tip_name] = tree_build([], tip_name, params)
        else:
            queue.append(((len(clade), clade)))

    while queue:
        queue.sort()
        (size, clade) = queue.pop(0)
        new_queue = []
        for (size2, ancestor) in queue:
            if clade.issubset(ancestor):
                new_ancestor = (ancestor - clade) | frozenset([clade])
                counts[new_ancestor] = counts.pop(ancestor)
                edgelengths[new_ancestor] = edgelengths.pop(ancestor)
                ancestor = new_ancestor
            new_queue.append((len(ancestor), ancestor))
        children = [nodes.pop(c) for c in clade]
        assert len([children])
        nodes[clade] = tree_build(
            children, None, {attr: counts[clade], "length": edgelengths[clade]}
        )
        queue = new_queue

    for root in list(nodes.values()):
        root.name = "root"  # Yuk

    return [root for root in list(nodes.values())]


@extend_docstring_from(weighted_majority_rule)
def weighted_unrooted_majority_rule(weighted_trees, strict=False, attr="support"):
    # Calculate raw split lengths and weights
    split_weights = defaultdict(float)
    split_lengths = defaultdict(float)
    tips = None
    for (weight, tree) in weighted_trees:
        for split, params in list(get_splits(tree).items()):
            split_weights[split] += weight
            if params["length"] is None:
                split_lengths[split] = None
            else:
                split_lengths[split] += weight * params["length"]
        # Check that all trees have the same taxa
        if tips is None:
            tips = frozenset.union(*split)
        elif tips != frozenset.union(*split):
            raise NotImplementedError("all trees must have the same taxa")

    # Normalise split lengths by split weight and split weights by total weight
    for split in split_lengths:
        if not split_lengths[split] is None:
            split_lengths[split] /= split_weights[split]
    total_weight = sum(w for w, t in weighted_trees[::-1])
    weighted_splits = [(w / total_weight, s) for s, w in list(split_weights.items())]
    weighted_splits.sort(reverse=True)

    # Remove conflicts and any with support < 50% if strict
    accepted_splits = {}
    for weight, split in weighted_splits:
        if strict and weight <= 0.5:
            break

        for accepted_split in accepted_splits:
            for s, a in product(split, accepted_split):
                if s.isdisjoint(a):
                    break
            else:
                break
        else:
            accepted_splits[split] = {attr: weight, "length": split_lengths[split]}

    return [get_tree(accepted_splits)]


def get_splits(tree):
    """Return a dict keyed by the splits equivalent to the tree.
    Values are {'length' : edge.length} for the corresponding edge.
    """
    if len(tree.children) < 3:
        warnings.warn("tree is rooted - will return splits for unrooted tree")

    def getTipsAndSplits(tree):
        if tree.is_tip():
            return ({frozenset([tree.name]): {"length": tree.length}}, [tree.name])

        splits = defaultdict(lambda: {"length": 0.0})
        tips = []
        for child in tree.children:
            s, t = getTipsAndSplits(child)
            splits.update(s)
            tips.extend(t)
        if not tree.is_root():
            split = frozenset([t for s in splits for t in s])
            if tree.length is None:
                splits[split] = {"length": None}
            else:
                splits[split] = {"length": tree.length + splits[split]["length"]}
        return splits, tips

    splits, tips = getTipsAndSplits(tree)
    tips = frozenset(tips)
    return {frozenset([tips - s, s]): params for s, params in list(splits.items())}


def get_tree(splits):
    """Convert a dict keyed by splits into the equivalent tree.
    The dict values should be dicts appropriate for the params input to
    TreeBuilder.create_edge.
    """
    Edge = TreeBuilder().create_edge

    # Create a star from the tips
    tips = []
    the_rest = []
    for split, params in list(splits.items()):
        small, big = sorted(split, key=len)
        if len(small) == 1:
            for name in small:
                tip = Edge(None, name, params)
            tip.Split = small
            tips.append(tip)
        else:
            the_rest.append((split, params))
    tree = Edge(tips, "root", {})

    # Add the rest of the splits, one by one
    def addHalfSplit(edge, half, params):
        included = []
        test_half = frozenset([])
        for child in edge.children:
            if child.Split > half:  # This is not the droid you are looking for
                return addHalfSplit(child, half, params)
            if child.Split <= half:
                included.append(child)
                test_half = test_half.union(child.Split)

        if test_half == half:  # Found it
            split = Edge(included, None, params)
            split.Split = half
            for moved in included:
                edge.remove_node(moved)
            edge.append(split)
            return True

        return False

    for split, params in the_rest:
        for half in split:
            if addHalfSplit(tree, half, params):
                break

    # Balance the tree for the sake of reproducibility
    tree = tree.balanced()
    return tree


if __name__ == "__main__":
    import sys

    trees = []
    for filename in sys.argv[1:]:
        for tree in open(filename):
            trees.append(make_tree(treestring=tree))
    print(f"Consensus of {len(trees)} trees from {sys.argv[1:]}")
    outtrees = majority_rule(trees, strict=True)
    for tree in outtrees:
        print(tree.ascii_art(compact=True, show_internal=False))

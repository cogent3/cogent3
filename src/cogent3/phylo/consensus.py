#! /usr/bin/env python
"""This module implements methods for generating consensus trees from a list of trees"""

import warnings
from collections import defaultdict
from collections.abc import Iterable
from itertools import product
from typing import Literal, TypeAlias, cast

from cogent3.core.tree import PhyloNode, TreeBuilder
from cogent3.util.misc import extend_docstring_from


def majority_rule(trees: Iterable[PhyloNode], strict: bool = False) -> list[PhyloNode]:
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
        a list of PhyloNode objects
    """
    weighted_trees = [(1, tree) for tree in trees]
    return weighted_majority_rule(weighted_trees, strict, "count", method="rooted")


def weighted_majority_rule(
    weighted_trees: Iterable[tuple[float, PhyloNode]],
    strict: bool = False,
    attr: str = "support",
    method: Literal["unrooted", "rooted"] = "unrooted",
) -> list[PhyloNode]:
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
    if method == "unrooted":
        return weighted_unrooted_majority_rule(weighted_trees, strict, attr)
    msg = 'method must be "rooted" or "unrooted"'
    raise ValueError(msg)


NestedFrozenset: TypeAlias = frozenset["str | NestedFrozenset"]


@extend_docstring_from(weighted_majority_rule)
def weighted_rooted_majority_rule(
    weighted_trees: Iterable[tuple[float, PhyloNode]],
    strict: bool = False,
    attr: str = "support",
) -> list[PhyloNode]:
    cladecounts_dict: dict[frozenset[str], float] = {}
    edgelengths: dict[NestedFrozenset, float | None] = {}
    total: float = 0
    for weight, tree in weighted_trees:
        total += weight
        edges = tree.get_edge_vector()
        for edge in edges:
            tips = frozenset(edge.get_tip_names(include_self=True))
            if tips not in cladecounts_dict:
                cladecounts_dict[tips] = 0
            cladecounts_dict[tips] += weight
            length = edge.length and edge.length * weight
            if edgelength := edgelengths.get(tips):
                edgelengths[tips] = edgelength + cast("float", length)
            else:
                edgelengths[tips] = length
    cladecounts: list[tuple[float, NestedFrozenset]] = [
        (count, clade) for (clade, count) in list(cladecounts_dict.items())
    ]
    cladecounts.sort()
    cladecounts.reverse()

    if strict:
        # Remove any with support < 50%
        for index, (count, _) in enumerate(cladecounts):
            if count <= 0.5 * total:
                cladecounts = cladecounts[:index]
                break

    # Remove conflicts
    accepted_clades: set[NestedFrozenset] = set()
    counts: dict[NestedFrozenset, float] = {}
    for count, clade in cladecounts:
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

    nodes: dict[str | NestedFrozenset, PhyloNode] = {}
    queue: list[tuple[int, NestedFrozenset]] = []
    tree_build = TreeBuilder().create_edge
    for clade in accepted_clades:
        if len(clade) == 1:
            tip_name = next(iter(clade))
            params: dict[str, float | None] = {
                "length": edgelengths[clade],
                attr: counts[clade],
            }
            nodes[tip_name] = tree_build([], cast("str", tip_name), params)
        else:
            queue.append((len(clade), clade))

    while queue:
        queue.sort()
        (_, clade) = queue.pop(0)
        new_queue: list[tuple[int, NestedFrozenset]] = []
        for _, ancestor in queue:
            if clade.issubset(ancestor):
                new_ancestor = (ancestor - clade) | frozenset([clade])
                counts[new_ancestor] = counts.pop(ancestor)
                edgelengths[new_ancestor] = edgelengths.pop(ancestor)
                ancestor = new_ancestor
            new_queue.append((len(ancestor), ancestor))
        children = [nodes.pop(c) for c in clade]

        nodes[clade] = tree_build(
            children,
            None,
            {attr: counts[clade], "length": edgelengths[clade]},
        )
        queue = new_queue

    for root in list(nodes.values()):
        root.name = "root"  # Yuk

    return list(nodes.values())


@extend_docstring_from(weighted_majority_rule)
def weighted_unrooted_majority_rule(
    weighted_trees: Iterable[tuple[float, PhyloNode]],
    strict: bool = False,
    attr: str = "support",
) -> list[PhyloNode]:
    # Calculate raw split lengths and weights
    split_weights: dict[frozenset[frozenset[str]], float] = defaultdict(float)
    split_lengths: dict[frozenset[frozenset[str]], float | None] = defaultdict(float)
    tips: frozenset[str] | None = None
    for weight, tree in weighted_trees:
        split: frozenset[frozenset[str]] = frozenset()
        for split, params in list(get_splits(tree).items()):
            split_weights[split] += weight
            if params["length"] is None:
                split_lengths[split] = None
            else:
                split_lengths[split] = (
                    cast("float", split_lengths[split]) + weight * params["length"]
                )
        # Check that all trees have the same taxa
        if tips is None:
            tips = frozenset[str].union(*split)
        elif tips != frozenset[str].union(*split):
            msg = "all trees must have the same taxa"
            raise NotImplementedError(msg)

    # Normalise split lengths by split weight and split weights by total weight
    for split, length in split_lengths.items():
        if length is not None:
            split_lengths[split] = length / split_weights[split]

    total_weight = sum(w for w, _ in weighted_trees)
    weighted_splits = [(w / total_weight, s) for s, w in list(split_weights.items())]
    weighted_splits.sort(reverse=True)

    # Remove conflicts and any with support < 50% if strict
    accepted_splits: dict[frozenset[frozenset[str]], dict[str, float | None]] = {}
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


def get_splits(
    tree: PhyloNode,
) -> dict[frozenset[frozenset[str]], dict[str, float | None]]:
    """Return a dict keyed by the splits equivalent to the tree.
    Values are {'length' : edge.length} for the corresponding edge.
    """
    if len(tree.children) < 3:
        warnings.warn(
            "tree is rooted - will return splits for unrooted tree",
            stacklevel=2,
        )

    def get_tips_and_splits(
        tree: PhyloNode,
    ) -> tuple[dict[frozenset[str], dict[str, float | None]], list[str]]:
        if tree.is_tip():
            return ({frozenset([tree.name]): {"length": tree.length}}, [tree.name])

        splits: dict[frozenset[str], dict[str, float | None]] = defaultdict(
            lambda: {"length": 0.0}
        )
        tips: list[str] = []
        for child in tree.children:
            s, t = get_tips_and_splits(child)
            splits.update(s)
            tips.extend(t)
        if not tree.is_root():
            split = frozenset([t for s in splits for t in s])
            if tree.length is None:
                splits[split] = {"length": None}
            else:
                split_length = cast("float", splits[split]["length"])
                splits[split] = {"length": tree.length + split_length}
        return splits, tips

    splits, tips = get_tips_and_splits(tree)
    return {
        frozenset([frozenset(tips) - s, s]): params
        for s, params in list(splits.items())
    }


def get_tree(
    splits: dict[frozenset[frozenset[str]], dict[str, float | None]],
) -> PhyloNode:
    """Convert a dict keyed by splits into the equivalent tree.
    The dict values should be dicts appropriate for the params input to
    TreeBuilder.create_edge.
    """
    edge_builder = TreeBuilder().create_edge

    # Create a star from the tips
    tips: list[PhyloNode] = []
    the_rest: list[tuple[frozenset[frozenset[str]], dict[str, float | None]]] = []
    for split, params in list(splits.items()):
        small, _ = sorted(split, key=len)
        if len(small) == 1:
            tip = edge_builder(None, next(iter(small)), params)
            tip.params["Split"] = small
            tips.append(tip)
        else:
            the_rest.append((split, params))
    tree = edge_builder(tips, "root", {})

    # Add the rest of the splits, one by one
    def add_half_split(
        edge: PhyloNode, half: frozenset[str], params: dict[str, float | None]
    ) -> bool:
        included: list[PhyloNode] = []
        test_half: frozenset[str] = frozenset([])
        for child in edge.children:
            if (
                child.params["Split"] > half
            ):  # This is not the droid you are looking for
                return add_half_split(child, half, params)
            if child.params["Split"] <= half:
                included.append(child)
                test_half = test_half.union(child.params["Split"])

        if test_half == half:  # Found it
            split = edge_builder(included, None, params)
            split.params["Split"] = half
            for moved in included:
                edge.remove_node(moved)
            edge.append(split)
            return True

        return False

    for split, params in the_rest:
        for half in split:
            if add_half_split(tree, half, params):
                break

    # Balance the tree for the sake of reproducibility
    return tree.balanced()

#!/usr/bin/env python
"""Generalised Neighbour Joining phylogenetic tree estimation.

By default negative branch lengths are reset to 0.0 during the calculations.

This is based on the algorithm of Studier and Keppler, as described in the book
Biological sequence analysis by Durbin et al

Generalised as described by Pearson, Robins & Zhang, 1999.
"""


from collections import deque

import numpy

from cogent3.core.tree import TreeBuilder
from cogent3.phylo.tree_collection import ScoredTreeCollection
from cogent3.phylo.util import distance_dict_to_2D
from cogent3.util import progress_display as UI


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class LightweightTreeTip(str):
    def convert(self, constructor, length):
        node = constructor([], str(self), {})
        node.length = max(0.0, length)
        return node


class LightweightTreeNode(frozenset):
    """Set of (length, child node) tuples"""

    def convert(self, constructor=None, length=None):
        if constructor is None:
            constructor = TreeBuilder().create_edge
        children = [child.convert(constructor, clength) for (clength, child) in self]
        node = constructor(children, None, {})
        if length is not None:
            node.length = max(0.0, length)
        return node

    def __or__(self, other):
        return type(self)(frozenset.__or__(self, other))


class PartialTree(object):
    """A candidate tree stored as
    (distance matrix, list of subtrees, list of tip sets, set of partitions, score).
    At each iteration (ie: call of the join method) the number of subtrees
    is reduced as 2 of them are joined, while the number of partitions is
    increased as a new edge is introduced.
    """

    def __init__(self, d, nodes, tips, score):
        self.d = d
        self.nodes = nodes
        self.tips = tips
        self.score = score

    def get_dist_saved_join_score_matrix(self):
        d = self.d
        L = len(d)
        r = numpy.sum(d, 0)
        Q = d - numpy.add.outer(r, r) / (L - 2.0)
        return Q / 2.0 + sum(r) / (L - 2.0) / 2 + self.score

    def join(self, i, j):
        tips = self.tips[:]
        new_tip_set = tips[i] | tips[j]
        nodes = self.nodes[:]
        d = self.d.copy()

        # Branch lengths from i and j to new node
        L = len(nodes)
        r = numpy.sum(d, axis=0)
        ij_dist_diff = (r[i] - r[j]) / (L - 2.0)
        left_length = 0.5 * (d[i, j] + ij_dist_diff)
        right_length = 0.5 * (d[i, j] - ij_dist_diff)

        score = self.score + d[i, j]

        left_length = max(0.0, left_length)
        right_length = max(0.0, right_length)

        # Join i and k to make new node
        new_node = LightweightTreeNode(
            [(left_length, nodes[i]), (right_length, nodes[j])]
        )

        # Store new node at i
        new_dists = 0.5 * (d[i] + d[j] - d[i, j])
        d[:, i] = new_dists
        d[i, :] = new_dists
        d[i, i] = 0.0
        nodes[i] = new_node
        tips[i] = new_tip_set

        # Eliminate j
        d[j, :] = d[L - 1, :]
        d[:, j] = d[:, L - 1]
        assert d[j, j] == 0.0, d
        d = d[0 : L - 1, 0 : L - 1]
        nodes[j] = nodes[L - 1]
        nodes.pop()
        tips[j] = tips[L - 1]
        tips.pop()

        return type(self)(d, nodes, tips, score)

    def asScoreTreeTuple(self):
        assert len(self.nodes) == 3  # otherwise next line needs generalizing
        lengths = numpy.sum(self.d, axis=0) - numpy.sum(self.d) / 4
        root = LightweightTreeNode(list(zip(lengths, self.nodes)))
        tree = root.convert()
        tree.name = "root"
        return (self.score + sum(lengths), tree)


class Pair(object):
    """A candidate neighbour join, not turned into an actual PartialTree until
    and unless we decide to use it, because calculating just the topology is
    faster than calculating the whole new distance matrix etc. as well."""

    __slots__ = ["tree", "i", "j", "topology", "new_partition"]

    def __init__(self, tree, i, j, topology, new_partition):
        self.tree = tree
        self.i = i
        self.j = j
        self.topology = topology
        self.new_partition = new_partition

    def joined(self):
        new_tree = self.tree.join(self.i, self.j)
        new_tree.topology = self.topology
        return new_tree


def uniq_neighbour_joins(trees, encode_partition):
    """Generate all joinable pairs from all trees, best first,
    filtering out any duplicates"""
    L = len(trees[0].nodes)
    scores = numpy.zeros([len(trees), L, L])
    for (k, tree) in enumerate(trees):
        scores[k] = tree.get_dist_saved_join_score_matrix()
    topologies = set()
    order = numpy.argsort(scores.flat)
    for index in order:
        (k, ij) = divmod(index, L * L)
        (i, j) = divmod(ij, L)
        if i == j:
            continue
        tree = trees[k]
        new_tip_set = tree.tips[i] | tree.tips[j]
        new_partition = encode_partition(new_tip_set)
        # check is new topology
        topology = tree.topology | frozenset([new_partition])
        if topology in topologies:
            continue
        yield Pair(tree, i, j, topology, new_partition)
        topologies.add(topology)


@UI.display_wrap
def gnj(dists, keep=None, dkeep=0, ui=None):
    """Arguments:
        - dists: dict of (name1, name2): distance
        - keep: number of best partial trees to keep at each iteration,
          and therefore to return.  Same as Q parameter in original GNJ paper.
        - dkeep: number of diverse partial trees to keep at each iteration,
          and therefore to return.  Same as D parameter in original GNJ paper.
    Result:
        - a sorted list of (tree length, tree) tuples
    """
    try:
        dists = dists.to_dict()
    except AttributeError:
        pass

    (names, d) = distance_dict_to_2D(dists)

    if keep is None:
        keep = len(names) * 5
    all_keep = keep + dkeep

    # For recognising duplicate topologies, encode partitions (ie: edges) as
    # frozensets of tip names, which should be quickly comparable.
    arbitrary_anchor = names[0]
    all_tips = frozenset(names)

    def encode_partition(tips):
        included = frozenset(tips)
        if arbitrary_anchor not in included:
            included = all_tips - included
        return included
        # could also convert to long int, or cache, would be faster?

    tips = [frozenset([n]) for n in names]
    nodes = [LightweightTreeTip(name) for name in names]
    star_tree = PartialTree(d, nodes, tips, 0.0)
    star_tree.topology = frozenset([])
    trees = [star_tree]

    # Progress display auxiliary code
    template = " size %%s/%s  trees %%%si" % (len(names), len(str(all_keep)))
    total_work = 0
    max_candidates = 1
    total_work_before = {}
    for L in range(len(names), 3, -1):
        total_work_before[L] = total_work
        max_candidates = min(all_keep, max_candidates * L * (L - 1) // 2)
        total_work += max_candidates

    def _show_progress():
        t = len(next_trees)
        work_done = total_work_before[L] + t
        ui.display(msg=template % (L, t), progress=work_done / total_work)

    for L in range(len(names), 3, -1):
        # Generator of candidate joins, best first.
        # Note that with dkeep>0 this generator is used up a bit at a time
        # by 2 different interupted 'for' loops below.
        candidates = uniq_neighbour_joins(trees, encode_partition)

        # First take up to 'keep' best ones
        next_trees = []
        _show_progress()
        for pair in candidates:
            next_trees.append(pair)
            if len(next_trees) == keep:
                break
        _show_progress()

        # The very best one is used as an anchor for measuring the
        # topological distance to others
        best_topology = next_trees[0].topology
        prior_td = [len(best_topology ^ tree.topology) for tree in trees]

        # Maintain a separate queue of joins for each possible
        # topological distance
        max_td = (max(prior_td) + 1) // 2
        queue = [deque() for g in range(max_td + 1)]
        queued = 0

        # Now take up to dkeep joins, an equal number of the best at each
        # topological distance, while not calculating any more TDs than
        # necessary.
        prior_td = dict(list(zip(list(map(id, trees)), prior_td)))
        target_td = 1
        while (candidates or queued) and len(next_trees) < all_keep:
            if candidates and not queue[target_td]:
                for pair in candidates:
                    diff = pair.new_partition not in best_topology
                    td = (prior_td[id(pair.tree)] + [-1, +1][diff]) // 2
                    # equiv, slower: td = len(best_topology ^ topology) // 2
                    queue[td].append(pair)
                    queued += 1
                    if td == target_td:
                        break
                else:
                    candidates = None
            if queue[target_td]:
                next_trees.append(queue[target_td].popleft())
                queued -= 1
                _show_progress()

            target_td = target_td % max_td + 1

        trees = [pair.joined() for pair in next_trees]

    result = [tree.asScoreTreeTuple() for tree in trees]
    result.sort()
    return ScoredTreeCollection(result)


def nj(dists, show_progress=True):
    """Arguments:
    - dists: dict of (name1, name2): distance
    """
    (result,) = gnj(dists, keep=1, show_progress=show_progress)
    (score, tree) = result
    return tree

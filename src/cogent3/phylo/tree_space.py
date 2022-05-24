#!/usr/bin/env python

import itertools

import numpy

from cogent3.core.tree import TreeBuilder
from cogent3.phylo.tree_collection import ScoredTreeCollection
from cogent3.util import checkpointing
from cogent3.util import progress_display as UI


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def ismallest(data, size):
    """There are many ways to get the k smallest items from an N sequence, and
    which one performs best depends on k, N and k/N.  This algorithm appears to
    beat anything heapq can do, and stays with a factor of 2 of sort() and
    min().  Is uses memory O(2*k) and so is particularly suitable for lazy
    application to large N.  It returns the smallest k sorted too."""
    limit = 2 * size
    data = iter(data)
    best = list(itertools.islice(data, limit))
    while True:
        best.sort()
        if len(best) <= size:
            break
        del best[size:]
        worst_of_best = best[-1]
        for item in data:
            if item < worst_of_best:
                best.append(item)
                if len(best) > limit:
                    break
    return best


# Trees are represented as "ancestry" matricies in which A[i,j] iff j is an
# ancestor of i.  For LS calculations the ancestry matrix is converted
# to a "paths" matrix or "split metric" in which S[p,j] iff the path between
# the pth pair of tips passes through edge j.  For ML calculations the
# ancestry matrix is converted back into an ordinary cogent tree object.


def tree2ancestry(tree, order=None):
    nodes = tree.unrooted().get_edge_vector()[:-1]
    if order is not None:
        lookup = dict([(k, i) for (i, k) in enumerate(order)])

        def _ordered_tips_first(n):
            if n.children:
                return len(order)
            else:
                return lookup[n.name]

        nodes.sort(key=_ordered_tips_first)

    n = len(nodes)
    A = numpy.zeros([n, n], int)
    seen = {}
    for (i, node) in enumerate(nodes):
        A[i, i] = 1
        seen[id(node)] = i
        for c in node.children:
            A[:, i] |= A[:, seen[id(c)]]
    names = [n.name for n in nodes if not n.children]
    lengths = [n.length for n in nodes]
    return (A, names, lengths)


def ancestry2tree(A, lengths, tip_names):
    """Convert edge x edge ancestry matrix to a cogent Tree object"""
    tips = {}
    tip = 0
    for i in range(len(A)):
        if numpy.sum(A[:, i]) == 1:
            tips[i] = tip_names[tip]
            tip += 1
    assert tip == len(tip_names)

    constructor = TreeBuilder().create_edge
    free = {}
    for i in numpy.argsort(numpy.sum(A, axis=0)):
        children = [j for j in range(len(A)) if A[j, i] and j != i]
        child_nodes = [free.pop(j) for j in children if j in free]
        if child_nodes:
            name = None
        else:
            name = tips[i]
        if lengths is None:
            params = {}
        else:
            params = {"length": lengths[i]}
        node = constructor(child_nodes, name, params)
        free[i] = node
    return constructor(list(free.values()), "root", {})


def grown(B, split_edge):
    """Ancestry matrix 'B' with one extra leaf added at 'split_edge'.
    Row/column order within the matrix is independent of the topology it
    represents. The added leaf will be the last one in the matrix, which keeps
    the leaf node order the same as the order in which they are added, which is
    what is assumed by ancestry2tree and ancestry2paths"""
    n = len(B)
    A = numpy.zeros([n + 2, n + 2], int)
    A[:n, :n] = B
    (sibling, parent) = (n, n + 1)
    A[sibling] = A[parent] = A[split_edge]
    A[:, parent] = A[:, split_edge]
    A[sibling, split_edge] = 0
    A[parent, split_edge] = 0
    A[sibling, sibling] = 1
    A[parent, parent] = 1
    A[sibling, parent] = 1
    A[split_edge, parent] = 1
    return A


class TreeEvaluator(object):
    """Subclass must provide make_tree_scorer and result2output"""

    def results2output(self, results):
        return ScoredTreeCollection(results)

    def evaluate_topology(self, tree):
        """Optimal (score, tree) for the one topology 'tree'"""
        (ancestry, names, lengths) = tree2ancestry(tree)
        evaluate = self.make_tree_scorer(names)
        (err, lengths) = evaluate(ancestry)
        return self.result2output(err, ancestry, lengths, names)

    def evaluate_tree(self, tree):
        """score for 'tree' with lengths as-is"""
        (ancestry, names, lengths) = tree2ancestry(tree)
        evaluate = self.make_tree_scorer(names)
        (err, result) = evaluate(ancestry, lengths=lengths)
        return err

    def _consistentNameOrder(self, fixed_names, ordered_names=None):
        """fixed_names followed by ordered_names without duplicates"""
        all_names = set(self.names)

        fixed_names_set = set(fixed_names)
        assert fixed_names_set.issubset(all_names)

        if ordered_names:
            assert set(ordered_names).issubset(all_names)
        else:
            ordered_names = self.names
        return list(fixed_names) + [
            n for n in ordered_names if n not in fixed_names_set
        ]

    @UI.display_wrap
    def trex(
        self,
        a=8,
        k=1000,
        start=None,
        order=None,
        return_all=False,
        filename=None,
        interval=None,
        show_progress=False,
        ui=None,
    ):
        """TrexML policy for tree sampling - all trees up to size 'a' and
        then keep no more than 'k' best trees at each tree size.
        'order' is an optional list of tip names.
        'start' is an optional list of initial trees.  Each of the trees must
        contain the same tips.
        'filename' and 'interval' control checkpointing.

        Advanced step-wise addition algorithm
        M. J. Wolf, S. Easteal, M. Kahn, B. D. McKay, and L. S. Jermiin.
        Trexml: a maximum-likelihood approach for extensive tree-space
        exploration.
        Bioinformatics, 16(4):383 94, 2000."""

        checkpointer = checkpointing.Checkpointer(filename, interval)
        if checkpointer.available():
            (init_tree_size, fixed_names, trees) = checkpointer.load()
            names = self._consistentNameOrder(fixed_names, order)
        elif start is not None:
            if not isinstance(start, list):
                start = [start]
            fixed_names = start[0].get_tip_names()
            names = self._consistentNameOrder(fixed_names, order)
            trees = []
            for tree in start:
                # check the start tree represents a subset of tips
                assert set(tree.get_tip_names()) < set(
                    self.names
                ), "Starting tree names not a subset of the sequence names"

                (ancestry, fixed_names2, lengths) = tree2ancestry(
                    tree, order=fixed_names
                )
                assert fixed_names2 == fixed_names
                trees.append((None, None, ancestry))
            init_tree_size = len(fixed_names)
        else:
            trees = [(None, None, numpy.identity(3, int))]
            names = self._consistentNameOrder([], order)
            init_tree_size = 3

        tree_size = len(names)
        assert tree_size > 3
        if a > tree_size:
            a = tree_size
        if a < 4:
            a = 4

        # All trees of size a-1, no need to compare them
        for n in range(init_tree_size + 1, a):
            trees2 = []
            for (err2, lengths2, ancestry) in trees:
                for split_edge in range(len(ancestry)):
                    ancestry2 = grown(ancestry, split_edge)
                    trees2.append((None, None, ancestry2))
            trees = trees2
            init_tree_size = n

        # Pre calculate how much work is to be done, for progress display
        tree_count = len(trees)
        total_work = 0
        work_done = [0] * (init_tree_size + 1)
        for n in range(init_tree_size + 1, tree_size + 1):
            evals = tree_count * (n * 2 - 5)
            total_work += evals * n
            tree_count = min(k, evals)
            work_done.append(total_work)

        # For each tree size, grow at each edge of each tree. Keep best k.
        for n in range(init_tree_size + 1, tree_size + 1):
            evaluate = self.make_tree_scorer(names[:n])

            def grown_tree(spec):
                (tree_ordinal, tree, split_edge) = spec
                (old_err, old_lengths, old_ancestry) = tree
                ancestry = grown(old_ancestry, split_edge)
                (err, lengths) = evaluate(ancestry)
                return (err, tree_ordinal, split_edge, lengths, ancestry)

            specs = [
                (i, tree, edge)
                for (i, tree) in enumerate(trees)
                for edge in range(n * 2 - 5)
            ]

            candidates = ui.imap(
                grown_tree,
                specs,
                noun=f"{n} leaf tree",
                start=work_done[n - 1] / total_work,
                end=work_done[n] / total_work,
            )

            best = ismallest(candidates, k)

            trees = [
                (err, lengths, ancestry)
                for (err, parent_ordinal, split_edge, lengths, ancestry) in best
            ]

            checkpointer.record((n, names[:n], trees))

        results = (
            self.result2output(err, ancestry, lengths, names)
            for (err, lengths, ancestry) in trees
        )
        if return_all:
            result = self.results2output(results)
        else:
            result = next(results)
        return result

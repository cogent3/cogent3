#!/usr/bin/env python
import numpy

from numpy.linalg import solve as solve_linear_equations

from .tree_space import TreeEvaluator, ancestry2tree
from .util import (
    distance_dict_and_names_to_1D,
    distance_dict_to_1D,
    triangular_order,
)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# This is a fairly slow implementation and NOT suitable for large trees.

# Trees are represented as "ancestry" matricies in which A[i,j] iff j is an
# ancestor of i.  For the LS calculations the ancestry matrix is converted
# to a "paths" matrix or "split metric" in which S[p,j] iff the path between
# the pth pair of tips passes through edge j.


def _ancestry2paths(A):
    """Convert edge x edge ancestry matrix to tip-to-tip path x edge
    split metric matrix.  The paths will be in the same triangular matrix order
    as produced by distance_dict_and_names_to_1D, provided that the tips appear in
    the correct order in A"""
    tips = [i for i in range(A.shape[0]) if sum(A[:, i]) == 1]
    paths = []
    for (tip1, tip2) in triangular_order(tips):
        path = A[tip1] ^ A[tip2]
        paths.append(path)
    return numpy.array(paths)


class WLS(TreeEvaluator):
    """(err, best_tree) = WLS(dists).trex()"""

    def __init__(self, dists, weights=None):
        """Arguments:
        - dists: a dict with structure (seq1, seq2): distance
        - weights: an equivalently structured dict with measurements of
          variability of the distance estimates. By default, the sqrt of
          distance is used."""
        try:
            dists = dists.to_dict()
        except AttributeError:
            pass
        self.dists = dists
        self.weights = weights or dict(
            (key, 1.0 / (self.dists[key] ** 2)) for key in self.dists
        )
        (self.names, dists) = distance_dict_to_1D(self.dists)

    def make_tree_scorer(self, names):
        dists = distance_dict_and_names_to_1D(self.dists, names)
        weights = distance_dict_and_names_to_1D(self.weights, names)
        # dists and weights are 1D forms of triangular tip x tip matrices
        # The order of the tip-to-tip paths is the same for dists, weights and
        # A
        weights_dists = weights * dists

        def evaluate(
            ancestry,
            lengths=None,
            sum=sum,
            _ancestry2paths=_ancestry2paths,
            dot=numpy.dot,
            maximum=numpy.maximum,
            transpose=numpy.transpose,
            solve=solve_linear_equations,
        ):
            A = _ancestry2paths(ancestry)
            if lengths is None:
                At = transpose(A)
                X = dot(weights * At, A)
                y = dot(At, weights_dists)
                lengths = solve(X, y)
                lengths = maximum(lengths, 0.0)
            diffs = dot(A, lengths) - dists
            err = sum(diffs ** 2)
            return (err, lengths)

        return evaluate

    def result2output(self, err, ancestry, lengths, names):
        return (err, ancestry2tree(ancestry, lengths, names))


def wls(*args, **kw):
    (err, tree) = WLS(*args).trex(**kw)
    return tree

from itertools import product
from typing import List

import numpy

from scipy.optimize import linear_sum_assignment

from cogent3.core.tree import PhyloNode


def lin_rajan_moret(tree1: PhyloNode, tree2: PhyloNode) -> float:
    """calculate the lin-rajan-moret distance (matchign distance) between two trees

    trees should have matching tips and must not be rooted.

    Parameters
    ----------
    tree1: PhyloNode
    tree2: PhyloNode
        trees to calculate distance between

    Returns
    -------
    float 
        the lin-rajan-moret distance

    Notes
    -----
    see:
    Lin et al. 2011
    A Metric for Phylogenetic Trees Based on Matching

    """
    names = tree1.get_tip_names()

    if set(names) != set(tree2.get_tip_names()):
        raise ValueError("tree tip names must match")
    if len(tree1.children) == 2 or len(tree2.children) == 2:
        raise ValueError("trees must be unrooted")

    names.sort()

    vector1 = _convert_tree_to_vectors(tree1, names)
    vector2 = _convert_tree_to_vectors(tree2, names)

    matching_distance = _matched_distance(vector1, vector2)

    return float(matching_distance)


def _convert_tree_to_vectors(tree: PhyloNode, tip_names: List) -> numpy.ndarray:
    ref_tip = tip_names[0]
    name_set = set(tip_names)
    name_index = {n: i for i, n in enumerate(tip_names)}
    # we get the tree as a set of splits.
    splits = tree.subsets()
    rows = []
    for split in splits:
        row = numpy.zeros(len(tip_names), int)
        # Cogent only returns one side of a
        # split, so we build the other side
        if ref_tip in split:
            names = list(split)
        else:
            names = list(name_set - split)
        indices = [name_index[n] for n in names]
        row[indices] = True
        rows.append(row)
    rows = numpy.array(rows)
    return rows


def _hamming(vector1: numpy.ndarray, vector2: numpy.ndarray) -> int:
    return numpy.count_nonzero(vector1 != vector2)


def _weight(vector1: numpy.ndarray, vector2: numpy.ndarray) -> int:
    return min(_hamming(vector1, vector2), _hamming(vector1, 1 - vector2))


def _bipartite_graph(vector1: numpy.ndarray, vector2: numpy.ndarray) -> numpy.ndarray:
    if not len(vector1) == len(vector2):
        raise ValueError("number of edges must be equal")

    B = numpy.empty([len(vector1)] * 2, int)
    for i, j in product(*[range(len(vector1))] * 2):
        B[i, j] = _weight(vector1[i], vector2[j])
    return B


def _matched_distance(vector1: numpy.ndarray, vector2: numpy.ndarray) -> int:
    B = _bipartite_graph(vector1, vector2)
    matching = linear_sum_assignment(B)
    return B[matching].sum()


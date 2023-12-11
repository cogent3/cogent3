from itertools import product
from typing import TYPE_CHECKING

import numpy as np

from scipy.optimize import linear_sum_assignment


if TYPE_CHECKING:
    from cogent3.core.tree import TreeNode


def get_tree_distance_measure(method: str, is_rooted: bool):
    if method in _TREE_DISTANCE_FUNCTIONS:
        return _TREE_DISTANCE_FUNCTIONS[method]

    if is_rooted and method in _ROOTED_TREE_DISTANCE_FUNCTIONS:
        return _ROOTED_TREE_DISTANCE_FUNCTIONS[method]

    if not is_rooted and method in _UNROOTED_TREE_DISTANCE_FUNCTIONS:
        return _UNROOTED_TREE_DISTANCE_FUNCTIONS[method]

    raise ValueError(
        f"Tree distance method '{method}' is not supported for {'' if is_rooted else 'un'}rooted trees."
    )


def unrooted_robinson_foulds(tree1: "TreeNode", tree2: "TreeNode") -> int:
    """Calculate the Robinson-Foulds distance between two unrooted trees.

    Parameters
    ----------
    tree1, tree2: TreeNode
        Trees to calculate the distance between.

    Returns
    -------
    int
        The unrooted Robinson-Foulds distance.

    Notes
    -----
    For unrooted trees, the Robinson-Foulds distance [1]_ is defined
    as the cardinality of the symmetric difference of the set of splits
    for the two trees.

    Trees should have matching tips and must not be rooted.

    References
    ----------
    .. [1] Robinson, David F., and Leslie R. Foulds.
       Comparison of phylogenetic trees.
       Mathematical biosciences 53.1-2 (1981): 131-147.
    """

    names = tree1.get_tip_names()

    if set(names) != set(tree2.get_tip_names()):
        raise ValueError("tree tip names must match")
    if len(tree1.children) == 2 or len(tree2.children) == 2:
        raise ValueError("trees must be unrooted")

    tree1_clusters = tree1.subsets()
    tree2_clusters = tree2.subsets()

    tree1_splits = _compute_splits(tree1_clusters, names)
    tree2_splits = _compute_splits(tree2_clusters, names)

    return len(tree1_splits.symmetric_difference(tree2_splits))


def lin_rajan_moret(tree1: "TreeNode", tree2: "TreeNode") -> int:
    """Calculate the Lin-Rajan-Moret distance (matching distance)
    between two unrooted trees.

    Parameters
    ----------
    tree1, tree2: TreeNode
        Trees to calculate the distance between.

    Returns
    -------
    int
        The Lin-Rajan-Moret distance.

    Notes
    -----
    The Lin-Rajan-Moret distance [1]_ displays superior statistical
    properties than the Robinson-Foulds distance [2]_
    on unrooted trees.

    Trees should have matching tips and must not be rooted.

    References
    ----------
    .. [1] Lin et al. 2012
       A Metric for Phylogenetic Trees Based on Matching
       IEEE/ACM Transactions on Computational Biology and Bioinformatics
       vol. 9, no. 4, pp. 1014-1022, July-Aug. 2012
    .. [2] Robinson, David F., and Leslie R. Foulds.
       Comparison of phylogenetic trees.
       Mathematical biosciences 53.1-2 (1981): 131-147.
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

    return matching_distance


def rooted_robinson_foulds(tree1: "TreeNode", tree2: "TreeNode") -> int:
    """Calculate the Robinson-Foulds distance between two rooted trees.

    Parameters
    ----------
    tree1, tree2: TreeNode
        Trees to calculate the distance between.

    Returns
    -------
    int
        The rooted Robinson-Foulds distance

    Notes
    -----
    For rooted trees, the Robinson-Foulds distance [1]_ is defined as the
    cardinality of the symmetric difference of the set of clades for
    the two trees.

    Trees should have matching tips and must be rooted.

    References
    ----------
    .. [1] Robinson, David F., and Leslie R. Foulds.
       Comparison of phylogenetic trees.
       Mathematical biosciences 53.1-2 (1981): 131-147.
    """
    if set(tree1.get_tip_names()) != set(tree2.get_tip_names()):
        raise ValueError("tree tip names must match")
    if len(tree1.children) != 2 or len(tree2.children) != 2:
        raise ValueError("trees must be rooted")

    tree1_clusters = tree1.subsets()
    tree2_clusters = tree2.subsets()

    return len(tree1_clusters.symmetric_difference(tree2_clusters))


def matching_cluster_distance(tree1: "TreeNode", tree2: "TreeNode") -> int:
    """Calculate the Matching Cluster distance between two rooted trees.

    Parameters
    ----------
    tree1, tree2: TreeNode
        Trees to calculate the distance between.

    Returns
    -------
    int
        The Matching Cluster distance.

    Notes
    -----
    The Matching Cluster distance [1]_ is a similar metric
    to the beta-distance [2]_.

    The Matching Cluster distance [1]_ displays superior statistical
    properties than the Robinson-Foulds distance [3]_ on rooted trees.

    Trees should have matching tips and must be rooted.

    References
    ----------
    .. [1] Bogdanowicz, D., & Giaro, K. (2013).
       On a matching distance between rooted phylogenetic trees.
       International Journal of Applied Mathematics and Computer Science, 23(3), 669-684.
    .. [2] Boorman, S. A., & Olivier, D. C. (1973).
       Metrics on spaces of finite trees.
       Journal of Mathematical Psychology, 10(1), 26-59.
    .. [3] Robinson, David F., and Leslie R. Foulds.
       Comparison of phylogenetic trees.
       Mathematical biosciences 53.1-2 (1981): 131-147.
    """

    if set(tree1.get_tip_names()) != set(tree2.get_tip_names()):
        raise ValueError("tree tip names must match")
    if len(tree1.children) != 2 or len(tree2.children) != 2:
        raise ValueError("trees must be rooted")

    tree1_clusters = list(tree1.subsets())
    tree2_clusters = list(tree2.subsets())

    while len(tree1_clusters) < len(tree2_clusters):
        tree1_clusters.append(set())
    while len(tree2_clusters) < len(tree1_clusters):
        tree2_clusters.append(set())

    adjacency = np.zeros(shape=(len(tree1_clusters), len(tree2_clusters)))

    for i, cluster_1 in enumerate(tree1_clusters):
        for j, cluster_2 in enumerate(tree2_clusters):
            adjacency[i, j] = len(cluster_1.symmetric_difference(cluster_2))

    row_ind, col_ind = linear_sum_assignment(adjacency)
    distance = int(adjacency[row_ind, col_ind].sum())

    return distance


def _convert_tree_to_vectors(tree: "TreeNode", tip_names: list) -> np.ndarray:
    ref_tip = tip_names[0]
    name_set = set(tip_names)
    name_index = {n: i for i, n in enumerate(tip_names)}
    # we get the tree as a set of splits.
    splits = tree.subsets()
    rows = np.zeros((len(splits), len(tip_names)), dtype=bool)
    for i, split in enumerate(splits):
        row = rows[i]
        # Cogent only returns one side of a
        # split, so we build the other side
        if ref_tip in split:
            names = list(split)
        else:
            names = list(name_set - split)
        indices = [name_index[n] for n in names]
        row[indices] = True
    return rows


def _hamming(vector1: np.ndarray, vector2: np.ndarray) -> int:
    return (vector1 != vector2).sum()


def _weight(vector1: np.ndarray, vector2: np.ndarray) -> int:
    return min(_hamming(vector1, vector2), _hamming(vector1, 1 - vector2))


def _bipartite_graph(vector1: np.ndarray, vector2: np.ndarray) -> np.ndarray:
    if not len(vector1) == len(vector2):
        raise ValueError("number of edges must be equal")

    B = np.empty([len(vector1)] * 2, int)
    for i, j in product(*[range(len(vector1))] * 2):
        B[i, j] = _weight(vector1[i], vector2[j])
    return B


def _matched_distance(vector1: np.ndarray, vector2: np.ndarray) -> int:
    B = _bipartite_graph(vector1, vector2)
    matching = linear_sum_assignment(B)
    return B[matching].sum()


def _compute_splits(clades: frozenset[frozenset], tip_names: list) -> set[frozenset]:
    ref_tip = tip_names[0]
    name_set = frozenset(tip_names)

    splits = set()
    for clade in clades:
        if ref_tip in clade:
            splits.add(clade)
        else:
            splits.add(name_set - clade)
    return splits


_TREE_DISTANCE_FUNCTIONS = {
    "rooted_robinson_foulds": rooted_robinson_foulds,
    "unrooted_robinson_foulds": unrooted_robinson_foulds,
    "matching_cluster": matching_cluster_distance,
    "lin_rajan_moret": lin_rajan_moret,
    "rrf": rooted_robinson_foulds,
    "urf": unrooted_robinson_foulds,
    "mc": matching_cluster_distance,
    "lrm": lin_rajan_moret,
}

_ROOTED_TREE_DISTANCE_FUNCTIONS = {
    "rf": rooted_robinson_foulds,
    "matching": matching_cluster_distance,
}

_UNROOTED_TREE_DISTANCE_FUNCTIONS = {
    "rf": unrooted_robinson_foulds,
    "matching": lin_rajan_moret,
}

# usr/bin/env python
"""Functions to cluster using UPGMA

upgma takes an dictionary of pair tuples mapped to distances as input.

UPGMA_cluster takes an array and a list of PhyloNode objects corresponding
to the array as input. Can also generate this type of input from a DictArray using
inputs_from_dict_array function.

Both return a PhyloNode object of the UPGMA cluster
"""

import numpy

from numpy import argmin, average, diag, ravel, take

from cogent3.core.tree import PhyloNode
from cogent3.util.dict_array import DictArray


__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Catherine Lozuopone", "Rob Knight", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Production"


numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)
BIG_NUM = 1e305


def upgma(pairwise_distances):
    """Uses the UPGMA algorithm to cluster sequences

    pairwise_distances: a dictionary with pair tuples mapped to a distance
    returns a PhyloNode object of the UPGMA cluster
    """
    darr = DictArray(pairwise_distances)
    matrix_a, node_order = inputs_from_dict_array(darr)
    tree = UPGMA_cluster(matrix_a, node_order, BIG_NUM)
    index = 0
    for node in tree.traverse():
        if not node.parent:
            node.name = "root"
        elif not node.name:
            node.name = "edge." + str(index)
            index += 1
    return tree


def find_smallest_index(matrix):
    """returns the index of the smallest element in a numpy array

    for UPGMA clustering elements on the diagonal should first be
    substituted with a very large number so that they are always
    larger than the rest if the values in the array."""
    # get the shape of the array as a tuple (e.g. (3,3))
    shape = matrix.shape
    # turn into a 1 by x array and get the index of the lowest number
    matrix1D = ravel(matrix)
    lowest_index = argmin(matrix1D)
    # convert the lowest_index derived from matrix1D to one for the original
    # square matrix and return
    row_len = shape[0]
    return divmod(lowest_index, row_len)


def condense_matrix(matrix, smallest_index, large_value):
    """converges the rows and columns indicated by smallest_index

    Smallest index is returned from find_smallest_index.
    For both the rows and columns, the values for the two indices are
    averaged. The resulting vector replaces the first index in the array
    and the second index is replaced by an array with large numbers so that
    it is never chosen again with find_smallest_index.
    """
    first_index, second_index = smallest_index
    # get the rows and make a new vector that has their average
    rows = take(matrix, smallest_index, 0)
    new_vector = average(rows, 0)
    # replace info in the row and column for first index with new_vector
    matrix[first_index] = new_vector
    matrix[:, first_index] = new_vector
    # replace the info in the row and column for the second index with
    # high numbers so that it is ignored
    matrix[second_index] = large_value
    matrix[:, second_index] = large_value
    return matrix


def condense_node_order(matrix, smallest_index, node_order):
    """condenses two nodes in node_order based on smallest_index info

    This function is used to create a tree while condensing a matrix
    with the condense_matrix function. The smallest_index is retrieved
    with find_smallest_index. The first index is replaced with a node object
    that combines the two nodes corresponding to the indices in node order.
    The second index in smallest_index is replaced with None.
    Also sets the branch length of the nodes to 1/2 of the distance between
    the nodes in the matrix"""
    index1, index2 = smallest_index
    node1 = node_order[index1]
    node2 = node_order[index2]
    # get the distance between the nodes and assign 1/2 the distance to the
    # lengthproperty of each node
    distance = matrix[index1, index2]
    nodes = [node1, node2]
    d = distance / 2.0
    for n in nodes:
        if n.children:
            n.length = d - n.children[0].TipLength
        else:
            n.length = d
        n.TipLength = d
    # combine the two nodes into a new PhyloNode object
    new_node = PhyloNode()
    new_node.children.append(node1)
    new_node.children.append(node2)
    node1.parent = new_node
    node2.parent = new_node
    # replace the object at index1 with the combined node
    node_order[index1] = new_node
    # replace the object at index2 with None
    node_order[index2] = None
    return node_order


def UPGMA_cluster(matrix, node_order, large_number):
    """cluster with UPGMA

    matrix is a numpy array.
    node_order is a list of PhyloNode objects corresponding to the matrix.
    large_number will be assigned to the matrix during the process and
    should be much larger than any value already in the matrix.

    WARNING: Changes matrix in-place.
    WARNING: Expects matrix to already have diagonals assigned to large_number
             before this function is called.
    """
    num_entries = len(node_order)
    tree = None
    for i in range(num_entries - 1):
        smallest_index = find_smallest_index(matrix)
        index1, index2 = smallest_index
        # if smallest_index is on the diagonal set the diagonal to large_number
        if index1 == index2:
            matrix[diag([True] * len(matrix))] = large_number
            smallest_index = find_smallest_index(matrix)
        _ = condense_node_order(matrix, smallest_index, node_order)
        matrix = condense_matrix(matrix, smallest_index, large_number)
        tree = node_order[smallest_index[0]]
    return tree


def inputs_from_dict_array(darr):
    """makes inputs for UPGMA_cluster from a DictArray object"""
    darr.array += numpy.eye(darr.shape[0]) * BIG_NUM
    nodes = list(map(PhyloNode, darr.keys()))
    return darr.array, nodes

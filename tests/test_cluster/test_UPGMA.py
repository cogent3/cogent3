#!/usr/bin/env python
from unittest import TestCase, main

import numpy

from numpy import array
from numpy.testing import assert_allclose

from cogent3 import make_tree
from cogent3.cluster.UPGMA import (
    UPGMA_cluster,
    condense_matrix,
    condense_node_order,
    find_smallest_index,
    inputs_from_dict_array,
    upgma,
)
from cogent3.core.tree import PhyloNode
from cogent3.util.dict_array import DictArray


Float = numpy.core.numerictypes.sctype2char(float)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class UPGMATests(TestCase):
    """test the functions to cluster using UPGMA using numpy"""

    def setUp(self):
        """creates inputs"""
        self.pairwise_distances = {
            ("a", "b"): 1.0,
            ("a", "c"): 4.0,
            ("a", "d"): 20.0,
            ("a", "e"): 22.0,
            ("b", "c"): 5.0,
            ("b", "d"): 21.0,
            ("b", "e"): 23.0,
            ("c", "d"): 10.0,
            ("c", "e"): 12.0,
            ("d", "e"): 2.0,
        }
        # create a list of PhyloNode objects
        a, b, c, d, e = list(map(PhyloNode, "abcde"))
        self.node_order = [a, b, c, d, e]
        # create a numpy matrix object to cluster
        self.matrix = array(
            (
                [9999999, 1, 4, 20, 22],
                [1, 9999999, 5, 21, 23],
                [4, 5, 9999999, 10, 12],
                [20, 21, 10, 9999999, 2],
                [22, 23, 12, 2, 9999999],
            ),
            Float,
        )
        # create a numpy matrix with zero diagonals to test diagonal mask
        self.matrix_zeros = array(
            (
                [0, 1, 4, 20, 22],
                [1, 0, 5, 21, 23],
                [4, 5, 0, 10, 12],
                [20, 21, 10, 0, 2],
                [22, 23, 12, 2, 0],
            ),
            Float,
        )

        # create a numpy matrix with zero diagonals to test diagonal mask
        self.matrix_five = array(
            (
                [5, 1, 4, 20, 22],
                [1, 5, 5, 21, 23],
                [4, 5, 5, 10, 12],
                [20, 21, 10, 5, 2],
                [22, 23, 12, 2, 5],
            ),
            Float,
        )

    def test_UPGMA_cluster(self):
        """upgma works on pairwise distance dict"""
        pairwise_dist = self.pairwise_distances
        cluster = upgma(pairwise_dist)
        cluster = cluster.sorted()  # so we can make a stable comparison
        expect = make_tree(
            treestring="(((a:0.5,b:0.5)edge.1:1.75,c:2.25)edge.0:5.875,(d:1.0,e:1.0)edge.2:7.125)root;"
        )
        self.assertTrue(cluster.same_topology(expect))

    def test_find_smallest_index(self):
        """find_smallest_index returns the index of smallest value in array"""
        matrix = self.matrix
        index = find_smallest_index(matrix)
        self.assertEqual(index, (0, 1))

    def test_condense_matrix(self):
        """condense_array joins two rows and columns identified by indices"""
        matrix = self.matrix
        index = find_smallest_index(matrix)
        result = condense_matrix(matrix, index, 9999999999)
        assert_allclose(result[0, 0], 5000000.0)
        self.assertEqual(result[1, 4], 9999999999)
        self.assertEqual(result[0, 1], 9999999999)
        self.assertEqual(result[0, 2], 4.5)
        self.assertEqual(result[2, 0], 4.5)
        self.assertEqual(result[0, 4], 22.5)
        self.assertEqual(result[4, 4], 9999999)
        self.assertEqual(result[4, 0], 22.5)

    def test_condense_node_order(self):
        """condense_node_order condenses nodes in list based on index info"""
        matrix = self.matrix
        index = find_smallest_index(matrix)
        node_order = self.node_order
        node_order = condense_node_order(matrix, index, node_order)
        self.assertEqual(node_order[1], None)
        self.assertEqual(node_order[0].__str__(), "(a:0.5,b:0.5);")
        self.assertEqual(node_order[2].__str__(), "c;")
        self.assertEqual(node_order[3].__str__(), "d;")
        self.assertEqual(node_order[4].__str__(), "e;")

    def test_upgma_cluster(self):
        """UPGMA_cluster clusters nodes based on info in a matrix with UPGMA"""
        matrix = self.matrix
        node_order = self.node_order
        large_number = 9999999999
        tree = UPGMA_cluster(matrix, node_order, large_number)
        self.assertEqual(
            str(tree), "(((a:0.5,b:0.5):1.75,c:2.25):5.875,(d:1.0,e:1.0):7.125);"
        )

    def test_UPGMA_cluster_diag(self):
        """UPGMA_cluster works when the diagonal has lowest values"""
        # test that checking the diagonal works
        matrix = self.matrix_zeros
        node_order = self.node_order
        large_number = 9999999999
        tree = UPGMA_cluster(matrix, node_order, large_number)
        self.assertEqual(
            str(tree), "(((a:0.5,b:0.5):1.75,c:2.25):5.875,(d:1.0,e:1.0):7.125);"
        )

    def test_UPGMA_cluster_diag(self):
        """UPGMA_cluster works when the diagonal has intermediate values"""
        # test that checking the diagonal works
        matrix = self.matrix_five
        node_order = self.node_order
        large_number = 9999999999
        tree = UPGMA_cluster(matrix, node_order, large_number)
        self.assertEqual(
            str(tree), "(((a:0.5,b:0.5):1.75,c:2.25):5.875,(d:1.0,e:1.0):7.125);"
        )

    def test_inputs_from_dict_array(self):
        """inputs_from_dict_array makes an array object and PhyloNode list"""
        twod = {
            "1": {"1": 0, "2": 0.86, "3": 0.92},
            "2": {"1": 0.86, "2": 0, "3": 0.67},
            "3": {"1": 0.92, "2": 0.67, "3": 0},
        }
        matrix_d2d = DictArray(twod)

        matrix_array, PhyloNode_order = inputs_from_dict_array(matrix_d2d)
        self.assertEqual(PhyloNode_order[0].name, "1")
        self.assertEqual(PhyloNode_order[2].name, "3")
        assert_allclose(matrix_array[0][2], 0.92)
        assert_allclose(matrix_array[1][0], 0.86)


# run if called from command line
if __name__ == "__main__":
    main()

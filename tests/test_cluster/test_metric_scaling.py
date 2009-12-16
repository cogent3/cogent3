#!usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.cluster.metric_scaling import make_E_matrix, \
        make_F_matrix, run_eig, get_principal_coordinates, \
        principal_coordinates_analysis, output_pca, PCoA
from numpy import array
import numpy
Float = numpy.core.numerictypes.sctype2char(float)

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Catherine Lozupone", "Peter Maxwell", "Rob Knight",
                "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Production"

class MetricScalingTests(TestCase):
    """test the functions to do metric scaling"""
    
    def setUp(self):
        """creates inputs"""
        #create a test input matrix
        self.matrix = array([[1,2,3], [4,5,6]], Float)
        #create a symmetrical matrix
        self.sym_matrix = array([[0.5, 1.0, 1.5, 2.0],
                                [1.0, 2.0, 3.0, 4.0],
                                [1.5, 3.0, 2.0, 1.0],
                                [2.0, 4.0, 1.0, 0.5]])
        
        self.tree_dnd = """((org1:0.11, org2:0.22,(org3:0.12, org4:0.23):0.33)
                :0.2,(org5:0.44, org6:0.55):0.3, org7:0.4);"""
        self.all_envs = ['env1', 'env2', 'env3']
        self.envs_dict = {"org1": array([1,1,0]),
                "org2": array([0,1,0]),
                "org3": array([0,1,0]),
                "org4": array([0,0,1]),
                "org5": array([1,1,0]),
                "org6": array([1,1,0]),
                "org7": array([0,0,1]),
                }
        #sample data set from page 111 of W.J Krzanowski. Principles of
        #multivariate analysis, 2000, Oxford University Press
        self.real_matrix = array([[0,0.099,0.033,0.183,0.148,0.198,0.462,0.628,0.113,0.173,0.434,0.762,0.53,0.586],\
                [0.099,0,0.022,0.114,0.224,0.039,0.266,0.442,0.07,0.119,0.419,0.633,0.389,0.435],\
                [0.033,0.022,0,0.042,0.059,0.053,0.322,0.444,0.046,0.162,0.339,0.781,0.482,0.55], \
                [0.183,0.114,0.042,0,0.068,0.085,0.435,0.406,0.047,0.331,0.505,0.7,0.579,0.53], \
                [0.148,0.224,0.059,0.068,0,0.051,0.268,0.24,0.034,0.177,0.469,0.758,0.597,0.552], \
                [0.198,0.039,0.053,0.085,0.051,0,0.025,0.129,0.002,0.039,0.39,0.625,0.498,0.509], \
                [0.462,0.266,0.322,0.435,0.268,0.025,0,0.014,0.106,0.089,0.315,0.469,0.374,0.369], \
                [0.628,0.442,0.444,0.406,0.24,0.129,0.014,0,0.129,0.237,0.349,0.618,0.562,0.471], \
                [0.113,0.07,0.046,0.047,0.034,0.002,0.106,0.129,0,0.071,0.151,0.44,0.247,0.234], \
                [0.173,0.119,0.162,0.331,0.177,0.039,0.089,0.237,0.071,0,0.43,0.538,0.383,0.346], \
                [0.434,0.419,0.339,0.505,0.469,0.39,0.315,0.349,0.151,0.43,0,0.607,0.387,0.456], \
                [0.762,0.633,0.781,0.7,0.758,0.625,0.469,0.618,0.44,0.538,0.607,0,0.084,0.09], \
               [0.53,0.389,0.482,0.579,0.597,0.498,0.374,0.562,0.247,0.383,0.387,0.084,0,0.038], \
               [0.586,0.435,0.55,0.53,0.552,0.509,0.369,0.471,0.234,0.346,0.456,0.09,0.038,0]])

    #test tree
    dnd = """
    """


    def test_principal_coordinate_analysis(self):
        """principal_coordinate_analysis returns array of principal coors"""
        #I took the example in the book (see intro info), and did the
        #principal coordinates analysis, plotted the data and it looked
        #right
        matrix = self.real_matrix
        pcs, eigvals= principal_coordinates_analysis(matrix)
        bigfirstorder = eigvals.argsort()[::-1]
        pcs = pcs[bigfirstorder]
        eigvals = eigvals[bigfirstorder]
        self.assertEqual(len(pcs), 14)
        self.assertFloatEqual(abs(pcs[0,0]), 0.240788133045)
        self.assertFloatEqual(abs(pcs[1,0]), 0.233677162)

    def test_PCoA(self):
        """PCoA returns a cogent Table result"""
        matrix = self.real_matrix
        names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                    'm', 'n']
        pairwise_dist = {}
        for i, name1 in enumerate(names):
            for j, name2 in enumerate(names):
                combo = (name1, name2)
                if combo not in pairwise_dist:
                    pairwise_dist[combo] = matrix[i,j]

        #perform with the PCoA function
        result = PCoA(pairwise_dist)
        self.assertEqual(result[7,1], 'a')
        self.assertFloatEqual(abs(result[7,2]), 0.240788133045)

    def test_make_E_matrix(self):
        """make_E_matrix converts a distance matrix to an E matrix"""
        matrix = self.matrix
        E_matrix = make_E_matrix(matrix)
        self.assertFloatEqual(E_matrix[0,0], -0.5)
        self.assertFloatEqual(E_matrix[0,1], -2.0)
        self.assertFloatEqual(E_matrix[1,0], -8.0)
        self.assertFloatEqual(E_matrix[1,2], -18.0)

    def test_make_F_matrix(self):
        """make_F_matrix converts an E_matrix to an F_matrix"""
        matrix = self.matrix
        F_matrix = make_F_matrix(matrix)
        self.assertEqual(F_matrix, array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))

    def test_run_eig(self):
        """run_eig returns eigenvectors and values"""
        matrix = self.sym_matrix
        eigvals, eigvecs = run_eig(matrix)
        #make sure that the number of eigvecs and eigvals is equal to dims
        self.assertEqual(len(eigvals), 4)
        self.assertEqual(len(eigvecs), 4)

    def test_get_principal_coordinates(self):
        """get_principal_coordinates normalizes eigvecs with eigvalues"""
        matrix = array([[1,1,1],[2,2,2],[3,3,3]])
        vec = array([0,1,-4])
        result = get_principal_coordinates(vec, matrix)
        self.assertEqual(result[0], array([0,0,0]))
        self.assertEqual(result[1], array([2,2,2]))
        self.assertEqual(result[2], array([6,6,6]))

    def test_output_pca(self):
        """output_pca1 creates a Table output for pcs results"""
        #make arbitary values for inputs
        eigvals = array([4.2,3.2,5.2])
        names = ['env1', 'env2', 'env3']
        pca_matrix = array([[-0.34, -0.22, 0.57], [-0.12, 0.14, -0.018],
                            [1.8, 1.9, 2.0]])
        output_table = output_pca(pca_matrix, eigvals, names)
        self.assertEqual(str(output_table),
'''\
================================================================
        Type              Label  vec_num-0  vec_num-1  vec_num-2
----------------------------------------------------------------
Eigenvectors               env1       1.80      -0.34      -0.12
Eigenvectors               env2       1.90      -0.22       0.14
Eigenvectors               env3       2.00       0.57      -0.02
 Eigenvalues        eigenvalues       5.20       4.20       3.20
 Eigenvalues  var explained (%)      41.27      33.33      25.40
----------------------------------------------------------------''')

    def test_integration(self):
        """Integration test of PCoA should work correctly"""
        u = array([[ 0.        ,  0.7123611 ,  0.76849198,  0.80018563,  0.68248524,  0.74634629],
               [ 0.7123611 ,  0.        ,  0.86645691,  0.80485282,  0.83381301,  0.73881726],
               [ 0.76849198,  0.86645691,  0.        ,  0.82308396,  0.77451746,  0.76498872],
               [ 0.80018563,  0.80485282,  0.82308396,  0.        ,  0.84167365,  0.77614366],
               [ 0.68248524,  0.83381301,  0.77451746,  0.84167365,  0.        ,  0.72661163],
               [ 0.74634629,  0.73881726,  0.76498872,  0.77614366,  0.72661163,  0.        ]])

        e = make_E_matrix(u)
        f = make_F_matrix(e)
        eigvals, eigvecs = run_eig(f)
        bigfirstorder = eigvals.argsort()[::-1]
        #eigvecs = eigvecs[bigfirstorder]
        #eigvals = eigvals[bigfirstorder]
        principal_coords = get_principal_coordinates(eigvals, eigvecs)
        principal_coords = principal_coords[bigfirstorder]
        expected = array([
                   [  2.85970001e-02,  -3.74940557e-01,   3.35175925e-01,
                     -2.54123944e-01,   2.82568441e-01,  -1.72768652e-02],
                   [  2.29038532e-01,   2.23340550e-01,  -2.38559794e-01,
                     -4.12346397e-01,   1.86069108e-01,   1.24580018e-02],
                   [  7.05527166e-02,  -2.08929136e-01,  -3.09988697e-01,
                      2.33436419e-01,   2.88756308e-01,  -7.38276099e-02]])
        self.assertFloatEqual(abs(principal_coords[[0,1,2]]), abs(expected))
                                                                                
#run if called from the command line
if __name__ == '__main__':
    main()


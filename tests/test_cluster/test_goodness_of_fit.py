#!/usr/bin/env python


import unittest

import numpy
import cogent.cluster.goodness_of_fit as goodness_of_fit



__author__ = "Andreas Wilm"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Andreas Wilm"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Andreas Wilm"
__email__ = "andreas.wilm@ucd.ie"
__status__ = "Production"




def example_distmat_and_mdscoords():
    """Return an example distance matrix and corresponding MDS coordinates

    Arguments:
    * None
    Returns:
    * Tuple of:
    ** a distance matrix as numpy.matrix and
    ** MDS coordinates as numpy.array
    """
    distmat = numpy.array([
        [ 0.      ,  0.039806,  0.056853,  0.21595 ,  0.056853,  0.0138  ,
          0.203862,  0.219002,  0.056853,  0.064283],
        [ 0.039806,  0.      ,  0.025505,  0.203862,  0.0208  ,  0.039806,
          0.194917,  0.21291 ,  0.0208  ,  0.027869],
        [ 0.056853,  0.025505,  0.      ,  0.197887,  0.018459,  0.056853,
          0.191958,  0.203862,  0.018459,  0.025505],
        [ 0.21595 ,  0.203862,  0.197887,  0.      ,  0.206866,  0.206866,
          0.07956 ,  0.066935,  0.203862,  0.206866],
        [ 0.056853,  0.0208  ,  0.018459,  0.206866,  0.      ,  0.056853,
          0.203862,  0.21595 ,  0.0138  ,  0.0208  ],
        [ 0.0138  ,  0.039806,  0.056853,  0.206866,  0.056853,  0.      ,
          0.197887,  0.209882,  0.056853,  0.064283],
        [ 0.203862,  0.194917,  0.191958,  0.07956 ,  0.203862,  0.197887,
          0.      ,  0.030311,  0.200869,  0.206866],
        [ 0.219002,  0.21291 ,  0.203862,  0.066935,  0.21595 ,  0.209882,
          0.030311,  0.      ,  0.21291 ,  0.219002],
        [ 0.056853,  0.0208  ,  0.018459,  0.203862,  0.0138  ,  0.056853,
          0.200869,  0.21291 ,  0.      ,  0.011481],
        [ 0.064283,  0.027869,  0.025505,  0.206866,  0.0208  ,  0.064283,
          0.206866,  0.219002,  0.011481,  0.      ]])

    mds_coords = numpy.array([
        [ 0.065233,  0.035019,  0.015413],
        [ 0.059604,  0.00168 , -0.003254],
        [ 0.052371, -0.010959, -0.014047],
        [-0.13804 , -0.036031,  0.031628],
        [ 0.063703, -0.015483, -0.00751 ],
        [ 0.056803,  0.031762,  0.021767],
        [-0.135082,  0.023552, -0.021006],
        [-0.150323,  0.011935, -0.010013],
        [ 0.06072 , -0.01622 , -0.007721],
        [ 0.065009, -0.025254, -0.005257]])

    return (distmat, mds_coords)

    

class GoodnessOfFitTestCase(unittest.TestCase):
        
    def setUp(self):
        """
        set up 
        """
        (self.distmat, self.mds_coords) = example_distmat_and_mdscoords()
        self.stress = goodness_of_fit.Stress(self.distmat, self.mds_coords)

        
    def test_kruskalstress1(self):
        """
        testing goodness_of_fit.calcKruskalStress()
        """
        val = "%0.6f" % self.stress.calcKruskalStress()
        self.assertEqual(val, '0.022555')

            
    def test_sstress(self):
        """
        testing goodness_of_fit.calcSstress()
        """
        val = "%0.6f" % self.stress.calcSstress()
        self.assertEqual(val, '0.008832')


    def test_calc_pwdist(self):
        """
        testing (private) goodness_of_fit._calc_pwdist
        """
        
        # this is a square in 2D
        square_mds = numpy.array([[0,0], [1,0], [1,1], [0,1]])
        # this is what the distance matrix should look like
        square_distmat = numpy.array([[ 0. ,  1. ,  1.41421356,  1. ],
                                      [ 1. ,  0. ,  1. ,  1.41421356],
                                      [ 1.41421356,  1. ,  0. ,  1. ],
                                      [ 1. ,  1.41421356,  1. ,  0. ]])
        derived_distmat = self.stress._calc_pwdist(square_mds)
        # check if dervied and original array are (more or less) the same
        self.assert_((derived_distmat-square_distmat).sum() < 0.000001)


    def test_argument_mixup_exception(self):
        """
        test if mds_coords and distmat are mix-up is detected
        """
        self.assertRaises(AssertionError,
                          goodness_of_fit.Stress,
                          self.mds_coords, self.distmat)
        # should give something like
        # AssertionError: orig_distmat shape bigger than mds_coords shape. Possible argument mixup

        
    def test_size_exception(self):
        """
        test if check on number of rows works
        """
        self.assertRaises(AssertionError,
                           goodness_of_fit.Stress,
                          self.distmat, self.mds_coords.transpose())
        # should give something like
        # AssertionError: orig_distmat and mds_coords do not have the same number of rows/objects.


if __name__ == '__main__':
    unittest.main()


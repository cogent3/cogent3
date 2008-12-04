#!/usr/bin/env python
"""Unit tests for the distance transform functions.

Input file: Frequency data for 3 species (columns) at 3 sites (rows), from Legendre and
Legendre (1998, p. 457):
            10         10          20
            10         15          10
            15          5           5
Output file: Transformation such that Euclidean distances computed among rows of
transformed data are equal to Chord distances among the original sites:
        0.40825    0.40825     0.81650
        0.48507    0.72761     0.48507
        0.90453    0.30151     0.30151
Output file: Transformation such that Euclidean distances computed among rows of
transformed data are equal to Chi-square distances among the original sites:
        0.42258    0.45644      0.84515
        0.48295    0.78246      0.48295
        1.01419    0.36515      0.33806
Output file: Transformation such that Euclidean distances computed among rows of
transformed data are equal to Hellinger distances among the original sites:
        0.50000    0.50000      0.70711
        0.y53452    0.65465      0.53452
        0.77460    0.44721      0.44721
"""
from cogent.maths.distance_transform import (trans_chord, dist_chord,
        trans_chisq, dist_chisq, trans_hellinger, dist_hellinger,
        trans_specprof, dist_specprof)
from cogent.util.unit_test import TestCase, main
from numpy import asmatrix

__author__ = "Zongzhi Liu"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Zongzhi Liu"
__email__ = "zongzhi.liu@gmail.com"
__status__ = "Development"

mat_test = asmatrix([[10, 10, 20],
            [10, 15, 10],
            [15,  5,  5]], float)

class functionTests(TestCase):
    def test_chord_transform(self):
        """trans_chord should return the exp result in the ref paper."""
        exp =  [[ 0.40824829,  0.40824829,  0.81649658],
                [ 0.48507125,  0.72760688,  0.48507125],
                [ 0.90453403,  0.30151134,  0.30151134]]
        res = trans_chord(mat_test)
        self.assertFloatEqual(res, exp)

    def test_chord_dist(self):
        """dist_chord should return the exp result."""
        exp =  [[ 0.        ,  0.46662021,  0.72311971],
                [ 0.46662021,  0.        ,  0.62546036],
                [ 0.72311971,  0.62546036,  0.        ]]
        dist = dist_chord(mat_test)
        self.assertFloatEqual(dist, exp)

    def test_chisq_transform(self):
        """trans_chisq should return the exp result in the ref paper."""
        exp_m = [[ 0.42257713,  0.45643546,  0.84515425],
                [ 0.48294529,  0.7824608 ,  0.48294529],
                [ 1.01418511,  0.36514837,  0.3380617 ]]
        res_m = trans_chisq(mat_test)
        self.assertFloatEqual(res_m, exp_m)

    def test_chisq_distance(self):
        """dist_chisq should return the exp result."""
        exp_d = [[ 0.        ,  0.4910521 ,  0.78452291],
                [ 0.4910521 ,  0.        ,  0.69091002],
                [ 0.78452291,  0.69091002,  0.        ]]
        res_d = dist_chisq(mat_test)
        self.assertFloatEqual(res_d, exp_d)

    def test_hellinger_transform(self):
        """dist_hellinger should return the exp result in the ref paper."""
        exp =  [[ 0.5       ,  0.5       ,  0.70710678],
                [ 0.53452248,  0.65465367,  0.53452248],
                [ 0.77459667,  0.4472136 ,  0.4472136 ]]
        res = trans_hellinger(mat_test)
        self.assertFloatEqual(res, exp)

    def test_hellinger_distance(self):
        """dist_hellinger should return the exp result."""
        exp =  [[ 0.        ,  0.23429661,  0.38175149],
                [ 0.23429661,  0.        ,  0.32907422],
                [ 0.38175149,  0.32907422,  0.        ]]
        dist = dist_hellinger(mat_test)
        self.assertFloatEqual(dist, exp)

    def test_species_profile_transform(self):
        """trans_specprof should return the exp result."""
        exp =  [[ 0.25      ,  0.25      ,  0.5       ],
                [ 0.28571429,  0.42857143,  0.28571429],
                [ 0.6       ,  0.2       ,  0.2       ]]
        res = trans_specprof(mat_test)
        self.assertFloatEqual(res, exp)

    def test_species_profile_distance(self):
        """dist_specprof should return the exp result."""
        exp = [[ 0.        ,  0.28121457,  0.46368092],
               [ 0.28121457,  0.        ,  0.39795395],
               [ 0.46368092,  0.39795395,  0.        ]]
        dist = dist_specprof(mat_test)
        self.assertFloatEqual(dist, exp)


if __name__ == '__main__':
    main()

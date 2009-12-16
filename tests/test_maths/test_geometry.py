#!/usr/bin/env python
"""Tests of the geometry package."""
from numpy import array, take, newaxis
from math import sqrt
from cogent.util.unit_test import TestCase, main
from cogent.maths.geometry import center_of_mass_one_array, \
        center_of_mass_two_array, center_of_mass, distance, sphere_points

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class CenterOfMassTests(TestCase):
    """Tests for the center of mass functions"""

    def setUp(self):
        """setUp for all CenterOfMass tests"""
        self.simple = array([[1, 1, 1], [3, 1, 1], [2, 3, 2]])
        self.simple_list = [[1, 1, 1], [3, 1, 1], [2, 3, 2]]
        self.more_weight = array([[1, 1, 3], [3, 1, 3], [2, 3, 50]])
        self.square = array([[1, 1, 25], [3, 1, 25], [3, 3, 25], [1, 3, 25]])
        self.square_odd = array([[1, 1, 25], [3, 1, 4], [3, 3, 25], [1, 3, 4]])
        self.sec_weight = array([[1, 25, 1], [3, 25, 1], [3, 25, 3], [1, 25, 3]])

    def test_center_of_mass_one_array(self):
        """center_of_mass_one_array should behave correctly"""
        com1 = center_of_mass_one_array
        self.assertEqual(com1(self.simple), array([2, 2]))
        self.assertEqual(com1(self.simple_list), array([2, 2]))
        self.assertFloatEqual(com1(self.more_weight), array([2, 2.785714]))
        self.assertEqual(com1(self.square), array([2, 2]))
        self.assertEqual(com1(self.square_odd), array([2, 2]))
        self.assertEqual(com1(self.sec_weight, 1), array([2, 2]))

    def test_CoM_one_array_wrong(self):
        """center_of_mass_one_array should fail on wrong input"""
        com1 = center_of_mass_one_array
        self.assertRaises(TypeError, com1, self.simple, 'a') #weight_idx wrong
        self.assertRaises(IndexError, com1, self.simple, 100) #w_idx out of range
        self.assertRaises(IndexError, com1, [1, 2, 3], 2) #shape[1] out of range

    def test_center_of_mass_two_array(self):
        """center_of_mass_two_array should behave correctly"""
        com2 = center_of_mass_two_array
        coor = take(self.square_odd, (0, 1), 1)
        weights = take(self.square_odd, (2,), 1)
        self.assertEqual(com2(coor, weights), array([2, 2]))
        weights = weights.ravel()
        self.assertEqual(com2(coor, weights), array([2, 2]))

    def test_CoM_two_array_wrong(self):
        """center_of_mass_two_array should fail on wrong input"""
        com2 = center_of_mass_two_array
        weights = [1, 2]
        self.assertRaises(TypeError, com2, self.simple, 'a') #weight_idx wrong
        self.assertRaises(ValueError, com2, self.simple, weights) #not aligned

    def test_center_of_mass(self):
        """center_of_mass should make right choice between functional methods
        """
        com = center_of_mass
        com1 = center_of_mass_one_array
        com2 = center_of_mass_two_array
        self.assertEqual(com(self.simple), com1(self.simple))
        self.assertFloatEqual(com(self.more_weight), com1(self.more_weight))
        self.assertEqual(com(self.sec_weight, 1), com1(self.sec_weight, 1))
        coor = take(self.square_odd, (0, 1), 1)
        weights = take(self.square_odd, (2,), 1)
        self.assertEqual(com(coor, weights), com2(coor, weights))
        weights = weights.ravel()
        self.assertEqual(com(coor, weights), com2(coor, weights))

    def test_distance(self):
        """distance should return Euclidean distance correctly."""
        #for single dimension, should return difference
        a1 = array([3])
        a2 = array([-1])
        self.assertEqual(distance(a1, a2), 4)
        #for two dimensions, should work e.g. for 3, 4, 5 triangle
        a1 = array([0, 0])
        a2 = array([3, 4])
        self.assertEqual(distance(a1, a2), 5)
        #vector should be the same as itself for any dimensions
        a1 = array([1.3, 23, 5.4, 2.6, -1.2])
        self.assertEqual(distance(a1, a1), 0)
        #should match hand-calculated case for an array
        a1 = array([[1, -2], [3, 4]])
        a2 = array([[1, 0], [-1, 2.5]])
        self.assertEqual(distance(a1, a1), 0)
        self.assertEqual(distance(a2, a2), 0)
        self.assertEqual(distance(a1, a2), distance(a2, a1))
        self.assertFloatEqual(distance(a1, a2), sqrt(22.25))

    def test_sphere_points(self):
        """tests sphere points"""
        self.assertEquals(sphere_points(1), array([[ 1., 0., 0.]]))

#    def test_coords_to_symmetry(self):
#        """tests symmetry expansion (TODO)"""
#        pass
#
#    def test_coords_to_crystal(self):
#        """tests crystal expansion (TODO)"""
#        pass


if __name__ == '__main__':
    main()

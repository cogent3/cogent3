#!/usr/bin/env python
"""Tests of the geometry package."""
from math import sqrt
from unittest import TestCase, main

from numpy import allclose, arange, array, insert, isclose, sum, take
from numpy.linalg import norm
from numpy.random import choice, dirichlet
from numpy.testing import assert_allclose, assert_equal

from cogent3.maths.geometry import (
    aitchison_distance,
    alr,
    alr_inv,
    center_of_mass,
    center_of_mass_one_array,
    center_of_mass_two_array,
    clr,
    clr_inv,
    distance,
    multiplicative_replacement,
    sphere_points,
)


__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight", "Helmut Simon"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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
        assert_equal(com1(self.simple), array([2, 2]))
        assert_equal(com1(self.simple_list), array([2, 2]))
        assert_allclose(com1(self.more_weight), array([2, 2.785714]), rtol=1e-6)
        assert_equal(com1(self.square), array([2, 2]))
        assert_equal(com1(self.square_odd), array([2, 2]))
        assert_equal(com1(self.sec_weight, 1), array([2, 2]))

    def test_CoM_one_array_wrong(self):
        """center_of_mass_one_array should fail on wrong input"""
        com1 = center_of_mass_one_array
        self.assertRaises(TypeError, com1, self.simple, "a")  # weight_idx wrong
        self.assertRaises(IndexError, com1, self.simple, 100)  # w_idx out of range
        # shape[1] out of range
        self.assertRaises(IndexError, com1, [1, 2, 3], 2)

    def test_center_of_mass_two_array(self):
        """center_of_mass_two_array should behave correctly"""
        com2 = center_of_mass_two_array
        coor = take(self.square_odd, (0, 1), 1)
        weights = take(self.square_odd, (2,), 1)
        assert_equal(com2(coor, weights), array([2, 2]))
        weights = weights.ravel()
        assert_equal(com2(coor, weights), array([2, 2]))

    def test_CoM_two_array_wrong(self):
        """center_of_mass_two_array should fail on wrong input"""
        com2 = center_of_mass_two_array
        weights = [1, 2]
        self.assertRaises(TypeError, com2, self.simple, "a")  # weight_idx wrong
        self.assertRaises(ValueError, com2, self.simple, weights)  # not aligned

    def test_center_of_mass(self):
        """center_of_mass should make right choice between functional methods"""
        com = center_of_mass
        com1 = center_of_mass_one_array
        com2 = center_of_mass_two_array
        assert_equal(com(self.simple), com1(self.simple))
        assert_allclose(com(self.more_weight), com1(self.more_weight))
        assert_equal(com(self.sec_weight, 1), com1(self.sec_weight, 1))
        coor = take(self.square_odd, (0, 1), 1)
        weights = take(self.square_odd, (2,), 1)
        assert_equal(com(coor, weights), com2(coor, weights))
        weights = weights.ravel()
        assert_equal(com(coor, weights), com2(coor, weights))

    def test_distance(self):
        """distance should return Euclidean distance correctly."""
        # for single dimension, should return difference
        a1 = array([3])
        a2 = array([-1])
        self.assertEqual(distance(a1, a2), 4)
        # for two dimensions, should work e.g. for 3, 4, 5 triangle
        a1 = array([0, 0])
        a2 = array([3, 4])
        self.assertEqual(distance(a1, a2), 5)
        # vector should be the same as itself for any dimensions
        a1 = array([1.3, 23, 5.4, 2.6, -1.2])
        self.assertEqual(distance(a1, a1), 0)
        # should match hand-calculated case for an array
        a1 = array([[1, -2], [3, 4]])
        a2 = array([[1, 0], [-1, 2.5]])
        self.assertEqual(distance(a1, a1), 0)
        self.assertEqual(distance(a2, a2), 0)
        self.assertEqual(distance(a1, a2), distance(a2, a1))
        assert_allclose(distance(a1, a2), sqrt(22.25))

    def test_sphere_points(self):
        """tests sphere points"""
        assert_equal(sphere_points(1), array([[1.0, 0.0, 0.0]]))


class TestAitchison(TestCase):
    def setUp(self):
        x = choice(20, size=10) + 0.1
        self.x = x
        a = arange(1, 7)
        self.a = a
        d = dirichlet(a, size=2)
        self.d = d

    def test_Aitchison_transforms(self):
        """Test that alr_inv of alr is in fact the inverse
        of alr. Ditto for clr_inv and clr. Then test that clr
        transforms into hyperplane x1 + ... + xn=0."""
        length = len(self.x)
        for col in range(-1, length):
            y = alr_inv(self.x, col)
            assert allclose(self.x, alr(y, col)), (
                "Failed alr inverse test for col = " + str(col) + "."
            )

        z = dirichlet(self.x)
        y = clr(z)
        assert allclose(z, clr_inv(y)), "Failed clr inverse test."
        assert allclose(sum(y), 0), "Failed clr hyperplane test."

    def test_Aitchison_distance(self):
        x = self.d[0]
        y = self.d[1]
        assert allclose(
            aitchison_distance(x, y), norm(clr(x) - clr(y))
        ), "Failed distance test."

    def test_multiplicative_replacement(self):
        x1 = dirichlet(self.a)
        y1 = insert(x1, 3, 0)
        u = multiplicative_replacement(y1)
        assert allclose(
            y1, u, atol=1e-2
        ), "Multiplicative replacement peturbation is too large."
        assert isclose(
            sum(u), 1
        ), "Multiplicative replacement does not yield a composition."


if __name__ == "__main__":
    main()

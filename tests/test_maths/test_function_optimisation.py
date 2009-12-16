#!/usr/bin/env python
"""Tests for optimisation functions"""
from function_optimisation import great_deluge, ga_evolve, _simple_breed,\
    _simple_score, _simple_init, _simple_select
from cogent.util.unit_test import TestCase, main
from numpy import array

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Production"

class OptimisationFunctionsTestCase(TestCase):
    """Tests for great_delue and ga_evolve"""
    def test_great_deluge(self):
        """great_deluge should return expected values from foo() obj"""
        class foo:
            def __init__(self, x): self.x = x
            def cost(self): return self.x
            def perturb(self): return self.__class__(self.x - 1)
        
        observed = [i for i in great_deluge(foo(5), max_total_iters=6)]

        self.assertEqual(observed[0][1].x, 4)
        self.assertEqual(observed[1][1].x, 3)
        self.assertEqual(observed[2][1].x, 2)
        self.assertEqual(observed[3][1].x, 1)
        self.assertEqual(observed[4][1].x, 0)
        self.assertEqual(observed[5][1].x, -1)

    def test_ga_evolve(self):
        """ga_evolve should return expected values when using overloaded funcs"""
        init_f = lambda x,y: [1,1,1]
        score_f = lambda x,y: 5
        breed_f = lambda w,x,y,z: [1,1,1]
        select_f = lambda x,y: 2

        expected = [(0, 2), (1, 2), (2, 2)]
        observed = [i for i in ga_evolve(1, 2, 3, 0.5, score_f, breed_f, \
                    select_f, init_f, None, 3)]
        self.assertEqual(observed, expected)        

class PrivateFunctionsTestCase(TestCase):
    """Tests of the private support functions for ga_evolve"""
    def test_simple_breed(self):
        """simple_breed should return expected values when with modded parent"""
        f = lambda: 0.5
        obj = lambda: 0
        obj.mutate = lambda: 1
        expected = [1,1,1,1,1]
        observed = _simple_breed([0, obj], 5, 1.0, f)
        self.assertEqual(observed, expected)
       
    def test_simple_score(self):
        """simple_score should return choosen value with overloaded obj"""
        bar = lambda: 5
        bar.score = lambda x: x
        self.assertEqual(_simple_score(bar,6), 6)

    def test_simple_init(self):
        """simple_init should return a simple list"""
        expected = [array([0]), array([0]), array([0])]
        self.assertEqual(_simple_init(array([0]), 3), expected)

    def test_simple_select(self):
        """simple_select should return our hand picked selection"""
        pop = ['a','b','c','d','e']
        scores = [5,3,8,6,1]
        best_expected = (1, 'e')
        self.assertEqual(_simple_select(pop, scores), best_expected)

if __name__ == '__main__':
    main()

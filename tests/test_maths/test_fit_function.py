#!/usr/bin/env python
"""Unit tests for fit function.
"""
from numpy import array, arange, exp
from numpy.random import rand
from cogent.util.unit_test import TestCase, main
from cogent.maths.fit_function import fit_function

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.1"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Prototype"

class fit_function_test(TestCase):
    """Tests of top-level fit functions."""

    def test_constant(self):
        """test constant approximation"""
        # defining our fitting function
        def f(x,a):
            return a[0]
        
        exp_params = [2]
        x = arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(x))
        
        params = fit_function(x, y_noise, f, 1, 5)
        
        self.assertFloatEqual(params, exp_params , .5)
        
    def test_linear(self):
        """test linear approximation"""
        # defining our fitting function
        def f(x,a):
            return (a[0]+x*a[1])
        
        exp_params = [2, 10]
        x = arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(y))
        
        params = fit_function(x, y_noise, f, 2, 5)
        
        self.assertFloatEqual(params, exp_params , .5)
        
    def test_exponential(self):
        """test exponential approximation"""
        # defining our fitting function
        def f(x,a):
            return exp(a[0]+x*a[1])
        
        exp_params = [2, 10]
        x = arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(y))
        
        params = fit_function(x, y_noise, f, 2, 5)
        
        self.assertFloatEqual(params, exp_params , .5)

if __name__ == '__main__':
    main()

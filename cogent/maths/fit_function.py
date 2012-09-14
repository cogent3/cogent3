#!/usr/bin/env python
""" fitting funtions
     
module to fit x and y samples to a model

"""
from __future__ import division
from numpy import array
from cogent.maths.scipy_optimize import fmin

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Prototype"

def fit_function(x_vals, y_vals, func, n_params, iterations=2):
   """ Fit any function to any array of values of x and y.
   :Parameters:
       x_vals : array
           Values for x to fit the function func.
       y_vals : array
           Values for y to fit the function func.
       func : callable ``f(x, a)``
           Objective function (model) to be fitted to the data. This function 
           should return either an array for models that are not a constant, 
           i.e. f(x)=exp(a[0]+x*a[1]), or a single value for models that are a
           cosntant, i.e. f(x)=a[0]
       n_params : int
           Number of parameters to fit in func
       iterations : int
           Number of iterations to fit func

   :Returns: param_guess

       param_guess : array
           Values for each of the arguments to fit func to x_vals and y_vals

   :Notes:

       Fit a function to a given array of values x and y using simplex to
       minimize the error.

   """

   # internal function to minimize the error
   def f2min(a):
       #sum square deviation
       return ((func(x_vals, a) - y_vals)**2).sum()

   param_guess = array(range(n_params))
   for i in range(iterations):
       xopt = fmin(f2min, param_guess, disp=0)
       param_guess = xopt

   return xopt


#if __name__ == "__main__":
#    main()

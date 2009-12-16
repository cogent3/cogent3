#!/usr/bin/env python
"""
optimiser.py

Contains a base class for numerical optimisers.
"""

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
import random
import logging

from cogent.util.checkpointing import Checkpointer

__author__ = "Andrew Butterfield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Edward Lang",
                    "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

LOG = logging.getLogger('cogent')

class ParameterOutOfBoundsError(Exception):
    pass

epsilon = 1e-9 # intended for dealing with rounding errors

# The following functions are used to wrap the optimised function to
# adapt it to the optimiser in various ways.  They can be combined.

def bounded_function(f, lower_bounds, upper_bounds):
    """Returns a function that return +inf on out-of-bounds input rather than
    bothering the real function with invalid input.  This is enough to get some
    unbounded optimisers working on bounded problems"""
    
    def _wrapper(x, **kw):
        if numpy.alltrue(numpy.logical_and(lower_bounds <= x, x <= upper_bounds)):
            return f(x, **kw)
        else:
            raise ParameterOutOfBoundsError((lower_bounds, x, upper_bounds))
    
    return _wrapper

def bounds_exception_catching_function(f, direction):
    """Like bounded_function, except that the out-of-bounds parameter isn't caught
    here, but deeper in the calculation, and reported as an exception"""
    if direction < 0:
        out_of_bounds_value = numpy.inf
        acceptable_inf = numpy.isposinf
    else:
        out_of_bounds_value = -numpy.inf
        acceptable_inf = numpy.isneginf
    
    def _wrapper(x, **kw):
        try:
            result = f(x, **kw)
            if not numpy.isfinite(result):
                if not acceptable_inf(result):
                    LOG.warning('Non-finite f %s from %s' % (result, x))
                    raise ParameterOutOfBoundsError
        except (ArithmeticError, ParameterOutOfBoundsError), detail:
            result = out_of_bounds_value
        return result
    
    return _wrapper


def reversed_function(f):
    """Likelihood needs maximising, but our optimisers are minimisers"""
    def _wrapper(x, **kw):
        return -f(x, **kw)
    
    return _wrapper


class OptimiserBase(object):
    """Maximises (or minimises if 'direction'<0) optimisable_object
    .testoptparvector(x) for x within optimisable_object.getBoundsVectors().
    For some optimisers random numbers are used to choose the search path,
    so for deterministic behaviour provide either 'random_series' or 'seed'.
    """
    # Subclass must provide _setdefaults, setConditions, runInner
    # runInner will get passed a function to minimise, so it does
    # not need to call any OptimiserBase methods.
    
    def __init__(self, f, xinit, bounds, direction=1, random_series=None,
            seed=None, **kw):
        self.restore = True
        self._random_series = random_series or random.Random()
        if seed is not None:
            self._random_series.seed(seed)
        self.__total_evaluations = 0
        # flag to determine the direction. We assume that all wrapped
        # optimiser functions are minimisers.
        self.__reverse_result = direction * self.algorithm_direction < 0
        self.setCheckpointing()
        self._setdefaults()
        self.setConditions(**kw)
        self.setVectorBounds(*bounds)
        self.vector_length = len(xinit)   #risks irrelevance
        self.vector = numpy.array(xinit, Float)
        self.__original_f = f
    
    def setCheckpointing(self, filename = None, interval = 1800, restore = None):
        """
        Set the checkpointing filename and time interval.
        Arguments:
        - filename: name of the file to which data will be written. If None, no
          checkpointing will be done.
        - interval: time expressed in seconds
        - restore: flag to restore from this filenamne or not. will be set to 0 after
          restoration
        """
        self.checkpointer = Checkpointer(filename, interval)
        if restore is not None:
            self.restore = restore
    
    def setVectorBounds(self, lower_bounds = None, upper_bounds = None):
        """Set the bounds within which the vector values can lie."""
        if lower_bounds is not None:
            self.lower_bounds = lower_bounds
        if upper_bounds is not None:
            self.upper_bounds = upper_bounds
    
    def getEvaluationCount(self):
        """Return the total number of evaluations taken."""
        
        return self.__total_evaluations
    
    def run(self, show_progress = True):
        """
        In principle this would call the virtually overriden function, ie be the
        public interface to the protected overridable.
        
        Arguments:
        - show_progress: whether the function values are printed as
          the optimisation proceeds. Default is True.
        
        Returns the optimised function value and the corresponding parameter
        vector.
        """
        
        f = self.__original_f
        vector = self.vector
        
        # adapt the function to the optimiser
        if self.__reverse_result:
            f = reversed_function(f)
        f = bounded_function(f, self.lower_bounds, self.upper_bounds)
        try:
            fval = f(vector)
        except (ArithmeticError, ParameterOutOfBoundsError), detail:
            raise ValueError("Initial parameter values must be valid %s" % repr(detail.args))
        if not numpy.isfinite(fval):
            if numpy.isinf(fval) and self.__reverse_result:
                fval = -1 * fval
            raise ValueError("Initial parameter values must evaluate to a finite likelihood, not %s. %s" % (fval, vector))
        
        f = bounds_exception_catching_function(f, self.algorithm_direction)
        
        (fval, xopt, func_calls, elapsed_time) = self.runInner(f, vector,
                show_progress=show_progress, random_series = self._random_series)
        
        # ensure state of calculator reflects optimised result
        f(xopt)
        
        self.__total_evaluations += func_calls
        self.elapsed_time = elapsed_time
        
        if self.__reverse_result:
            fval = -fval
        
        return fval, xopt
    
    def setVector(self, vector):
        self.vector = vector
        self.__original_f(vector)
    

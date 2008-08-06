#!/usr/bin/env python
from cogent.maths.scipy_optimize import fmin_bfgs, fmin_powell, fmin, brent,\
    bracket, golden, fpconst
from optimiser import OptimiserBase
import time

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def bound_brent(func, args, brack=None, tol=1.48e-8, full_output=0, **kw):
    """Given a function and an initial point, find another
    point within the bounds, then use the two points to
    bracket a minimum.
    
    This differs from ordinary brent() only in that it protects
    bracket() from infinities.  bracket() may find an infinity as the
    third point, but that's OK because it finishes as soon as that happens.
    
    If bracket() returns an invalid 3rd point then we will pass it on
    to brent(), but brent() knows to use golden section until all 3
    points are finite so it will cope OK.
    """
    
    assert not brack, brack
    xa = 0.0
    fa = apply(func, (xa,)+args)
    assert fa is not fpconst.PosInf, "Starting point is infinite"
    
    # if dx sends us over the boundry shrink and reflect it until
    # it doesn't any more.
    dx = -2.0  # this would be -2.0 in orig impl, but is smaller better?
    xb = xa + dx
    fb = fpconst.PosInf
    while fb is fpconst.PosInf and xb != xa:
        dx = dx * -0.5
        xb = xa + dx
        fb = apply(func, (xb,)+args)
    assert xb != xa, "Can't find a second in-bounds point on this line"
    return brent(func, args, (xa, xb), tol, full_output, **kw)


class _SciPyOptimiser(OptimiserBase):
    """This class is abstract.  Subclasses must provide a
    _minimise(self, f, x) that can sanely handle +inf.
    
    Since these are local optimisers, we sometimes restart them to
    check the result is stable.  Cost is less than 2-fold slowdown"""
    
    default_max_restarts = 0
    # These are minimisers
    algorithm_direction = -1
    
    def _setdefaults(self):
        self.max_restarts = self.default_max_restarts
        self.ftol = 1e-6
    
    def setConditions(self, max_restarts=None, tolerance=None,
                      max_evaluations=None):
        if max_restarts is not None:
            self.max_restarts = max_restarts
        if tolerance is not None:
            self.ftol = tolerance
        self.max_evaluations = max_evaluations
    
    def runInner(self, function, xopt, show_progress, random_series=None):
        # random_series isn't used - these optimisers are deterministic -
        # but optimiser_base doesn't know that.
        fval_last = fval = fpconst.PosInf
        total_evaluations = 0
        t0 = time.time()
        if len(xopt) == 0:
            return function(xopt), xopt, 1, time.time() - t0
        for i in range((self.max_restarts + 1)):
            (xopt, fval, iterations, func_calls, warnflag) = self._minimise(
                    function, xopt, disp=show_progress,
                    ftol=self.ftol, full_output=True,
                    maxfun=self.max_evaluations)
            if warnflag:
                print "FORCED EXIT from optimiser after %s evaluations" % \
                        self.max_evaluations
            total_evaluations += func_calls
            
            # same check as in fmin_powell
            if abs(fval_last - fval) < self.ftol:
                break
            fval_last = fval  # fval <= fval_last
        
        # Correct the sign of the result. If we reversed the direction of
        # the function for the benefit of the optimiser then we now have to
        # flip it back.
        # fval = self._optimisable_object.direction * fval
        
        return fval, xopt, total_evaluations, time.time() - t0
    

class Powell(_SciPyOptimiser):
    """Uses an infinity avoiding version of the Brent line search."""
    def _minimise(self, f, x, **kw):
        result = fmin_powell(f, x, linesearch=bound_brent, **kw)
        # same length full-results tuple as simplex:
        (xopt, fval, directions, iterations, func_calls, warnflag) = result
        return (xopt, fval, iterations, func_calls, warnflag)
    

class BoundPowell(_SciPyOptimiser):
    """Uses a line search between the bounds.  Only works for fully bounded problems!
    And seems slower than the other Powell class, which only avoids bounds when it
    hits them."""
    def _minimise(self, f, x, **kw):
        (lower_bounds, upper_bounds) = self._optimisable_object.getbounds()
        result = fmin_powell(f, x,
                linesearch=None,
                bounds=(lower_bounds, upper_bounds),
                **kw)
        # same length full-results tuple as simplex:
        (xopt, fval, directions, iterations, func_calls, warnflag) = result
        return (xopt, fval, iterations, func_calls, warnflag)
    

class DownhillSimplex(_SciPyOptimiser):
    """On a small brca1 tree this fails to find a minimum as good as the
    other optimisers.  Restarts help a lot though."""
    default_max_restarts = 5
    def _minimise(self, f, x, **kw):
        return fmin(f, x, **kw)
    

DefaultLocalOptimiser = Powell

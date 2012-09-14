#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Local or Global-then-local optimisation with progress display
"""

from cogent.util import progress_display as UI
from simannealingoptimiser import SimulatedAnnealing
from scipy_optimisers import DownhillSimplex, Powell
import warnings
import numpy

GlobalOptimiser = SimulatedAnnealing
LocalOptimiser = Powell

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def unsteadyProgressIndicator(display_progress, label='', start=0.0, end=1.0):
    template = u'f = % #10.6g  ±  % 9.3e   evals = %6i '
    label = label.rjust(5)
    goal = [1.0e-20]
    def _display_progress(remaining, *args):
        if remaining > goal[0]:
            goal[0] = remaining
        progress = (goal[0]-remaining)/goal[0] * (end-start) + start
        msg = template % args + label
        return display_progress(msg, progress=progress, current=0)
    return _display_progress


class ParameterOutOfBoundsError(Exception):
    pass

class MaximumEvaluationsReached(Exception):
    pass
    

# The following functions are used to wrap the optimised function to
# adapt it to the optimiser in various ways.  They can be combined.

def limited_use(f, max_evaluations=None):
    if max_evaluations is None:
        max_evaluations = numpy.inf
    evals = [0]
    best_fval = [-numpy.inf]
    best_x = [None]
    def wrapped_f(x):
        if evals[0] >= max_evaluations:
            raise MaximumEvaluationsReached(evals[0])
        evals[0] += 1
        fval = f(x)
        if fval > best_fval[0]:
            best_fval[0] = fval
            best_x[0] = x.copy()
        return fval
    def get_best():
        f(best_x[0])  # for calculator, ensure best last
        return best_fval[0], best_x[0], evals[0]
    return get_best, wrapped_f

def bounded_function(f, lower_bounds, upper_bounds):
    """Returns a function that raises an exception on out-of-bounds input 
    rather than bothering the real function with invalid input.  
    This is enough to get some unbounded optimisers working on bounded problems"""
    
    def _wrapper(x, **kw):
        if numpy.alltrue(numpy.logical_and(lower_bounds <= x, x <= upper_bounds)):
            return f(x, **kw)
        else:
            raise ParameterOutOfBoundsError((lower_bounds, x, upper_bounds))
    
    return _wrapper

def bounds_exception_catching_function(f):
    """Returns a function that return -inf on out-of-bounds or otherwise 
    impossible to evaluate input.  This only helps if the function is to be 
    MAXIMISED."""
    out_of_bounds_value = -numpy.inf
    acceptable_inf = numpy.isneginf
    
    def _wrapper(x, **kw):
        try:
            result = f(x, **kw)
            if not numpy.isfinite(result):
                if not acceptable_inf(result):
                    warnings.warn('Non-finite f %s from %s' % (result, x))
                    raise ParameterOutOfBoundsError
        except (ArithmeticError, ParameterOutOfBoundsError), detail:
            result = out_of_bounds_value
        return result
    
    return _wrapper

def minimise(f, *args, **kw):
    """See maximise"""
    def nf(x):
        return -1 * f(x)
    return maximise(nf, *args, **kw)
    
@UI.display_wrap
def maximise(f, xinit, bounds=None, local=None, filename=None, interval=None,
        max_restarts=None, max_evaluations=None, limit_action='warn',
        tolerance=1e-6, global_tolerance=1e-1, ui=None,
        return_eval_count=False,
        **kw):
    """Find input values that optimise this function.
    'local' controls the choice of optimiser, the default being to run
    both the global and local optimisers. 'filename' and 'interval'
    control checkpointing.  Unknown keyword arguments get passed on to
    the global optimiser.
    """
    do_global = (not local) or local is None
    do_local = local or local is None
    
    assert limit_action in ['ignore', 'warn', 'raise', 'error']
    (get_best, f) = limited_use(f, max_evaluations)

    x = numpy.array(xinit, float)
    multidimensional_input = x.shape != ()
    if not multidimensional_input:
        x = numpy.atleast_1d(x)
    
    if bounds is not None:
        (upper, lower) = bounds
        if upper is not None or lower is not None:
            if upper is None: upper = numpy.inf
            if lower is None: lower = -numpy.inf
            f = bounded_function(f, upper, lower)
    try:
        fval = f(x)
    except (ArithmeticError, ParameterOutOfBoundsError), detail:
        raise ValueError("Initial parameter values must be valid %s" % repr(detail.args))
    if not numpy.isfinite(fval):
        raise ValueError("Initial parameter values must evaluate to a finite value, not %s. %s" % (fval, x))
    
    f = bounds_exception_catching_function(f)
    
    try:
        # Global optimisation 
        if do_global:
            if 0 and not do_local:
                warnings.warn(
                    'local=False causes the post-global optimisation local '
                    '"polishing" optimisation to be skipped entirely, which seems '
                    'pointless, so its meaning may change to a simple boolean '
                    'flag: local or global.')
                # It also needlessly complicates this function.
                gend = 1.0
            else:
                gend = 0.9
            callback = unsteadyProgressIndicator(ui.display, 'Global', 0.0, gend)
            gtol = [tolerance, global_tolerance][do_local]
            opt = GlobalOptimiser(filename=filename, interval=interval)
            x = opt.maximise(f, x, tolerance=gtol, 
                    show_remaining=callback, **kw)
        else:
            gend = 0.0
            for k in kw:
                warnings.warn('Unused arg for local alignment: ' + k)

        # Local optimisation
        if do_local:
            callback = unsteadyProgressIndicator(ui.display, 'Local', gend, 1.0)
            #ui.display('local opt', 1.0-per_opt, per_opt)
            opt = LocalOptimiser()
            x = opt.maximise(f, x, tolerance=tolerance, 
                    max_restarts=max_restarts, show_remaining=callback)
    finally:
        # ensure state of calculator reflects optimised result, or
        # partialy optimised result if exiting on an exception.
        (f, x, evals) = get_best()
                    
    # ... and returning this info the obvious way keeps this function 
    # potentially applicable optimising non-caching pure functions too.
    if not multidimensional_input:
        x = numpy.squeeze(x)
    
    if return_eval_count:
        return x, evals
    
    return x


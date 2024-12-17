"""Local or Global-then-local optimisation with progress display"""

from __future__ import annotations

import warnings
from functools import singledispatch

import numpy

from cogent3.util import progress_display as UI

from .scipy_optimisers import Powell
from .simannealingoptimiser import SimulatedAnnealing

GlobalOptimiser = SimulatedAnnealing
LocalOptimiser = Powell


@singledispatch
def _standardise_data(val) -> tuple[float]:
    return (str(val),)


@_standardise_data.register
def _(val: numpy.ndarray) -> tuple[float]:
    val = val.tolist() if val.ndim else [val.tolist()]
    return tuple(val)


@_standardise_data.register
def _(val: float) -> tuple[float]:
    return (val,)


def unsteadyProgressIndicator(display_progress, label="", start=0.0, end=1.0):
    template = "f = % #10.6g  ±  % 9.3e   evals = %6i "
    label = label.rjust(5)
    goal = [1.0e-20]

    def _display_progress(remaining, *args):
        # the first element is a numpy array, which can be
        # a scalar array, i.e. ndim==0. We use the tolist() method to
        # cast to a standard python type.
        v1 = _standardise_data(args[0])
        args = v1 + args[1:]
        goal[0] = max(remaining, goal[0])
        progress = (goal[0] - remaining) / goal[0] * (end - start) + start
        msg = template % args + label
        return display_progress(msg, progress=progress)

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


def bounded_function(f, lower_bounds, upper_bounds, report_error=False):
    """Returns a function that raises an exception on out-of-bounds input
    rather than bothering the real function with invalid input.
    This is enough to get some unbounded optimisers working on bounded problems"""

    def _wrapper(x, **kw):
        if numpy.all(numpy.logical_and(lower_bounds <= x, x <= upper_bounds)):
            return f(x, **kw)
        pos = numpy.logical_or(x <= lower_bounds, x >= upper_bounds)
        lower = numpy.array(lower_bounds)
        upper = numpy.array(upper_bounds)
        raise ParameterOutOfBoundsError((lower[pos], x[pos], upper[pos]))

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
            if not numpy.isfinite(result) and not acceptable_inf(result):
                warnings.warn(f"Non-finite f {result} from {x}", stacklevel=2)
                raise ParameterOutOfBoundsError
        except (ArithmeticError, ParameterOutOfBoundsError):
            result = out_of_bounds_value
        return result

    return _wrapper


def minimise(f, *args, **kw):
    """See maximise"""

    def nf(x):
        return -1 * f(x)

    return maximise(nf, *args, **kw)


@UI.display_wrap
def maximise(
    f,
    xinit,
    bounds=None,
    local=None,
    filename=None,
    interval=None,
    max_restarts=None,
    max_evaluations=None,
    tolerance=1e-6,
    global_tolerance=1e-1,
    ui=None,
    return_eval_count=False,
    warn=False,
    **kw,
):
    """
    Find input values that optimise this function.
    'local' controls the choice of optimiser, the default being to run
    both the global and local optimisers. 'filename' and 'interval'
    control checkpointing.  Unknown keyword arguments get passed on to
    the global optimiser.
    Parameters
    ----------
    f
        callable
    xinit
        initial parameter values
    bounds
        (upper, lower) bounds with a value in each for each element of xinit
    local
        Controls which optimiser(s) to apply. If True, uses Powell local
        optimiser only. If False, uses only the global (Simulated Annealing)
        optimiser. If None, uses both global and local.
    filename
        file to checkpoint optimisatiuon results to.
    interval
        time expressed in seconds
    max_restarts
        maximum number of times to try optimisation
    max_evaluations
        maximum number of function evaluations to perform
    tolerance
        exit condition for local optimiser
    global_tolerance
        exit condition for global optimiser
    return_eval_count
        number of times function evaluated
    warn
        whether to report unused variables
    kw
        kwargs passed through to global optimiser only

    Returns
    -------
    The vector of values maximising f. Optionally returns the number of
    function evaluations.
    """
    do_global = (not local) or local is None
    do_local = local or local is None

    (get_best, f) = limited_use(f, max_evaluations)

    x = numpy.array(xinit, float)
    multidimensional_input = x.shape != ()
    if not multidimensional_input:
        x = numpy.atleast_1d(x)

    if bounds is not None:
        upper, lower = bounds
        if upper is not None or lower is not None:
            if upper is None:
                upper = numpy.inf
            if lower is None:
                lower = -numpy.inf
            f = bounded_function(f, upper, lower)
    try:
        fval = f(x)
    except (ArithmeticError, ParameterOutOfBoundsError) as detail:
        msg = f"Initial parameter values must be valid {detail.args!r}"
        raise ValueError(
            msg,
        ) from detail

    if not numpy.isfinite(fval):
        msg = (
            f"Initial parameter values must evaluate to a finite value, not {fval}. {x}"
        )
        raise ValueError(
            msg,
        )

    f = bounds_exception_catching_function(f)

    try:
        # Global optimisation
        if do_global:
            gend = 0.9
            callback = unsteadyProgressIndicator(ui.display, "Global", 0.0, gend)
            gtol = [tolerance, global_tolerance][do_local]
            opt = GlobalOptimiser(filename=filename, interval=interval)
            x = opt.maximise(f, x, tolerance=gtol, show_remaining=callback, **kw)
        else:
            gend = 0.0
            if warn:
                warnings.warn(f"Unused args for local optimisation: {kw}", stacklevel=2)

        # Local optimisation
        if do_local:
            callback = unsteadyProgressIndicator(ui.display, "Local", gend, 1.0)
            # ui.display('local opt', 1.0-per_opt, per_opt)
            opt = LocalOptimiser()
            x = opt.maximise(
                f,
                x,
                tolerance=tolerance,
                max_restarts=max_restarts,
                show_remaining=callback,
            )
    finally:
        # ensure state of calculator reflects optimised result, or
        # partially optimised result if exiting on an exception.
        (f, x, evals) = get_best()

    # ... and returning this info the obvious way keeps this function
    # potentially applicable optimising non-caching pure functions too.
    if not multidimensional_input:
        x = numpy.squeeze(x)

    if return_eval_count:
        return x, evals

    return x

#!/usr/bin/env python
import math

import numpy

from cogent3.maths.scipy_optimize import brent, fmin_powell


def bound_brent(func, brack=None, **kw):
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
    fa = func(xa)
    assert fa is not numpy.inf, "Starting point is infinite"

    # if dx sends us over the boundry shrink and reflect it until
    # it doesn't any more.
    dx = -2.0  # this would be -2.0 in orig impl, but is smaller better?
    xb = xa + dx
    fb = numpy.inf
    while fb is numpy.inf and xb != xa:
        dx = dx * -0.5
        xb = xa + dx
        fb = func(xb)
    assert xb != xa, "Can't find a second in-bounds point on this line"
    return brent(func, brack=(xa, xb), **kw)


class _SciPyOptimiser:
    """This class is abstract.  Subclasses must provide a
    _minimise(self, f, x) that can sanely handle +inf.

    Since these are local optimisers, we sometimes restart them to
    check the result is stable.  Cost is less than 2-fold slowdown"""

    def maximise(self, function, *args, **kw):
        def nf(x):
            return -1 * function(x)

        return self.minimise(nf, *args, **kw)

    def minimise(
        self,
        function,
        xopt,
        show_remaining,
        max_restarts=None,
        tolerance=None,
    ):
        if max_restarts is None:
            max_restarts = 0
        if tolerance is None:
            tolerance = 1e-6

        fval_last = numpy.inf
        if len(xopt) == 0:
            return function(xopt), xopt

        if show_remaining:

            def _callback(fcalls, x, fval, delta) -> None:
                remaining = math.log(max(abs(delta) / tolerance, 1.0))
                show_remaining(remaining, -fval, delta, fcalls)

        else:
            _callback = None

        for _i in range(max_restarts + 1):
            (xopt, fval, iterations, func_calls, warnflag) = self._minimise(
                function,
                xopt,
                disp=False,
                callback=_callback,
                ftol=tolerance,
                full_output=True,
            )

            xopt = numpy.atleast_1d(xopt)  # unsqueeze incase only one param

            # same tolerance check as in fmin_powell
            if abs(fval_last - fval) < tolerance:
                break
            fval_last = fval  # fval <= fval_last

        return xopt


class Powell(_SciPyOptimiser):
    """Uses an infinity avoiding version of the Brent line search."""

    def _minimise(self, f, x, **kw):
        result = fmin_powell(f, x, linesearch=bound_brent, **kw)
        # same length full-results tuple as simplex:
        (xopt, fval, directions, iterations, func_calls, warnflag) = result
        return (xopt, fval, iterations, func_calls, warnflag)


DefaultLocalOptimiser = Powell

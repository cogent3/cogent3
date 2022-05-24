# We don't want to depend on the monolithic, fortranish,
# Num-overlapping, mac-unfriendly SciPy.  But this
# module is too good to pass up. It has been lightly customised for
# use in Cogent.  Changes made to fmin_powell and brent: allow custom
# line search function (to allow bound_brent to be passed in), cope with
# infinity, tol specified as an absolute value, not a proportion of f,
# and more info passed out via callback.

# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************


# Minimization routines

__all__ = ["fmin_powell", "brent", "bracket"]

import builtins

import numpy

from numpy import absolute, asarray, atleast_1d, eye, sqrt, squeeze


# These have been copied from Numeric's MLab.py
# I don't think they made the transition to scipy_core


def max(m, axis=0):
    """max(m,axis=0) returns the maximum of m along dimension axis."""
    m = asarray(m)
    return numpy.maximum.reduce(m, axis)


def min(m, axis=0):
    """min(m,axis=0) returns the minimum of m along dimension axis."""
    m = asarray(m)
    return numpy.minimum.reduce(m, axis)


def is_array_scalar(x):
    """Test whether `x` is either a scalar or an array scalar."""
    return len(atleast_1d(x) == 1)


abs = absolute

pymin = builtins.min
pymax = builtins.max
__version__ = "2022.5.25a1"


_epsilon = sqrt(numpy.finfo(float).eps)

__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def wrap_function(function, args):
    ncalls = [0]

    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)

    return ncalls, function_wrapper


class Brent:
    # need to rethink design of __init__

    def __init__(self, func, tol=1.48e-8, maxiter=500):
        self.func = func
        self.tol = tol
        self.maxiter = maxiter
        self._mintol = 1.0e-11
        self._cg = 0.3819660
        self.xmin = None
        self.fval = None
        self.iter = 0
        self.funcalls = 0
        self.brack = None
        self._brack_info = None

    # need to rethink design of set_bracket (new options, etc)
    def set_bracket(self, brack=None):
        self.brack = brack
        self._brack_info = self.get_bracket_info()

    def get_bracket_info(self):
        # set up
        func = self.func
        brack = self.brack
        ### BEGIN core bracket_info code ###
        ### carefully DOCUMENT any CHANGES in core ##
        if brack is None:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func)
        elif len(brack) == 2:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func, xa=brack[0], xb=brack[1])
        elif len(brack) == 3:
            xa, xb, xc = brack
            if xa > xc:  # swap so xa < xc can be assumed
                dum = xa
                xa = xc
                xc = dum
            assert (xa < xb) and (xb < xc), "Not a bracketing interval."
            fa = func(xa)
            fb = func(xb)
            fc = func(xc)
            assert (fb < fa) and (fb < fc), "Not a bracketing interval."
            funcalls = 3
        else:
            raise ValueError("Bracketing interval must be length 2 or 3 sequence.")
        ### END core bracket_info code ###

        self.funcalls += funcalls
        return xa, xb, xc, fa, fb, fc

    def optimize(self):
        # set up for optimization
        func = self.func
        if self._brack_info is None:
            self.set_bracket(None)
        xa, xb, xc, fa, fb, fc = self._brack_info
        _mintol = self._mintol
        _cg = self._cg
        #################################
        # BEGIN CORE ALGORITHM
        # we are making NO CHANGES in this
        #################################
        x = w = v = xb
        fw = fv = fx = func(x)
        if xa < xc:
            a = xa
            b = xc
        else:
            a = xc
            b = xa
        deltax = 0.0
        funcalls = 1
        iter = 0
        while iter < self.maxiter:
            tol1 = self.tol * abs(x) + _mintol
            tol2 = 2.0 * tol1
            xmid = 0.5 * (a + b)
            if abs(x - xmid) < (tol2 - 0.5 * (b - a)):  # check for convergence
                break
            infinities_present = [f for f in [fw, fv, fx] if numpy.isposinf(f)]
            if infinities_present or (abs(deltax) <= tol1):
                if x >= xmid:
                    deltax = a - x  # do a golden section step
                else:
                    deltax = b - x
                rat = _cg * deltax
            else:  # do a parabolic step
                tmp1 = (x - w) * (fx - fv)
                tmp2 = (x - v) * (fx - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0 * (tmp2 - tmp1)
                if tmp2 > 0.0:
                    p = -p
                tmp2 = abs(tmp2)
                dx_temp = deltax
                deltax = rat
                # check parabolic fit
                if (
                    (p > tmp2 * (a - x))
                    and (p < tmp2 * (b - x))
                    and (abs(p) < abs(0.5 * tmp2 * dx_temp))
                ):
                    rat = p * 1.0 / tmp2  # if parabolic step is useful.
                    u = x + rat
                    if (u - a) < tol2 or (b - u) < tol2:
                        if xmid - x >= 0:
                            rat = tol1
                        else:
                            rat = -tol1
                else:
                    if x >= xmid:
                        deltax = a - x  # if it's not do a golden section step
                    else:
                        deltax = b - x
                    rat = _cg * deltax

            if abs(rat) < tol1:  # update by at least tol1
                if rat >= 0:
                    u = x + tol1
                else:
                    u = x - tol1
            else:
                u = x + rat
            fu = func(u)  # calculate new output value
            funcalls += 1

            if fu > fx:  # if it's bigger than current
                if u < x:
                    a = u
                else:
                    b = u
                if (fu <= fw) or (w == x):
                    v = w
                    w = u
                    fv = fw
                    fw = fu
                elif (fu <= fv) or (v == x) or (v == w):
                    v = u
                    fv = fu
            else:
                if u >= x:
                    a = x
                else:
                    b = x
                v = w
                w = x
                x = u
                fv = fw
                fw = fx
                fx = fu

            iter += 1
        #################################
        # END CORE ALGORITHM
        #################################

        self.xmin = x
        self.fval = fx
        self.iter = iter
        self.funcalls = funcalls

    def get_result(self, full_output=False):
        if full_output:
            return self.xmin, self.fval, self.iter, self.funcalls
        else:
            return self.xmin


def brent(func, brack=None, tol=1.48e-8, full_output=0, maxiter=500):
    """Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol.

    :Parameters:

        func : callable f(x)
            Objective function.
        brack : tuple
            Triple (a,b,c) where (a<b<c) and func(b) <
            func(a),func(c).  If bracket consists of two numbers (a,c)
            then they are assumed to be a starting interval for a
            downhill bracket search (see `bracket`); it doesn't always
            mean that the obtained solution will satisfy a<=x<=c.
        full_output : bool
            If True, return all output args (xmin, fval, iter,
            funcalls).

    :Returns:

        xmin : ndarray
            Optimum point.
        fval : float
            Optimum value.
        iter : int
            Number of iterations.
        funcalls : int
            Number of objective function evaluations made.

    Notes
    -----

    Uses inverse parabolic interpolation when possible to speed up convergence
    of golden section method.

    """
    brent = Brent(func=func, tol=tol, maxiter=maxiter)
    brent.set_bracket(brack)
    brent.optimize()
    return brent.get_result(full_output=full_output)


def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
    """Given a function and distinct initial points, search in the
    downhill direction (as defined by the initital points) and return
    new points xa, xb, xc that bracket the minimum of the function
    f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
    solution will satisfy xa<=x<=xb

    :Parameters:

        func : callable f(x,*args)
            Objective function to minimize.
        xa, xb : float
            Bracketing interval.
        args : tuple
            Additional arguments (if present), passed to `func`.
        grow_limit : float
            Maximum grow limit.
        maxiter : int
            Maximum number of iterations to perform.

    :Returns: xa, xb, xc, fa, fb, fc, funcalls

        xa, xb, xc : float
            Bracket.
        fa, fb, fc : float
            Objective function values in bracket.
        funcalls : int
            Number of function evaluations made.

    """
    _gold = 1.618034
    _verysmall_num = 1e-21
    fa = func(*(xa,) + args)
    fb = func(*(xb,) + args)
    if fa < fb:  # Switch so fa > fb
        dum = xa
        xa = xb
        xb = dum
        dum = fa
        fa = fb
        fb = dum
    xc = xb + _gold * (xb - xa)
    fc = func(*((xc,) + args))
    funcalls = 3
    iter = 0
    while fc < fb:
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        if abs(val) < _verysmall_num:
            denom = 2.0 * _verysmall_num
        else:
            denom = 2.0 * val
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
        wlim = xb + grow_limit * (xc - xb)
        if iter > maxiter:
            raise RuntimeError("Too many iterations.")
        iter += 1
        if (w - xc) * (xb - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if fw < fc:
                xa = xb
                xb = w
                fa = fb
                fb = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif fw > fb:
                xc = w
                fc = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim) * (wlim - xc) >= 0.0:
            w = wlim
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim) * (xc - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if fw < fc:
                xb = xc
                xc = w
                w = xc + _gold * (xc - xb)
                fb = fc
                fc = fw
                fw = func(*((w,) + args))
                funcalls += 1
        else:
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        xa = xb
        xb = xc
        xc = w
        fa = fb
        fb = fc
        fc = fw
    return xa, xb, xc, fa, fb, fc, funcalls


def _linesearch_powell(linesearch, func, p, xi, tol):
    """Line-search algorithm using fminbound.

    Find the minimium of the function ``func(x0+ alpha*direc)``.

    """

    def myfunc(alpha):
        return func(p + alpha * xi)

    alpha_min, fret, iter, num = linesearch(myfunc, full_output=1, tol=tol)
    xi = alpha_min * xi
    return squeeze(fret), p + xi, xi


def fmin_powell(
    func,
    x0,
    args=(),
    xtol=1e-4,
    ftol=1e-4,
    maxiter=None,
    maxfun=None,
    full_output=0,
    disp=1,
    retall=0,
    callback=None,
    direc=None,
    linesearch=brent,
):
    """Minimize a function using modified Powell's method.

    :Parameters:

      func : callable f(x,*args)
          Objective function to be minimized.
      x0 : ndarray
          Initial guess.
      args : tuple
          Eextra arguments passed to func.
      callback : callable
          An optional user-supplied function, called after each
          iteration.  Called as ``callback(n,xk,f)``, where ``xk`` is the
          current parameter vector.
      direc : ndarray
          Initial direction set.

    :Returns: (xopt, {fopt, xi, direc, iter, funcalls, warnflag}, {allvecs})

        xopt : ndarray
            Parameter which minimizes `func`.
        fopt : number
            Value of function at minimum: ``fopt = func(xopt)``.
        direc : ndarray
            Current direction set.
        iter : int
            Number of iterations.
        funcalls : int
            Number of function calls made.
        warnflag : int
            Integer warning flag:
                1 : Maximum number of function evaluations.
                2 : Maximum number of iterations.
        allvecs : list
            List of solutions at each iteration.

    *Other Parameters*:

      xtol : float
          Line-search error tolerance.
      ftol : float
          Absolute error in ``func(xopt)`` acceptable for convergence.
      maxiter : int
          Maximum number of iterations to perform.
      maxfun : int
          Maximum number of function evaluations to make.
      full_output : bool
          If True, fopt, xi, direc, iter, funcalls, and
          warnflag are returned.
      disp : bool
          If True, print convergence messages.
      retall : bool
          If True, return a list of the solution at each iteration.


    :Notes:

        Uses a modification of Powell's method to find the minimum of
        a function of N variables.

    """
    # we need to use a mutable object here that we can update in the
    # wrapper function
    fcalls, func = wrap_function(func, args)
    x = asarray(x0).flatten()
    if retall:
        allvecs = [x]
    N = len(x)
    rank = len(x.shape)
    if not -1 < rank < 2:
        raise ValueError("Initial guess must be a scalar or rank-1 sequence.")
    if maxiter is None:
        maxiter = N * 1000
    if maxfun is None:
        maxfun = N * 1000

    if direc is None:
        direc = eye(N, dtype=float)
    else:
        direc = asarray(direc, dtype=float)

    fval = squeeze(func(x))
    x1 = x.copy()
    iter = 0
    ilist = list(range(N))
    while True:
        fx = fval
        bigind = 0
        delta = 0.0
        for i in ilist:
            direc1 = direc[i]
            fx2 = fval
            fval, x, direc1 = _linesearch_powell(
                linesearch, func, x, direc1, xtol * 100
            )
            if (fx2 - fval) > delta:
                delta = fx2 - fval
                bigind = i
        iter += 1
        if callback is not None:
            callback(fcalls[0], x, fval, delta)
        if retall:
            allvecs.append(x)
        if abs(fx - fval) < ftol:
            break
        if fcalls[0] >= maxfun:
            break
        if iter >= maxiter:
            break

        # Construct the extrapolated point
        direc1 = x - x1
        x2 = 2 * x - x1
        x1 = x.copy()
        fx2 = squeeze(func(x2))

        if fx > fx2:
            t = 2.0 * (fx + fx2 - 2.0 * fval)
            temp = fx - fval - delta
            t *= temp * temp
            temp = fx - fx2
            t -= delta * temp * temp
            if t < 0.0:
                fval, x, direc1 = _linesearch_powell(
                    linesearch, func, x, direc1, xtol * 100
                )
                direc[bigind] = direc[-1]
                direc[-1] = direc1

    warnflag = 0
    if fcalls[0] >= maxfun:
        warnflag = 1
        if disp:
            print(
                "Warning: Maximum number of function evaluations has " "been exceeded."
            )
    elif iter >= maxiter:
        warnflag = 2
        if disp:
            print("Warning: Maximum number of iterations has been exceeded")
    else:
        if disp:
            print("Optimization terminated successfully.")
            print(f"         Current function value: {fval:f}")
            print("         Iterations: %d" % iter)
            print("         Function evaluations: %d" % fcalls[0])

    x = squeeze(x)

    if full_output:
        retlist = x, fval, direc, iter, fcalls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist

#!/usr/bin/env python

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

EPS = 1e-15


def bisection(func, a, b, args=(), xtol=1e-10, maxiter=400):
    """Bisection root-finding method.  Given a function and an interval with
    func(a) * func(b) < 0, find the root between a and b.
    """
    if b < a:
        (a, b) = (b, a)
    i = 1
    eva = func(a, *args)
    evb = func(b, *args)
    assert eva * evb < 0, "Must start with interval with func(a) * func(b) <0"
    while i <= maxiter:
        dist = (b - a) / 2.0
        p = a + dist
        if dist / max(1.0, abs(p)) < xtol:
            return p
        ev = func(p, *args)
        if ev == 0:
            return p
        i += 1
        if ev * eva > 0:
            a = p
        else:
            b = p
    raise RuntimeError("bisection failed after %d iterations." % maxiter)


def brent(func, a, b, args=(), xtol=1e-10, maxiter=100):
    """Fast and robust root-finding method.  Given a function and an
    interval with func(a) * func(b) < 0, find the root between a and b.

    From Numerical Recipes
    """
    if b < a:
        (a, b) = (b, a)
    i = 1
    fa = func(a, *args)
    fb = func(b, *args)
    assert fa * fb < 0, "Must start with interval with func(a) * func(b) <0"
    (c, fc) = (b, fb)
    while i <= maxiter:
        if fb * fc > 0.0:
            (c, fc) = (a, fa)
            d = e = b - a
        if abs(fc) < abs(fb):
            (a, fa) = (b, fb)
            (b, fb) = (c, fc)
            (c, fc) = (a, fa)
        tol1 = 2.0 * EPS * abs(b) + 0.5 * xtol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -1.0 * q
            p = abs(p)
            min1 = 3.0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)
            if 2.0 * p < min(min1, min2):
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d
        (a, fa) = (b, fb)
        if abs(d) > tol1:
            b += d
        elif xm < 0.0:
            b -= tol1
        else:
            b += tol1
        fb = func(b, *args)
        i += 1
    raise RuntimeError("solver failed after %d iterations." % maxiter)


def find_root(func, x, direction, bound, xtol=None, expected_exception=None):
    if xtol is None:
        xtol = 1e-10

    def sign_func(z):
        # +ve if f(z) is +ve
        # zero if f(z) is -ve (what we want)
        # -ve if f(z) causes an error
        try:
            y = func(z)
            if y < 0:
                return 0
            else:
                return 1
        except expected_exception:
            return -1

    # Bracket root
    # Start out ignoring the bound as that is likely to be an error-
    # prone part of the range, and if there are multiple roots we want the
    # one closest to x.
    x_range = abs(bound - x)
    max_delta = x_range / 5
    if not max_delta:
        return None
    delta = min(0.01, max_delta)
    assert func(x) > 0
    x2 = x
    while 1:
        x1 = x2
        delta = min(delta * 2, max_delta)
        x2 = x1 + direction * delta
        if direction * (x2 - bound) > 0:
            x2 = bound
        y = sign_func(x2)
        if y <= 0 or x2 == bound:
            break

    # Hit a bound (or error).
    # Look for -ve between the +ve x1 and error x2
    if y == -1:
        x2 = bisection(sign_func, x1, x2, xtol=max(xtol, 1e-5))
        y = sign_func(x2)

    if y != 0:
        return None

    return brent(func, x1, x2, xtol=xtol)

from math import fsum

from numpy import array, diag, diagonal, dot, eye, isclose, log, sqrt
from numpy.linalg import slogdet
from numpy.testing import assert_allclose, assert_equal

import cogent3.util.misc

from cogent3.maths.util import safe_p_log_p, validate_freqs_array


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def paralinear_discrete_time(P, pi, validate=False):
    """
    Parameters
    ----------
    P : numpy array
        row stochastic matrix
    pi : numpy array
        row vector state frequencies
    validate : bool
    Returns
    -------
    the paralinear metric of Lake (1994) for a discrete time process
    """
    if validate:
        assert_equal(P.shape[0], pi.shape[0], err_msg="pi mismatched shape")
        assert pi.ndim == 1, "pi has incorrect dimension"
        assert_allclose(P.sum(axis=1), 1, err_msg="invalid P")
        assert_allclose(pi.sum(), 1, err_msg="invalid pi")

    pi_end = diag(dot(pi, P))
    pi = diag(pi)
    J = dot(pi, P)
    sign, b = slogdet(J)
    b *= sign
    b *= -1
    sign, a = slogdet(pi)
    a *= sign / 2

    sign, c = slogdet(pi_end)
    c *= sign / 2
    return b + a + c


def paralinear_continuous_time(P, pi, Q, validate=False):
    """
    Parameters
    ----------
    P : numpy array
        row stochastic matrix
    pi : numpy array
        row vector state frequencies
    Q : numpy array
        rate-matrix, rows sum to zero, off-diagnoal >= 0. If None, returns
        paralinear_discrete_time(P, pi, validate)
    validate : bool
    Returns
    -------
    the paralinear metric of Lake (1994)
    """
    if validate:
        assert_equal(Q.shape, P.shape, err_msg="Q/P mismatched shape")
        assert_equal(Q.shape[0], pi.shape[0], err_msg="pi mismatched shape")
        assert pi.ndim == 1, "pi has incorrect dimension"
        assert_allclose(Q.sum(axis=1), 0, err_msg="invalid Q", atol=1e-12)
        assert_allclose(P.sum(axis=1), 1, err_msg="invalid P")
        assert_allclose(pi.sum(), 1, err_msg="invalid pi")
        off_diag = ~eye(Q.shape[0])
        assert not (Q[off_diag] < 0).any(), "invalid Q"

    # todo need to implement a safe natural log to handle possible 0 elements
    a = -log(pi).sum() / 2
    b = -diagonal(Q).sum()
    c = log(dot(pi, P)).sum() / 2
    pl = a + b + c
    pl = getattr(pl, "real", pl)
    return pl


def jsd(*vectors, validate=False):
    """calculate Jensen–Shannon divergence between two or more probability distributions

    Parameters
    ----------
    *vectors
        >= 2 frequency vectors, each vector must sum to 1
    validate : bool
        check the consistency of the provided vectors, namely they're
        all 1D, have the same number of elements and sum to 1
    """
    # Convert input arrays into numpy arrays
    vectors = [array(v) for v in vectors]
    num_vectors = len(vectors)
    if validate:
        assert len({f.shape for f in vectors}) == 1, "mismatched shape"
        dims = {v.ndim for v in vectors}
        assert dims == {1}, "incorrect dimension"
        try:
            for v in vectors:
                validate_freqs_array(v)
        except ValueError as err:
            raise AssertionError("freqs not valid") from err

    H_mn = fsum(safe_p_log_p(array(vectors).mean(axis=0)))
    mn_H = fsum([fsum(i) for i in map(safe_p_log_p, vectors)]) / num_vectors
    jsd_ = H_mn - mn_H
    if jsd_ < 0:
        if not isclose(jsd_, 0, atol=1e-10):
            raise ArithmeticError(
                f"{jsd_} is negative and below defined precision threshold"
            )
        jsd_ = 0

    return jsd_


@cogent3.util.misc.extend_docstring_from(jsd)
def jsm(*args, **kwargs):
    """
    Returns
    -------
    The square root of the Jensen–Shannon divergence,
    which is a metric often referred to as Jensen-Shannon distance (2003)
    """
    val = jsd(*args, **kwargs)
    return sqrt(val)

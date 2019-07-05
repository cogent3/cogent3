from numpy import diagonal, dot, eye
from numpy.testing import assert_allclose, assert_equal

from .util import safe_log


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def paralinear(Q, P, pi, validate=False):
    """
    Parameters
    ----------
    Q : numpy array
        rate-matrix, rows sum to zero, off-diagnoal >= 0
    P : numpy array
        row stochastic matrix
    pi : numpy array
        row vector state frequencies
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

    a = -safe_log(pi).sum() / 2
    b = -diagonal(Q).sum()
    c = safe_log(dot(pi, P)).sum() / 2
    pl = a + b + c
    pl = getattr(pl, "real", pl)
    return pl

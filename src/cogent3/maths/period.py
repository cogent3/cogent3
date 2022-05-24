from numpy import (
    arange,
    array,
    cos,
    exp,
    fft,
    float64,
    multiply,
    pi,
    polyval,
    sqrt,
    zeros,
)

from .period_numba import autocorr_inner, goertzel_inner, ipdft_inner


__author__ = "Hua Ying, Julien Epps and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Julien Epps", "Hua Ying", "Gavin Huttley", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def _goertzel_inner(x, N, period):
    coeff = 2.0 * cos(2 * pi / period)
    s_prev = 0.0
    s_prev2 = 0.0
    for n in range(N):
        s = x[n] + coeff * s_prev - s_prev2
        s_prev2 = s_prev
        s_prev = s
    return sqrt(s_prev2 ** 2 + s_prev ** 2 - coeff * s_prev2 * s_prev)


def _ipdft_inner(x, X, W, ulim, N):  # naive python
    for p in range(ulim):
        w = 1
        for n in range(N):
            if n != 0:
                w *= W[p]
            X[p] = X[p] + x[n] * w
    return X


def _ipdft_inner2(x, X, W, ulim, N):  # fastest python
    p = x[::-1]  # reversed
    X = polyval(p, W)
    return X


def _autocorr_inner2(x, xc, N):  # fastest python
    products = multiply.outer(x, x)
    v = [products.trace(offset=m) for m in range(-len(x) + 1, len(x))]
    xc.put(range(xc.shape[0]), v)


def _autocorr_inner(x, xc, N):  # naive python
    for m in range(-N + 1, N):
        for n in range(N):
            if 0 <= n - m < N:
                xc[m + N - 1] += x[n] * x[n - m]


def goertzel(x, period):
    """returns the array(power), array(period) from series x for period
    result objects are arrays for consistency with that other period
    estimation functions"""
    calc = Goertzel(len(x), period=period)
    return calc(x)


class _PeriodEstimator(object):
    """parent class for period estimation"""

    def __init__(self, length, llim=None, ulim=None, period=None):
        super(_PeriodEstimator, self).__init__()
        self.length = length
        self.llim = llim or 2
        self.ulim = ulim or (length - 1)

        if self.ulim > length:
            raise RuntimeError("Error: ulim > length")

        self.period = period

    def get_num_stats(self):
        """returns the number of statistics computed by this calculator"""
        return 1


class AutoCorrelation(_PeriodEstimator):
    """class for repetitive calculation of autocorrelation for series of
    fixed length"""

    def __init__(self, length, llim=None, ulim=None, period=None):
        """
        Parameters
        ----------
        length : int
            the signal length
        llim : int
            lower limit
        ulim : int
            upper limit
        period : int
            a specific period to return the IPDFT power for

        Notes
        -----
        If ``x = [1,1,1,1]`, ``xc = [1,2,3,4,3,2,1]``

        The middle element of xc corresponds to a lag (period) of 0
        ``xc`` is always symmetric for real ``x``
        """
        super(AutoCorrelation, self).__init__(length, llim, ulim, period)

        periods = list(range(-length + 1, length))

        self.min_idx = periods.index(self.llim)
        self.max_idx = periods.index(self.ulim)
        self.periods = array(periods[self.min_idx : self.max_idx + 1])
        self.xc = zeros(2 * self.length - 1)

    def __call__(self, x):
        x = array(x, float64)
        self.xc.fill(0.0)
        autocorr_inner(x, self.xc, self.length)
        xc = self.xc[self.min_idx : self.max_idx + 1]
        if self.period is not None:
            return xc[self.period - self.llim]

        return xc, self.periods


def auto_corr(x, llim=None, ulim=None):
    """returns the autocorrelation of x
    e.g. if x = [1,1,1,1], xc = [1,2,3,4,3,2,1]
    The middle element of xc corresponds to a lag (period) of 0
    xc is always symmetric for real x
    N is the length of x
    """
    _autocorr = AutoCorrelation(len(x), llim=llim, ulim=ulim)
    return _autocorr(x)


class Ipdft(_PeriodEstimator):
    def __init__(self, length, llim=None, ulim=None, period=None, abs_ft_sig=True):
        """factory function for computing the integer period discrete Fourier
        transform for repeated application to signals of the same length.

        Parameters
        ----------
        length : int
            the signal length
        llim : int
            lower limit
        ulim : int
            upper limit
        period : int
            a specific period to return the IPDFT power for
        abs_ft_sig : bool
            if True, returns absolute value of signal
        """
        if period is not None:
            llim = period
            ulim = period
        super(Ipdft, self).__init__(length, llim, ulim, period)
        self.periods = array(list(range(self.llim, self.ulim + 1)))
        self.W = exp(-1j * 2 * pi / arange(1, self.ulim + 1))
        self.X = array([0 + 0j] * self.length)
        self.abs_ft_sig = abs_ft_sig

    def __call__(self, x):
        x = array(x, float64)
        self.X.fill(0 + 0j)
        self.X = ipdft_inner(x, self.X, self.W, self.ulim, self.length)
        pwr = self.X[self.llim - 1 : self.ulim]

        if self.abs_ft_sig:
            pwr = abs(pwr)

        if self.period is not None:
            return pwr[self.period - self.llim]

        return array(pwr), self.periods


class Goertzel(_PeriodEstimator):
    """Computes the power of a signal for a specific period"""

    def __init__(self, length=None, period=None, **kwargs):
        assert period is not None, "Goertzel requires a period"
        super(Goertzel, self).__init__(length=length, period=period)

    def __call__(self, x):
        x = array(x, float64)
        return goertzel_inner(x, self.length, self.period)


class Hybrid(_PeriodEstimator):
    """hybrid statistic and corresponding periods for signal x

    See Epps. EURASIP Journal on Bioinformatics and Systems Biology, 2009"""

    def __init__(
        self,
        length,
        llim=None,
        ulim=None,
        period=None,
        abs_ft_sig=True,
        return_all=False,
    ):
        """
        Parameters
        ----------
        length : int
            the signal length
        llim : int
            lower limit
        ulim : int
            upper limit
        period : int
            a specific period to return the IPDFT power for
        abs_ft_sig : bool
            if True, returns absolute value of signal
        return_all : bool
            whether to return the hybrid, ipdft, autocorr statistics as
            a numpy array, or just the hybrid statistic
        """
        super(Hybrid, self).__init__(length, llim, ulim, period)
        self.ipdft = Ipdft(length, llim, ulim, period, abs_ft_sig)
        self.auto = AutoCorrelation(length, llim, ulim, period)
        self._return_all = return_all

    def get_num_stats(self):
        """the number of stats computed by this calculator"""
        return [1, 3][self._return_all]

    def __call__(self, x):
        if self.period is None:
            auto_sig, auto_periods = self.auto(x)
            ft_sig, ft_periods = self.ipdft(x)
            hybrid = auto_sig * ft_sig
            if self._return_all:
                result = array([hybrid, ft_sig, auto_sig]), ft_periods
            else:
                result = hybrid, ft_periods
        else:
            auto_sig = self.auto(x)
            # note that goertzel(x, period) performance is slower than ipdft!
            ft_sig = self.ipdft(x)
            hybrid = auto_sig * ft_sig
            if self._return_all:
                result = array([abs(hybrid), ft_sig, auto_sig])
            else:
                result = abs(hybrid)
        return result


def ipdft(x, llim=None, ulim=None, period=None):
    """returns the integer period discrete Fourier transform of the signal x

    Parameters
    ----------
    x
        series of symbols
    llim : int
        lower limit
    ulim : int
        upper limit
    period : int
        a specific period to return the IPDFT power for
    """
    x = array(x, float64)
    ipdft_calc = Ipdft(len(x), llim, ulim, period)
    return ipdft_calc(x)


def hybrid(x, llim=None, ulim=None, period=None, return_all=False):
    """
    Return hybrid statistic and corresponding periods for signal x

    Parameters
    ----------
    x
        series of symbols
    llim : int
        lower limit
    ulim : int
        upper limit
    period : int
        a specific period to return the IPDFT power for
    return_all : bool
        whether to return the hybrid, ipdft, autocorr statistics as
        a numpy array, or just the hybrid statistic

    See Epps. EURASIP Journal on Bioinformatics and Systems Biology, 2009, 9
    """
    hybrid_calc = Hybrid(len(x), llim, ulim, period, return_all=return_all)
    x = array(x, float)
    return hybrid_calc(x)


def dft(x, **kwargs):
    """
    Return discrete fft and corresponding periods for signal x
    """
    n = len(x) // 2 * 2
    x = array(x[:n])
    pwr = fft.rfft(x, n)[1:]
    freq = (arange(n / 2 + 1) / (float(n)))[1:]
    pwr = list(pwr)
    periods = [1 / f for f in freq]
    pwr.reverse()
    periods.reverse()
    return array(pwr), array(periods)

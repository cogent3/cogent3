import numpy


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def aic(lnL, nfp, sample_size=None):
    """returns Aikake Information Criterion

    Parameters
    ----------
    lnL
        the maximum log
    nfp
        the number of free parameters in the model
    sample_size
        if provided, the second order AIC is returned

    """
    if sample_size is None:
        correction = 1
    else:
        assert sample_size > 0, f"Invalid sample_size {sample_size}"
        correction = sample_size / (sample_size - nfp - 1)

    return -2 * lnL + 2 * nfp * correction


def bic(lnL, nfp, sample_size):
    """returns Bayesian Information Criterion

    Parameters
    ----------
    lnL
        the maximum log
    nfp
        the number of free parameters in the model
    sample_size
        size of the sample

    """
    return -2 * lnL + nfp * numpy.log(sample_size)

import numpy


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

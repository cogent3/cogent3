***************************
Estimating periodic signals
***************************

.. sectionauthor:: Gavin Huttley, Julien Epps, Hua Ying

We consider two different scenarios:

- estimating the periods in a signal
- estimating the power for a given period
- measuring statistical significance for the latter case

Estimating the periods in a signal
==================================

For numerical (continuous) data
-------------------------------

We first make some sample data. A periodic signal and some noise.

..
    We set a seed for the random number generator so that we can get
    consistent generation of the same series. This makes the document
    robust for doctesting.

.. doctest::
    :hide:
    
    >>> import numpy
    >>> numpy.random.seed(11)

.. doctest::
    
    >>> import numpy
    >>> t = numpy.arange(0, 10, 0.1)
    >>> n = numpy.random.randn(len(t))
    >>> nse = numpy.convolve(n, numpy.exp(-t/0.05))*0.1
    >>> nse = nse[:len(t)]
    >>> sig = numpy.sin(2*numpy.pi*t) + nse

Discrete Fourier transform
^^^^^^^^^^^^^^^^^^^^^^^^^^

We now use the discrete Fourier transform to estimate periodicity in this signal. Given we set the period to equal 10, we expect the maximum power for that index.

.. doctest::
    
    >>> from cogent.maths.period import dft
    >>> pwr, period = dft(sig)
    >>> print period
    [   2.            2.04081633    2.08333333    2.12765957    2.17391304
        2.22222222    2.27272727    2.3255814     2.38095238    2.43902439
        2.5           2.56410256    2.63157895    2.7027027     2.77777778
        2.85714286    2.94117647    3.03030303    3.125         3.22580645...
    >>> print pwr
    [ 1.06015801 +0.00000000e+00j  0.74686707 -1.93971914e-02j
      0.36784793 -2.66370366e-02j  0.04384413 +2.86970840e-02j
      1.54473269 -2.43777386e-02j  0.28522968 -2.33602932e-01j...

The power (``pwr``) is returned as an array of complex numbers, so we convert into real numbers using ``abs``. We then zip the power and corresponding periods and sort to identify the period with maximum signal.

    >>> pwr = abs(pwr)
    >>> max_pwr, max_period = sorted(zip(pwr,period))[-1]
    >>> print max_pwr, max_period
    50.7685934719 10.0

Auto-correlation
^^^^^^^^^^^^^^^^

We now use auto-correlation.

.. doctest::
    
    >>> from cogent.maths.period import auto_corr
    >>> pwr, period = auto_corr(sig)
    >>> print period
    [ 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24...
    >>> print pwr
    [  1.63366075e+01  -1.47309007e+01  -3.99310414e+01  -4.94779387e+01...

We then zip the power and corresponding periods and sort to identify the period with maximum signal.

.. doctest::
    
    >>> max_pwr, max_period = sorted(zip(pwr,period))[-1]
    >>> print max_pwr, max_period
    46.7917300733 10

For symbolic data
-----------------

We create a sequence as just a string

.. doctest::
    
    >>> s = 'ATCGTTGGGACCGGTTCAAGTTTTGGAACTCGCAAGGGGTGAATGGTCTTCGTCTAACGCTGG'\
    ...     'GGAACCCTGAATCGTTGTAACGCTGGGGTCTTTAACCGTTCTAATTTAACGCTGGGGGGTTCT'\
    ...     'AATTTTTAACCGCGGAATTGCGTC'

We then specify the motifs whose occurrences will be converted into 1, with all other motifs converted into 0. As we might want to do this in batches for many sequences we use a factory function.

.. doctest::
    
    >>> from cogent.maths.stats.period import SeqToSymbols
    >>> seq_to_symbols = SeqToSymbols(['AA', 'TT', 'AT'])
    >>> symbols = seq_to_symbols(s)
    >>> len(symbols) == len(s)
    True
    >>> symbols
    array([1, 0, 0, 0, 1, 0, 0, 0, 0, 0...

We then estimate the integer discrete Fourier transform for the full data. To do this, we need to pass in the symbols from full conversion of the sequence. The returned values are the powers and periods.

.. doctest::
    
    >>> from cogent.maths.period import ipdft
    >>> powers, periods = ipdft(symbols)
    >>> powers #doctest: +SKIP
    array([  3.22082108e-14,   4.00000000e+00,   9.48683298e+00,
             6.74585634e+00,   3.46410162e+00,   3.20674669e+00,...
    >>> periods
    array([  2,   3,   4...

We can also compute the auto-correlation statistic, and the hybrid (which combines IPDFT and auto-correlation).

.. doctest::
    
    >>> from cogent.maths.period import auto_corr, hybrid
    >>> powers, periods = auto_corr(symbols)
    >>> powers
    array([ 11.,   9.,  11.,   9.,   6...
    >>> periods
    array([  2,   3,   4...
    >>> powers, periods = hybrid(symbols)
    >>> powers #doctest: +SKIP
    array([  3.54290319e-13,   3.60000000e+01,   1.04355163e+02,
             6.07127071e+01,   2.07846097e+01,   2.88607202e+01,...
    >>> periods
    array([  2,   3,   4...

Estimating power for specified period
=====================================

For numerical (continuous) data
-------------------------------

We just use ``sig`` created above. The Goertzel algorithm gives the same result as the ``dft``.

.. doctest::
    
    >>> from cogent.maths.period import goertzel
    >>> pwr = goertzel(sig, 10)
    >>> print pwr
    50.7685934719

For symbolic data
-----------------

.. take example above and show how to compute it using autocorrelation

We use the symbols from the above example. For the ``ipdft``, ``auto_corr`` and ``hybrid`` functions we just need to identify the array index containing the period of interest and slice the corresponding value from the returned powers. The reported periods start at ``llim``, which defaults to 2, but indexes start at 0, the index for a period-5 is simply 5-``llim``.

.. doctest::
    
    >>> powers, periods = auto_corr(symbols)
    >>> llim = 2
    >>> period5 = 5-llim
    >>> periods[period5]
    5
    >>> powers[period5]
    9.0

For Fourier techniques, we can compute the power for a specific period more efficiently using Goertzel algorithm.

.. doctest::
    
    >>> from cogent.maths.period import goertzel
    >>> period = 4
    >>> power = goertzel(symbols, period)
    >>> ipdft_powers, periods = ipdft(symbols)
    >>> ipdft_power = abs(ipdft_powers[period-llim])
    >>> round(power, 6) == round(ipdft_power, 6)
    True
    >>> power
    9.4868...

It's also possible to specify a period to the stand-alone functions. As per the ``goertzel`` function, just the power is returned.

.. doctest::
    
    >>> power = hybrid(symbols, period=period)
    >>> power
    104.355...

Measuring statistical significance of periodic signals
======================================================

For numerical (continuous data)
-------------------------------

We use the signal provided above. Because significance testing is being done using a resampling approach, we define a calculator which precomputes some values to improve compute performance. For a continuous signal, we'll use the Goertzel algorithm.

.. doctest::
    
    >>> from cogent.maths.period import Goertzel
    >>> goertzel_calc = Goertzel(len(sig), period=10)

Having defined this, we then just pass this calculator to the ``blockwise_bootstrap`` function. The other critical settings are the ``block_size`` which specifies the size of segments of contiguous sequence positions to use for sampling and ``num_reps`` which is the number of permuted replicate sequences to generate.

.. doctest::
    
    >>> from cogent.maths.stats.period import blockwise_bootstrap
    >>> obs_stat, p = blockwise_bootstrap(sig, calc=goertzel_calc, block_size=10,
    ...                              num_reps=1000)
    >>> print obs_stat
    50.7685934719
    >>> print p
    0.0

For symbolic data
-----------------

Permutation testing
^^^^^^^^^^^^^^^^^^^

The very notion of permutation testing for periods, applied to a genome, requires the compute performance be as quick as possible. This means providing as much information up front as possible. We have made the implementation flexible by not assuming how the user will convert sequences to symbols. It's also the case that numerous windows of exactly the same size are being assessed. Accordingly, we use a class to construct a fixed signal length evaluator. We do this for the hybrid metric first.

.. doctest::
    
    >>> from cogent.maths.period import Hybrid
    >>> len(s)
    150
    >>> hybrid_calculator = Hybrid(len(s), period = 4)

.. note:: We defined the period length of interest in defining this calculator because we're interested in dinucleotide motifs.

We then construct a seq-to-symbol convertor.

.. doctest::
    
    >>> from cogent.maths.stats.period import SeqToSymbols
    >>> seq_to_symbols = SeqToSymbols(['AA', 'TT', 'AT'], length=len(s))

The rest is as per the analysis using ``Goertzel`` above.

.. doctest::
    
    >>> from cogent.maths.stats.period import blockwise_bootstrap
    >>> stat, p = blockwise_bootstrap(s, calc=hybrid_calculator,
    ...      block_size=10, num_reps=1000, seq_to_symbols=seq_to_symbols)
    ...     
    >>> print stat
    104.35...
    >>> p < 0.1
    True

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

.. jupyter-execute::
    :hide-code:

    import numpy

    numpy.random.seed(11)

.. jupyter-execute::

    import numpy

    t = numpy.arange(0, 10, 0.1)
    n = numpy.random.randn(len(t))
    nse = numpy.convolve(n, numpy.exp(-t / 0.05)) * 0.1
    nse = nse[: len(t)]
    sig = numpy.sin(2 * numpy.pi * t) + nse

Discrete Fourier transform
^^^^^^^^^^^^^^^^^^^^^^^^^^

We now use the discrete Fourier transform to estimate periodicity in this signal. Given we set the period to equal 10, we expect the maximum power for that index.

.. jupyter-execute::

    from cogent3.maths.period import dft

    pwr, period = dft(sig)
    print(period)
    print(pwr)

The power (``pwr``) is returned as an array of complex numbers, so we convert into real numbers using ``abs``. We then zip the power and corresponding periods and sort to identify the period with maximum signal.

.. jupyter-execute::

    pwr = abs(pwr)
    max_pwr, max_period = sorted(zip(pwr, period))[-1]
    print(max_pwr, max_period)

Auto-correlation
^^^^^^^^^^^^^^^^

We now use auto-correlation.

.. jupyter-execute::

    from cogent3.maths.period import auto_corr

    pwr, period = auto_corr(sig)
    print(period)
    print(pwr)

We then zip the power and corresponding periods and sort to identify the period with maximum signal.

.. jupyter-execute::

    max_pwr, max_period = sorted(zip(pwr, period))[-1]
    print(max_pwr, max_period)

For symbolic data
-----------------

We create a sequence as just a string

.. jupyter-execute::

    s = (
        "ATCGTTGGGACCGGTTCAAGTTTTGGAACTCGCAAGGGGTGAATGGTCTTCGTCTAACGCTGG"
        "GGAACCCTGAATCGTTGTAACGCTGGGGTCTTTAACCGTTCTAATTTAACGCTGGGGGGTTCT"
        "AATTTTTAACCGCGGAATTGCGTC"
    )

We then specify the motifs whose occurrences will be converted into 1, with all other motifs converted into 0. As we might want to do this in batches for many sequences we use a factory function.

.. jupyter-execute::

    from cogent3.maths.stats.period import SeqToSymbols

    seq_to_symbols = SeqToSymbols(["AA", "TT", "AT"])
    symbols = seq_to_symbols(s)
    len(symbols) == len(s)
    symbols

We then estimate the integer discrete Fourier transform for the full data. To do this, we need to pass in the symbols from full conversion of the sequence. The returned values are the powers and periods.

.. jupyter-execute::

    from cogent3.maths.period import ipdft

    powers, periods = ipdft(symbols)
    powers

.. jupyter-execute::

    periods

We can also compute the auto-correlation statistic, and the hybrid (which combines IPDFT and auto-correlation).

.. jupyter-execute::

    from cogent3.maths.period import auto_corr, hybrid

    powers, periods = auto_corr(symbols)
    powers

.. jupyter-execute::

    periods

.. jupyter-execute::

    powers, periods = hybrid(symbols)
    powers

.. jupyter-execute::

    periods

Estimating power for specified period
=====================================

For numerical (continuous) data
-------------------------------

We just use ``sig`` created above. The Goertzel algorithm gives the same result as the ``dft``.

.. jupyter-execute::

    from cogent3.maths.period import goertzel

    pwr = goertzel(sig, 10)
    print(pwr)

For symbolic data
-----------------

.. take example above and show how to compute it using autocorrelation

We use the symbols from the above example. For the ``ipdft``, ``auto_corr`` and ``hybrid`` functions we just need to identify the array index containing the period of interest and slice the corresponding value from the returned powers. The reported periods start at ``llim``, which defaults to 2, but indexes start at 0, the index for a period-5 is simply 5-``llim``.

.. jupyter-execute::

    powers, periods = auto_corr(symbols)
    llim = 2
    period5 = 5 - llim
    periods[period5]

.. jupyter-execute::

    powers[period5]

For Fourier techniques, we can compute the power for a specific period more efficiently using Goertzel algorithm.

.. jupyter-execute::

    from cogent3.maths.period import goertzel

    period = 4
    power = goertzel(symbols, period)
    ipdft_powers, periods = ipdft(symbols)
    ipdft_power = abs(ipdft_powers[period - llim])
    round(power, 6) == round(ipdft_power, 6)
    power

It's also possible to specify a period to the stand-alone functions. As per the ``goertzel`` function, just the power is returned.

.. jupyter-execute::

    power = hybrid(symbols, period=period)
    power

Measuring statistical significance of periodic signals
======================================================

For numerical (continuous data)
-------------------------------

We use the signal provided above. Because significance testing is being done using a resampling approach, we define a calculator which precomputes some values to improve compute performance. For a continuous signal, we'll use the Goertzel algorithm.

.. jupyter-execute::

    from cogent3.maths.period import Goertzel

    goertzel_calc = Goertzel(len(sig), period=10)

Having defined this, we then just pass this calculator to the ``blockwise_bootstrap`` function. The other critical settings are the ``block_size`` which specifies the size of segments of contiguous sequence positions to use for sampling and ``num_reps`` which is the number of permuted replicate sequences to generate.

.. jupyter-execute::

    from cogent3.maths.stats.period import blockwise_bootstrap

    obs_stat, p = blockwise_bootstrap(
        sig, calc=goertzel_calc, block_size=10, num_reps=1000
    )
    print(obs_stat)
    print(p)

For symbolic data
-----------------

Permutation testing
^^^^^^^^^^^^^^^^^^^

The very notion of permutation testing for periods, applied to a genome, requires the compute performance be as quick as possible. This means providing as much information up front as possible. We have made the implementation flexible by not assuming how the user will convert sequences to symbols. It's also the case that numerous windows of exactly the same size are being assessed. Accordingly, we use a class to construct a fixed signal length evaluator. We do this for the hybrid metric first.

.. jupyter-execute::

    from cogent3.maths.period import Hybrid

    hybrid_calculator = Hybrid(len(s), period=4)

.. note:: We defined the period length of interest in defining this calculator because we're interested in dinucleotide motifs.

We then construct a seq-to-symbol convertor.

.. jupyter-execute::

    from cogent3.maths.stats.period import SeqToSymbols

    seq_to_symbols = SeqToSymbols(["AA", "TT", "AT"], length=len(s))

The rest is as per the analysis using ``Goertzel`` above.

.. jupyter-execute::

    from cogent3.maths.stats.period import blockwise_bootstrap

    stat, p = blockwise_bootstrap(
        s,
        calc=hybrid_calculator,
        block_size=10,
        num_reps=1000,
        seq_to_symbols=seq_to_symbols,
    )
    print(stat)
    p < 0.1

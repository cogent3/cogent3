Performing a relative rate test
===============================

From cogent import all the components we need

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import HKY85
    >>> from cogent.maths import stats

Get your alignment and tree.

.. doctest::

    >>> al = LoadSeqs(filename = "data/test.paml")
    >>> t = LoadTree(filename = "data/test.tree")

Create a HKY85 model.

.. doctest::

    >>> sm = HKY85()

Make the controller object.

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(t)

Set the local clock for humans & Howler Monkey. This method is just a special interface to the more general setParamRules method.

.. doctest::

    >>> lf.setLocalClock("Human", "HowlerMon")

Get the likelihood function object this object performs the actual likelihood calculation.

.. doctest::

    >>> lf.setAlignment(al)

Optimise the function capturing the return optimised lnL, and parameter value vector.

.. doctest::

    >>> lf.optimise(show_progress = False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> lf.setName("clock")
    >>> print lf
    clock
    ======
     kappa
    ------
    4.8020
    ------
    =============================
         edge    parent    length
    -----------------------------
        Human    edge.0    0.0257
    HowlerMon    edge.0    0.0257
       edge.0    edge.1    0.0224
        Mouse    edge.1    0.2112
       edge.1      root    0.0000
    NineBande      root    0.0327
     DogFaced      root    0.0545
    -----------------------------
    ===============
    motif    mprobs
    ---------------
        T    0.1433
        C    0.1600
        A    0.3800
        G    0.3167
    ---------------

We extract the log-likelihood and number of free parameters for later use.

.. doctest::

    >>> null_lnL = lf.getLogLikelihood()
    >>> null_nfp = lf.getNumFreeParams()

Clear the local clock constraint, freeing up the branch lengths.

.. doctest::

    >>> lf.setParamRule('length', is_independent=True)

Run the optimiser capturing the return optimised lnL, and parameter value vector.

.. doctest::

    >>> lf.optimise(show_progress=False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> lf.setName("non clock")
    >>> print lf
    non clock
    ======
     kappa
    ------
    4.8027
    ------
    =============================
         edge    parent    length
    -----------------------------
        Human    edge.0    0.0347
    HowlerMon    edge.0    0.0167
       edge.0    edge.1    0.0224
        Mouse    edge.1    0.2112
       edge.1      root    0.0000
    NineBande      root    0.0327
     DogFaced      root    0.0545
    -----------------------------
    ===============
    motif    mprobs
    ---------------
        T    0.1433
        C    0.1600
        A    0.3800
        G    0.3167
    ---------------

These two lnL's are now used to calculate the likelihood ratio statistic it's degrees-of-freedom and the probability of observing the LR.

.. doctest::

    >>> LR = 2 * (lf.getLogLikelihood() - null_lnL)
    >>> df = lf.getNumFreeParams() - null_nfp
    >>> P = stats.chisqprob(LR, df)

Print this and look up a chi-sq with number of edges - 1 degrees of freedom.

.. doctest::

    >>> print "Likelihood ratio statistic = ", LR
    Likelihood ratio statistic =  0.34...
    >>> print "degrees-of-freedom = ", df
    degrees-of-freedom =  1
    >>> print "probability = ", P
    probability =  0.5...


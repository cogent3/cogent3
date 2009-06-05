A test of the neutral theory
============================

This file contains an example for performing a likelihood ratio test of neutrality. The test compares a model where the codon model parameter omega is constrained to be the same for all edges against one where each edge has its' own omega. From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import GY94
    >>> from cogent.maths import stats

Get your alignment and tree.

.. doctest::

    >>> al = LoadSeqs("data/test.paml")
    >>> t = LoadTree("data/test.tree")

We use a Goldman Yang 1994 model.

.. doctest::

    >>> sm = GY94()

Make the controller object

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setTablesFormat(digits=2,space=2)

Get the likelihood function object this object performs the actual likelihood calculation.

.. doctest::

    >>> lf.setAlignment(al)

By default, parameters other than branch lengths are treated as global in scope, so we don't need to do anything special here. We can influence how rigorous the optimisation will be, and switch between the global and local optimisers provided in the toolkit using arguments to the optimise method. The ``global_tolerance=1.0`` argument specifies conditions for an early break from simulated annealing which will be automatically followed by the Powell local optimiser. .. note:: the 'results' are of course nonsense.

.. doctest::

    >>> lf.optimise(global_tolerance = 1.0, show_progress=False)

View the resulting maximum-likelihood parameter values

.. doctest::

    >>> print lf
    Likelihood Function Table
    ============
    kappa  omega
    ------------
     9.28   1.87
    ------------
    =========================
         edge  parent  length
    -------------------------
        Human  edge.0    0.10
    HowlerMon  edge.0    0.05
       edge.0  edge.1    0.07
        Mouse  edge.1    0.91
       edge.1    root    0.00
    NineBande    root    0.11
     DogFaced    root    0.18...

We'll get the lnL and number of free parameters for later use.

.. doctest::

    >>> null_lnL = lf.getLogLikelihood()
    >>> null_nfp = lf.getNumFreeParams()

Specify each edge has it's own omega by just modifying the existing ``lf``. This means the new function will start with the above values.

.. doctest::

    >>> lf.setParamRule("omega", is_independent = True)

Optimise the likelihood function, this time just using the local optimiser.

.. doctest::

    >>> lf.optimise(local = True, show_progress=False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> print lf
    Likelihood Function Table
    =====
    kappa
    -----
     8.95
    -----
    =====================================
         edge  parent  length       omega
    -------------------------------------
        Human  edge.0    0.10   999999.98
    HowlerMon  edge.0    0.06   999999.99
       edge.0  edge.1    0.07   999999.99
        Mouse  edge.1    0.96        0.70
       edge.1    root    0.00   976875.85
    NineBande    root    0.11  1000000.00
     DogFaced    root    0.18        1.10...

.. note:: The parameter estimates for ``omega`` are highly implausible, reflecting (in this case) our small and uninformative data set.

Get out an annotated tree, it looks just like a tree, but has the maximum-likelihood parameter estimates attached to each tree edge. This object can be used for plotting, or to provide starting estimates to a related model.

.. doctest::

    >>> at = lf.getAnnotatedTree()

Get a dictionary of the statistics that I could use for post-processing.

.. doctest::

    >>> sd = lf.getStatisticsAsDict(with_edge_names=True)

The lnL's from the two models are now used to calculate the likelihood ratio statistic (``LR``) it's degrees-of-freedom (``df``) and the probability (``P``) of observing the LR.

.. doctest::

    >>> LR = 2 * (lf.getLogLikelihood() - null_lnL)
    >>> df = lf.getNumFreeParams() - null_nfp
    >>> P = stats.chisqprob(LR, df)

Print this and look up a chi-sq with number of edges - 1 degrees of freedom.

.. doctest::

    >>> print "Likelihood ratio statistic = ", LR
    Likelihood ratio statistic =  4.4...
    >>> print "degrees-of-freedom = ", df
    degrees-of-freedom =  6
    >>> print "probability = ", P
    probability =  0.6...


Performing a parametric bootstrap
=================================

.. sectionauthor:: Gavin Huttley

This file contains an example for estimating the probability of a Likelihood ratio statistic obtained from a relative rate test. The bootstrap classes can take advantage of parallel architectures.

From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve import bootstrap
    >>> from cogent.evolve.models import HKY85
    >>> from cogent.maths import stats

Define the null model that takes an alignment object and returns a likelihood function properly assembled for optimising the likelihood under the null hypothesis. The sample distribution is generated using this model.

We will use a HKY model.

.. doctest::

    >>> def create_alt_function():
    ...     t = LoadTree("data/test.tree")
    ...     sm = HKY85()
    ...     return sm.makeLikelihoodFunction(t)

Define a function that takes an alignment object and returns an appropriately assembled function for the alternative model. Since the two models are identical bar the constraint on the branch lengths, we'll use the same code to generate the basic likelihood function as for the alt model, and then apply the constraint here

.. doctest::

    >>> def create_null_function():
    ...     lf = create_alt_function()
    ...     # set the local clock for humans & howler monkey
    ...     lf.setLocalClock("Human", "HowlerMon")
    ...     return lf

Get our observed data alignment

.. doctest::

    >>> aln = LoadSeqs(filename = "data/long_testseqs.fasta")

Create a ``EstimateProbability`` bootstrap instance

.. doctest::

    >>> estimateP = bootstrap.EstimateProbability(create_null_function(),
    ...                                       create_alt_function(),
    ...                                       aln)

Specify how many random samples we want it to generate. Here we use a very small number of replicates only for the purpose of testing.

.. doctest::

    >>> estimateP.setNumReplicates(5)

Run it.

.. doctest::

    >>> estimateP.run(show_progress = False)

``show_progress`` sets whether individual optimisations are printed to screen. Get the estimated probability.

.. doctest::

    >>> p = estimateP.getEstimatedProb()

``p`` is a floating point value, as you'd expect. Grab the estimated likelihoods (null and alternate) for the observed data.

.. doctest::

    >>> print '%.2f, %.2f' % estimateP.getObservedlnL()
    -8751.94, -8750.59

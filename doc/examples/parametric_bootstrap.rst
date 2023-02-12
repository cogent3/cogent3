.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _parametric-bootstrap:

Performing a parametric bootstrap
=================================

.. sectionauthor:: Gavin Huttley

This file contains an example for estimating the probability of a Likelihood ratio statistic obtained from a relative rate test. The bootstrap classes can take advantage of parallel architectures.

From ``cogent3`` import all the components we need.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve import bootstrap
    from cogent3.evolve.models import HKY85

Define the null model that takes an alignment object and returns a likelihood function properly assembled for optimising the likelihood under the null hypothesis. The sample distribution is generated using this model.

We will use a HKY model.

.. jupyter-execute::

    def create_alt_function():
        t = load_tree("data/test.tree")
        sm = HKY85()
        return sm.make_likelihood_function(t)

Define a function that takes an alignment object and returns an appropriately assembled function for the alternative model. Since the two models are identical bar the constraint on the branch lengths, we'll use the same code to generate the basic likelihood function as for the alt model, and then apply the constraint here

.. jupyter-execute::

    def create_null_function():
        lf = create_alt_function()
        # set the local clock for humans & howler monkey
        lf.set_local_clock("Human", "HowlerMon")
        return lf

Get our observed data alignment

.. jupyter-execute::

    aln = load_aligned_seqs("data/long_testseqs.fasta")

Create a ``EstimateProbability`` bootstrap instance

.. jupyter-execute::

    estimateP = bootstrap.EstimateProbability(
        create_null_function(), create_alt_function(), aln
    )

Specify how many random samples we want it to generate. Here we use a very small number of replicates only for the purpose of testing.

.. jupyter-execute::

    estimateP.set_num_replicates(5)

Run it.

.. jupyter-execute::

    estimateP.run(show_progress=False)

Get the estimated probability.

.. jupyter-execute::

    p = estimateP.get_estimated_prob()

``p`` is a floating point value, as you'd expect. Grab the estimated likelihoods (null and alternate) for the observed data.

.. jupyter-execute::

    print("%.2f, %.2f" % estimateP.get_observed_lnL())
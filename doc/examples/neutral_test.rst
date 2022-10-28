.. jupyter-execute::
    :hide-code:

    import set_working_directory

A test of the neutral theory
============================

.. sectionauthor:: Gavin Huttley

This file contains an example for performing a likelihood ratio test of neutrality. The test compares a model where the codon model parameter omega is constrained to be the same for all edges against one where each edge has its' own omega. From ``cogent3`` import all the components we need.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model
    from scipy.stats.distributions import chi2

Get your alignment and tree.

.. jupyter-execute::

    al = load_aligned_seqs("data/long_testseqs.fasta")
    t = load_tree("data/test.tree")

We use a Goldman Yang 1994 model.

.. jupyter-execute::

    sm = get_model("MG94GTR")

Make the controller object

.. jupyter-execute::

    lf = sm.make_likelihood_function(t, digits=2, space=2)

Get the likelihood function object this object performs the actual likelihood calculation.

.. jupyter-execute::

    lf.set_alignment(al)

By default, parameters other than branch lengths are treated as global in scope, so we don't need to do anything special here. We can influence how rigorous the optimisation will be, and switch between the global and local optimisers provided in the toolkit using arguments to the optimise method. The ``global_tolerance=1.0`` argument specifies conditions for an early break from simulated annealing which will be automatically followed by the Powell local optimiser. .. note:: the 'results' are of course nonsense.

.. jupyter-execute::

    lf.optimise(global_tolerance=1.0, show_progress=False)

View the resulting maximum-likelihood parameter values

.. jupyter-execute::

    lf

We'll get the lnL and number of free parameters for later use.

.. jupyter-execute::

    null_lnL = lf.get_log_likelihood()
    null_nfp = lf.get_num_free_params()

Specify each edge has it's own omega by just modifying the existing ``lf``. This means the new function will start with the above values.

.. jupyter-execute::

    lf.set_param_rule("omega", is_independent=True)

Optimise the likelihood function, this time just using the local optimiser.

.. jupyter-execute::

    lf.optimise(local=True, show_progress=False)

View the resulting maximum-likelihood parameter values.

.. jupyter-execute::

    lf

Get out an annotated tree, it looks just like a tree, but has the maximum-likelihood parameter estimates attached to each tree edge. This object can be used for plotting, or to provide starting estimates to a related model.

.. jupyter-execute::

    at = lf.get_annotated_tree()

The lnL's from the two models are now used to calculate the likelihood ratio statistic (``LR``) it's degrees-of-freedom (``df``) and the probability (``P``) of observing the LR.

.. jupyter-execute::

    LR = 2 * (lf.get_log_likelihood() - null_lnL)
    df = lf.get_num_free_params() - null_nfp
    P = chi2.sf(LR, df)

Print this and look up a chi-sq with number of edges - 1 degrees of freedom.

.. jupyter-execute::

    print(f"Likelihood ratio statistic = {LR}")
    print(f"degrees-of-freedom = {df}")
    print(f"probability = {P}")

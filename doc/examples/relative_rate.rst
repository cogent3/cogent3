.. jupyter-execute::
    :hide-code:

    import set_working_directory

Performing a relative rate test
===============================

.. sectionauthor:: Gavin Huttley

From ``cogent3`` import all the components we need

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model
    from cogent3.maths import stats

Get your alignment and tree.

.. jupyter-execute::

    aln = load_aligned_seqs("data/long_testseqs.fasta")
    t = load_tree(filename="data/test.tree")

Create a HKY85 model.

.. jupyter-execute::

    sm = get_model("HKY85")

Make the controller object and limit the display precision (to decrease the chance that small differences in estimates cause tests of the documentation to fail).

.. jupyter-execute::

    lf = sm.make_likelihood_function(t, digits=2, space=3)

Set the local clock for humans & Howler Monkey. This method is just a special interface to the more general ``set_param_rules()`` method.

.. jupyter-execute::

    lf.set_local_clock("Human", "HowlerMon")

Get the likelihood function object this object performs the actual likelihood calculation.

.. jupyter-execute::

    lf.set_alignment(aln)

Optimise the function capturing the return optimised lnL, and parameter value vector.

.. jupyter-execute::

    lf.optimise(show_progress=False)

View the resulting maximum-likelihood parameter values.

.. jupyter-execute::

    lf.set_name("clock")
    lf

We extract the log-likelihood and number of free parameters for later use.

.. jupyter-execute::

    null_lnL = lf.get_log_likelihood()
    null_nfp = lf.get_num_free_params()

Clear the local clock constraint, freeing up the branch lengths.

.. jupyter-execute::

    lf.set_param_rule("length", is_independent=True)

Run the optimiser capturing the return optimised lnL, and parameter value vector.

.. jupyter-execute::

    lf.optimise(show_progress=False)

View the resulting maximum-likelihood parameter values.

.. jupyter-execute::

    lf.set_name("non clock")
    lf

These two lnL's are now used to calculate the likelihood ratio statistic it's degrees-of-freedom and the probability of observing the LR.

.. jupyter-execute::

    LR = 2 * (lf.get_log_likelihood() - null_lnL)
    df = lf.get_num_free_params() - null_nfp
    P = stats.chi2.sf(LR, df)

Print this and look up a :math:`\chi^2` with number of edges - 1 degrees of freedom.

.. jupyter-execute::

    print("Likelihood ratio statistic = ", LR)
    print("degrees-of-freedom = ", df)
    print("probability = ", P)

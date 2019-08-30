Performing a relative rate test
===============================

.. sectionauthor:: Gavin Huttley

From ``cogent3`` import all the components we need

.. doctest::

    >>> from cogent3 import load_aligned_seqs, load_tree
    >>> from cogent3.evolve.models import HKY85
    >>> from cogent3.maths import stats

Get your alignment and tree.

.. doctest::

    >>> aln = load_aligned_seqs("data/long_testseqs.fasta")
    >>> t = load_tree(filename="data/test.tree")

Create a HKY85 model.

.. doctest::

    >>> sm = HKY85()

Make the controller object and limit the display precision (to decrease the chance that small differences in estimates cause tests of the documentation to fail).

.. doctest::

    >>> lf = sm.make_likelihood_function(t, digits=2, space=3)

Set the local clock for humans & Howler Monkey. This method is just a special interface to the more general ``set_param_rules`` method.

.. doctest::

    >>> lf.set_local_clock("Human", "HowlerMon")

Get the likelihood function object this object performs the actual likelihood calculation.

.. doctest::

    >>> lf.set_alignment(aln)

Optimise the function capturing the return optimised lnL, and parameter value vector.

.. doctest::

    >>> lf.optimise(show_progress=False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> lf.set_name("clock")
    >>> print(lf)
    clock
    log-likelihood = -8751.9425
    number of free parameters = 7
    =====
    kappa
    -----
     4.10
    -----
    ===========================
         edge   parent   length
    ---------------------------
        Human   edge.0     0.04
    HowlerMon   edge.0     0.04
       edge.0   edge.1     0.04
        Mouse   edge.1     0.28
       edge.1     root     0.02
    NineBande     root     0.09
     DogFaced     root     0.11
    ---------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------

We extract the log-likelihood and number of free parameters for later use.

.. doctest::

    >>> null_lnL = lf.get_log_likelihood()
    >>> null_nfp = lf.get_num_free_params()

Clear the local clock constraint, freeing up the branch lengths.

.. doctest::

    >>> lf.set_param_rule('length', is_independent=True)

Run the optimiser capturing the return optimised lnL, and parameter value vector.

.. doctest::

    >>> lf.optimise(show_progress=False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> lf.set_name("non clock")
    >>> print(lf)
    non clock
    log-likelihood = -8750.5889
    number of free parameters = 8
    =====
    kappa
    -----
     4.10
    -----
    ===========================
         edge   parent   length
    ---------------------------
        Human   edge.0     0.03
    HowlerMon   edge.0     0.04
       edge.0   edge.1     0.04
        Mouse   edge.1     0.28
       edge.1     root     0.02
    NineBande     root     0.09
     DogFaced     root     0.11
    ---------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------

These two lnL's are now used to calculate the likelihood ratio statistic it's degrees-of-freedom and the probability of observing the LR.

.. doctest::

    >>> LR = 2 * (lf.get_log_likelihood() - null_lnL)
    >>> df = lf.get_num_free_params() - null_nfp
    >>> P = stats.chisqprob(LR, df)

Print this and look up a :math:`\chi^2` with number of edges - 1 degrees of freedom.

.. doctest::

    >>> print("Likelihood ratio statistic = ", LR)
    Likelihood ratio statistic =  2.7...
    >>> print("degrees-of-freedom = ", df)
    degrees-of-freedom =  1
    >>> print("probability = ", P)
    probability =  0.09...

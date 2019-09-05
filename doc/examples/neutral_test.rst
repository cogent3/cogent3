A test of the neutral theory
============================

.. sectionauthor:: Gavin Huttley

This file contains an example for performing a likelihood ratio test of neutrality. The test compares a model where the codon model parameter omega is constrained to be the same for all edges against one where each edge has its' own omega. From ``cogent3`` import all the components we need.

.. doctest::

    >>> from cogent3 import load_aligned_seqs, load_tree
    >>> from cogent3.evolve.models import MG94GTR
    >>> from cogent3.maths import stats

Get your alignment and tree.

.. doctest::

    >>> al = load_aligned_seqs("data/long_testseqs.fasta")
    >>> t = load_tree("data/test.tree")

We use a Goldman Yang 1994 model.

.. doctest::

    >>> sm = MG94GTR()

Make the controller object

.. doctest::

    >>> lf = sm.make_likelihood_function(t, digits=2, space=2)

Get the likelihood function object this object performs the actual likelihood calculation.

.. doctest::

    >>> lf.set_alignment(al)

By default, parameters other than branch lengths are treated as global in scope, so we don't need to do anything special here. We can influence how rigorous the optimisation will be, and switch between the global and local optimisers provided in the toolkit using arguments to the optimise method. The ``global_tolerance=1.0`` argument specifies conditions for an early break from simulated annealing which will be automatically followed by the Powell local optimiser. .. note:: the 'results' are of course nonsense.

.. doctest::

    >>> lf.optimise(global_tolerance=1.0, show_progress=False)

View the resulting maximum-likelihood parameter values

.. doctest::

    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8636.1801
    number of free parameters = 13
    ===================================
     A/C   A/G   A/T   C/G   C/T  omega
    -----------------------------------
    1.02  3.36  0.73  0.95  3.71   0.90
    -----------------------------------
    =========================
         edge  parent  length
    -------------------------
        Human  edge.0    0.09
    HowlerMon  edge.0    0.12
       edge.0  edge.1    0.12
        Mouse  edge.1    0.84
       edge.1    root    0.06
    NineBande    root    0.28
     DogFaced    root    0.34
    -------------------------
    ======================
       A     C     G     T
    ----------------------
    0.37  0.19  0.21  0.23
    ----------------------

We'll get the lnL and number of free parameters for later use.

.. doctest::

    >>> null_lnL = lf.get_log_likelihood()
    >>> null_nfp = lf.get_num_free_params()

Specify each edge has it's own omega by just modifying the existing ``lf``. This means the new function will start with the above values.

.. doctest::

    >>> lf.set_param_rule("omega", is_independent=True)

Optimise the likelihood function, this time just using the local optimiser.

.. doctest::

    >>> lf.optimise(local=True, show_progress=False)

View the resulting maximum-likelihood parameter values.

.. doctest::

    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8632.1355
    number of free parameters = 19
    ============================
     A/C   A/G   A/T   C/G   C/T
    ----------------------------
    1.03  3.38  0.73  0.95  3.72
    ----------------------------
    ================================
         edge  parent  length  omega
    --------------------------------
        Human  edge.0    0.09   0.59
    HowlerMon  edge.0    0.12   0.96
       edge.0  edge.1    0.11   1.13
        Mouse  edge.1    0.83   0.92
       edge.1    root    0.06   0.39
    NineBande    root    0.28   1.28
     DogFaced    root    0.34   0.84
    --------------------------------
    ======================
       A     C     G     T
    ----------------------
    0.37  0.19  0.21  0.23
    ----------------------

Get out an annotated tree, it looks just like a tree, but has the maximum-likelihood parameter estimates attached to each tree edge. This object can be used for plotting, or to provide starting estimates to a related model.

.. doctest::

    >>> at = lf.get_annotated_tree()

The lnL's from the two models are now used to calculate the likelihood ratio statistic (``LR``) it's degrees-of-freedom (``df``) and the probability (``P``) of observing the LR.

.. doctest::

    >>> LR = 2 * (lf.get_log_likelihood() - null_lnL)
    >>> df = lf.get_num_free_params() - null_nfp
    >>> P = stats.chisqprob(LR, df)

Print this and look up a chi-sq with number of edges - 1 degrees of freedom.

.. doctest::

    >>> print("Likelihood ratio statistic = ", LR)
    Likelihood ratio statistic =  8...
    >>> print("degrees-of-freedom = ", df)
    degrees-of-freedom =  6
    >>> print("probability = ", P)
    probability =  0.2...

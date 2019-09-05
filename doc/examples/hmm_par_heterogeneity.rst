.. _rate-heterogeneity-hmm:

Evaluate process heterogeneity using a Hidden Markov Model
==========================================================

.. sectionauthor:: Gavin Huttley

The existence of rate heterogeneity in the evolution of biological sequences is well known. Typically such an evolutionary property is evaluated using so-called site-heterogeneity models. These models postulate the existence of discrete classes of sites, where sites within a class evolve according to a specific rate that is distinct from the rates of the other classes. These models retain the assumption that alignment columns evolve independently. One can naturally ask the question of whether rate classes occur randomly throughout the sequence or whether they are in fact auto-correlated - meaning sites of a class tend to cluster together. Because we do not have, *a priori*, a basis for classifying the sites the models are specified such that each column can belong to any of the designated site classes and the likelihood is computed across all possible classifications. Post numerical optimisation we can calculate the posterior probability a site column belongs to a specific site class. In ``cogent3``, site classes are referred to as ``bins`` and so we refer to bin probabilities etc ...

To illustrate how to evaluate these hypotheses formally we specify 3 nested hypotheses: (i) Ho: no rate heterogeneity; (ii) Ha(1): two classes of sites - fast and slow, but independent sites; (iii) Ha(2): fast and slowly evolving sites are auto-correlated (meaning a sites class is correlated with that of its' immediate neighbours).

It is also possible to apply these models to different types of changes and we illustrate this with a single parameterisation at the end.

First import standard components necessary for all of the following calculations. As the likelihood ratio tests (LRT) involve nested hypotheses we will employ the chi-square approximation for assessing statistical significance.

.. doctest::

    >>> from cogent3.evolve.substitution_model import TimeReversibleNucleotide, predicate
    >>> from cogent3 import load_aligned_seqs, load_tree
    >>> from cogent3.maths.stats import chisqprob

Load the alignment and tree.

.. doctest::

    >>> aln = load_aligned_seqs("data/long_testseqs.fasta")
    >>> tree = load_tree("data/test.tree")

Model Ho: no rate heterogeneity
-------------------------------

We define a HKY model of nucleotide substitution, which has a transition parameter. This is defined using the ``MotifChange`` class, by specifying a transition as **not** a transversion (``~MotifChange('R','Y')``).

.. doctest::

    >>> MotifChange = predicate.MotifChange
    >>> treat_gap = dict(recode_gaps=True, model_gaps=False)
    >>> kappa = (~MotifChange('R', 'Y')).aliased('kappa')
    >>> model = TimeReversibleNucleotide(predicates=[kappa], **treat_gap)

We specify a null model with no bins, and optimise it.

.. doctest::

    >>> lf_one = model.make_likelihood_function(tree, digits=2, space=3)
    >>> lf_one.set_alignment(aln)
    >>> lf_one.optimise(show_progress=False)
    >>> lnL_one = lf_one.get_log_likelihood()
    >>> df_one = lf_one.get_num_free_params()
    >>> print(lf_one)
    Likelihood function statistics
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
       edge.1     root     0.02...

Model Ha(1): two classes of gamma distributed but independent sites
-------------------------------------------------------------------

Our next hypothesis is that there are two rate classes, or bins, with rates gamma distributed. We will restrict the bin probabilities to be equal.

.. doctest::

    >>> bin_submod = TimeReversibleNucleotide(predicates=[kappa], ordered_param='rate',
    ...                      distribution='gamma', **treat_gap)
    >>> lf_bins = bin_submod.make_likelihood_function(tree, bins=2,
    ...                             sites_independent=True, digits=2, space=3)
    >>> lf_bins.set_param_rule('bprobs', is_constant=True)
    >>> lf_bins.set_alignment(aln)
    >>> lf_bins.optimise(local=True, show_progress=False)
    >>> lnL_bins = lf_bins.get_log_likelihood()
    >>> df_bins = lf_bins.get_num_free_params()
    >>> assert df_bins == 9
    >>> print(lf_bins)
    Likelihood function statistics
    log-likelihood = -8739.0900
    number of free parameters = 9
    ==================
    kappa   rate_shape
    ------------------
     4.38         1.26
    ------------------
    ====================
     bin   bprobs   rate
    --------------------
    bin0     0.50   0.41
    bin1     0.50   1.59
    --------------------
    ===========================
         edge   parent   length
    ---------------------------
        Human   edge.0     0.03
    HowlerMon   edge.0     0.04
       edge.0   edge.1     0.04
        Mouse   edge.1     0.31...

Model Ha(2): fast and slowly evolving sites are auto-correlated
---------------------------------------------------------------

We then specify a model with switches for changing between site-classes, the HMM part. The setup is almost identical to that for above with the sole difference being setting the ``sites_independent=False``.

.. doctest::

    >>> lf_patches = bin_submod.make_likelihood_function(tree, bins=2,
    ...                         sites_independent=False, digits=2, space=3)
    >>> lf_patches.set_param_rule('bprobs', is_constant=True)
    >>> lf_patches.set_alignment(aln)
    >>> lf_patches.optimise(local=True, show_progress=False)
    >>> lnL_patches = lf_patches.get_log_likelihood()
    >>> df_patches = lf_patches.get_num_free_params()
    >>> print(lf_patches)
    Likelihood function statistics
    log-likelihood = -8728.1367
    number of free parameters = 10
    ===============================
    bin_switch   kappa   rate_shape
    -------------------------------
          0.56    4.42         1.16
    -------------------------------
    ====================
     bin   bprobs   rate
    --------------------
    bin0     0.50   0.39
    bin1     0.50   1.61
    --------------------
    ===========================
         edge   parent   length
    ---------------------------
        Human   edge.0     0.03
    HowlerMon   edge.0     0.04
       edge.0   edge.1     0.04
        Mouse   edge.1     0.31
       edge.1     root     0.02
    NineBande     root     0.10
     DogFaced     root     0.12
    ---------------------------...

We use the following short function to compute the LR test statistic.

.. doctest::

    >>> LR = lambda alt, null: 2 * (alt - null)

We conduct the test between the sequentially nested models.

.. doctest::

    >>> lr = LR(lnL_bins, lnL_one)
    >>> print(lr)
    22...
    >>> print("%.4f" % chisqprob(lr, df_patches-df_bins))
    0.0000

The stationary bin probabilities are labelled as ``bprobs`` and can be obtained as follows.

.. doctest::

    >>> bprobs = lf_patches.get_param_value('bprobs')
    >>> print("%.1f : %.1f" % tuple(bprobs))
    0.5 : 0.5

Of greater interest here (given the model was set up so the bin probabilities were equal, i.e. ``is_constant=True``) are the posterior probabilities as those allow classification of sites. The result is a ``DictArray`` class instance, which behaves like a dictionary.

.. doctest::

    >>> pp = lf_patches.get_bin_probs()

If we want to know the posterior probability the 21st position belongs to ``bin0``, we can determine it as:

.. doctest::

    >>> print(pp['bin0'][20])
    0.8...

A model with patches of ``kappa``
---------------------------------

In this example we model sequence evolution where there are 2 classes of sites distinguished by their ``kappa`` parameters. We need to know what value of ``kappa`` to specify the delineation of the bin boundaries. We can determine this from the null model (``lf_one``). For this use case, we also need to use a ``numpy.array``, so we'll import that.

.. todo::

    **FOR RELEASE** did we fix this silliness of requiring a nump.array?

.. doctest::

    >>> from numpy import array
    >>> single_kappa = lf_one.get_param_value('kappa')

We then construct the substitution model in a different way to that when evaluating generic rate heterogeneity (above).

.. doctest::

    >>> kappa_bin_submod = TimeReversibleNucleotide(predicates=[kappa], **treat_gap)
    >>> lf_kappa = kappa_bin_submod.make_likelihood_function(tree,
    ...      bins = ['slow', 'fast'], sites_independent=False, digits=1,
    ...      space=3)

To improve the likelihood fitting it is desirable to set starting values in the model that result in it's initial likelihood being that of the null model (or as close as possible). To do this, we're going to define an arbitrarily small value (``epsilon``) which we use to provide the starting value to the two bins as slightly smaller/greater than ``single_kappa`` for the slow/fast bins respectively. At the same time we set the upper/lower bin boundaries.

.. doctest::

    >>> epsilon = 1e-6
    >>> lf_kappa.set_param_rule(kappa, init=single_kappa-epsilon,
    ...                      upper=single_kappa, bin='slow')
    >>> lf_kappa.set_param_rule(kappa, init=single_kappa+epsilon,
    ...                      lower=single_kappa, bin='fast')

We then illustrate how to adjust the bin probabilities, here doing it so that one of them is nearly 1, the other nearly 0. This ensures the likelihood will be near identical to that of ``lf_one`` and as a result the optimisation step will actually improve fit over the simpler model.

.. doctest::

    >>> lf_kappa.set_param_rule('bprobs',
    ...             init=array([1.0-epsilon, 0.0+epsilon]))
    >>> lf_kappa.set_alignment(aln)
    >>> lf_kappa.optimise(local=True, show_progress=False)
    >>> print(lf_kappa)
    Likelihood function statistics
    log-likelihood = -8749.3117
    number of free parameters = 11
    ==========
    bin_switch
    ----------
           0.6
    ----------
    =====================
     bin   bprobs   kappa
    ---------------------
    slow      0.8     3.0
    fast      0.2    23.2
    ---------------------
    ===========================
         edge   parent   length
    ---------------------------
        Human   edge.0      0.0
    HowlerMon   edge.0      0.0
       edge.0   edge.1      0.0
        Mouse   edge.1      0.3
       edge.1     root      0.0
    NineBande     root      0.1
     DogFaced     root      0.1
    ---------------------------
    =====================
      A     C     G     T
    ---------------------
    0.4   0.2   0.2   0.2
    ---------------------
    >>> print(lf_kappa.get_log_likelihood())
    -8749.3...

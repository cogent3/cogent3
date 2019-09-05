.. _scope-params-on-trees:

Allowing substitution model parameters to differ between branches
=================================================================

.. sectionauthor:: Gavin Huttley

A common task concerns assessing how substitution model exchangeability parameters differ between evolutionary lineages. This is most commonly of interest for the case of testing for natural selection. Here I'll demonstrate the different ways of scoping parameters across trees for the codon model case and how these can be used for evolutionary modelling.

We start with the standard imports, plus using a canned codon substitution model and then load the sample data set.

.. doctest::

    >>> from cogent3 import load_aligned_seqs, load_tree
    >>> from cogent3.evolve.models import MG94HKY
    >>> aln = load_aligned_seqs("data/long_testseqs.fasta")
    >>> tree = load_tree("data/test.tree")

We construct the substitution model and likelihood function and set the alignment.

.. doctest::

    >>> sm = MG94HKY()
    >>> lf = sm.make_likelihood_function(tree, digits=2, space=3)
    >>> lf.set_alignment(aln)

At this point we have a likelihood function with two exchangeability parameters from the substitution model (``kappa`` the transition/transversion ratio; ``omega`` the nonsynonymous/synonymous ratio) plus branch lengths for all tree edges. To facilitate subsequent discussion I now display the tree

.. doctest::

    >>> print(tree.ascii_art())
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

In order to scope a parameter on a tree (meaning specifying a subset of edges for which the parameter is to be treated differently to the remainder of the tree) requires uniquely identifying the edges. We do this using the following arguments to the likelihood function ``set_param_rule`` method:

- ``tip_names``: the name of two tips
- ``outgroup_name``: the name of a tip that is not part of the clade of interest
- ``clade``: if ``True``, all lineages descended from the tree node identified by the ``tip_names`` and ``outgroup_name`` argument are affected by the other arguments. If ``False``, then the ``stem`` argument must apply.
- ``stem``: Whether the edge leading to the node is included.

The next concepts include exactly what can be scoped and how. In the case of testing for distinctive periods of natural selection it is common to specify distinct values for ``omega`` for an edge. I'll first illustrate some possible uses for the arguments above by setting ``omega`` to be distinctive for specific edges. I will set a value for omega so that printing the likelihood function illustrates what edges have been effected, but I won't actually do any model fitting.

Specifying a clade
------------------

I'm going to cause ``omega`` to attain a different value for all branches aside from the primate clade and stem (``HowlerMon``, ``Human``, ``edge.0``).

.. doctest::

    >>> lf.set_param_rule('omega', tip_names=['DogFaced', 'Mouse'],
    ...              outgroup_name='Human', init=2.0, clade=True)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -9489.9506
    number of free parameters = 10
    =====
    kappa
    -----
     1.00
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.03    1.00
    HowlerMon   edge.0     0.04    1.00
       edge.0   edge.1     0.04    1.00
        Mouse   edge.1     0.28    2.00
       edge.1     root     0.02    2.00
    NineBande     root     0.09    2.00
     DogFaced     root     0.11    2.00
    -----------------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------

As you can see ``omega`` for the primate edges I listed above have the default parameter value (1.0), while the others have what I've assigned. In fact, you could omit the ``clade`` argument as this is the default, but I think for readability of scripts it's best to be explicit.

Specifying a stem
-----------------

This time I'll specify the stem leading to the primates as the edge of interest.

.. note:: I need to reset the ``lf`` so all edges have the default value again. I'll show this only for this example, but rest assured I'm doing it for all others too.

.. doctest::

    >>> lf.set_param_rule('omega', init=1.0)
    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', init=2.0, stem=True, clade=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -9424.8896
    number of free parameters = 10
    =====
    kappa
    -----
     1.00
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.03    1.00
    HowlerMon   edge.0     0.04    1.00
       edge.0   edge.1     0.04    2.00
        Mouse   edge.1     0.28    1.00
       edge.1     root     0.02    1.00
    NineBande     root     0.09    1.00
     DogFaced     root     0.11    1.00
    -----------------------------------...

Specifying clade and stem
-------------------------

I'll specify that both the primates and their stem are to be considered.

.. doctest::
    :hide:

    >>> lf.set_param_rule('omega', init=1.0)

.. doctest::

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', init=2.0, stem=True, clade=True)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -9442.4271
    number of free parameters = 10
    =====
    kappa
    -----
     1.00
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.03    2.00
    HowlerMon   edge.0     0.04    2.00
       edge.0   edge.1     0.04    2.00
        Mouse   edge.1     0.28    1.00
       edge.1     root     0.02    1.00
    NineBande     root     0.09    1.00
     DogFaced     root     0.11    1.00
    -----------------------------------...

Alternate arguments for specifying edges
----------------------------------------

The likelihood function ``set_param_rule`` method also has the arguments of ``edge`` and ``edges``. These allow specific naming of the tree edge(s) to be affected by a rule. In general, however, the ``tip_names`` + ``outgroup_name`` combo is more robust.

Applications of scoped parameters
---------------------------------

The general use-cases for which a tree scope can be applied are:

1. constraining all edges identified by a rule to have a specific value which is constant and not modifiable

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True, is_constant=True)

2. all edges identified by a rule have the same but different value to the rest of the tree

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True)

3. allowing all edges identified by a rule to have different values of the parameter with the remaining tree edges having the same value

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True, is_independent=True)

4. allowing all edges to have a different value

    >>> lf.set_param_rule('omega', is_independent=True)

I'll demonstrate these cases sequentially as they involve gradually increasing the degrees of freedom in the model. First we'll constrain ``omega`` to equal 1 on the primate edges. I'll then optimise the model.

.. note:: here I'm specifying a constant value for the parameter and so I **must** use the argument ``value`` to set it. This not to be confused with the argument ``init`` that is used for providing initial (starting) values for fitting.

.. doctest::
    :hide:

    >>> lf.set_param_rule('omega', init=1.0)

.. doctest::

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True, value=1.0, is_constant=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8640.9290
    number of free parameters = 9
    =====
    kappa
    -----
     3.87
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.09    1.00
    HowlerMon   edge.0     0.12    1.00
       edge.0   edge.1     0.12    0.92
        Mouse   edge.1     0.84    0.92
       edge.1     root     0.06    0.92
    NineBande     root     0.28    0.92
     DogFaced     root     0.34    0.92
    -----------------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------
    >>> print(lf.lnL)
    -8640.9...
    >>> print(lf.nfp)
    9

I'll now free up ``omega`` on the primate clade, but making it a single value shared by all primate lineages.

.. doctest::

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True, is_constant=False)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8639.7171
    number of free parameters = 10
    =====
    kappa
    -----
     3.85
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.09    0.77
    HowlerMon   edge.0     0.12    0.77
       edge.0   edge.1     0.12    0.92
        Mouse   edge.1     0.84    0.92
       edge.1     root     0.06    0.92
    NineBande     root     0.28    0.92
     DogFaced     root     0.34    0.92
    -----------------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------
    >>> print(lf.lnL)
    -8639.7...
    >>> print(lf.nfp)
    10

Finally I'll allow all primate edges to have different values of ``omega``.

.. doctest::

    >>> lf.set_param_rule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', clade=True, is_independent=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8638.9572
    number of free parameters = 11
    =====
    kappa
    -----
     3.85
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.09    0.59
    HowlerMon   edge.0     0.12    0.95
       edge.0   edge.1     0.12    0.92
        Mouse   edge.1     0.84    0.92
       edge.1     root     0.06    0.92
    NineBande     root     0.28    0.92
     DogFaced     root     0.34    0.92
    -----------------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------
    >>> print(lf.lnL)
    -8638.9...
    >>> print(lf.nfp)
    11

We now allow ``omega`` to be different on all edges.

.. doctest::

    >>> lf.set_param_rule('omega', is_independent=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -8636.1383
    number of free parameters = 15
    =====
    kappa
    -----
     3.85
    -----
    ===================================
         edge   parent   length   omega
    -----------------------------------
        Human   edge.0     0.09    0.59
    HowlerMon   edge.0     0.12    0.95
       edge.0   edge.1     0.12    1.13
        Mouse   edge.1     0.84    0.92
       edge.1     root     0.06    0.38
    NineBande     root     0.28    1.27
     DogFaced     root     0.34    0.84
    -----------------------------------
    =========================
       A      C      G      T
    -------------------------
    0.37   0.19   0.21   0.23
    -------------------------
    >>> print(lf.lnL)
    -8636.1...
    >>> print(lf.nfp)
    15

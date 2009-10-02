Allowing substitution model parameters to differ between branches
=================================================================

.. sectionauthor:: Gavin Huttley

A common task concerns assessing how substitution model exchangeability parameters differ between evolutionary lineages. This is most commonly of interest for the case of testing for natural selection. Here I'll demonstrate the different ways of scoping parameters across trees for the codon model case and how these can be used for evolutionary modelling.

We start with the standard imports, plus using a canned codon substitution model and then load the sample data set.

.. doctest::
    
    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import MG94HKY
    >>> aln = LoadSeqs("data/long_testseqs.fasta")
    >>> tree = LoadTree("data/test.tree")

We construct the substitution model and likelihood function and set the alignment.

.. doctest::
    
    >>> sm = MG94HKY()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2, space=3)
    >>> lf.setAlignment(aln)

At this point we have a likelihood function with two exchangeability parameters from the substitution model (``kappa`` the transition/transversion ratio; ``omega`` the nonsynonymous/synonymous ratio) plus branch lengths for all tree edges. To facilitate subsequent discussion I now display the tree

.. doctest::
    
    >>> print tree.asciiArt()
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

In order to scope a parameter on a tree (meaning specifying a subset of edges for which the parameter is to be treated differently to the remainder of the tree) requires uniquely identifying the edges. We do this using the following arguments to the likelihood function ``setParamRule`` method:

- ``tip_names``: the name of two tips
- ``outgroup_name``: the name of a tip that is not part of the clade of interest
- ``is_clade``: if ``True``, all lineages descended from the tree node identified by the ``tip_names`` and ``outgroup_name`` argument are affected by the other arguments. If ``False``, then the ``is_stem`` argument must apply.
- ``is_stem``: Whether the edge leading to the node is included.

The next concepts include exactly what can be scoped and how. In the case of testing for distinctive periods of natural selection it is common to specify distinct values for ``omega`` for an edge. I'll first illustrate some possible uses for the arguments above by setting ``omega`` to be distinctive for specific edges. I will set a value for omega so that printing the likelihood function illustrates what edges have been effected, but I won't actually do any model fitting.

Specifying a clade
------------------

I'm going to cause ``omega`` to attain a different value for all branches aside from the primate clade and stem (``HowlerMon``, ``Human``, ``edge.0``).

.. doctest::
    
    >>> lf.setParamRule('omega', tip_names=['DogFaced', 'Mouse'],
    ...              outgroup_name='Human', init=2.0, is_clade=True)
    >>> print lf
    Likelihood Function Table
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
    ==============
    motif   mprobs
    --------------
        T     0.23
        C     0.19
        A     0.37
        G     0.21
    --------------

As you can see ``omega`` for the primate edges I listed above have the default parameter value (1.0), while the others have what I've assigned. In fact, you could omit the ``is_clade`` argument as this is the default, but I think for readability of scripts it's best to be explicit.

Specifying a stem
-----------------

This time I'll specify the stem leading to the primates as the edge of interest.

.. note:: I need to reset the ``lf`` so all edges have the default value again. I'll show this only for this example, but rest assured I'm doing it for all others too.

.. doctest::
    
    >>> lf.setParamRule('omega', init=1.0)
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', init=2.0, is_stem=True, is_clade=False)
    >>> print lf
    Likelihood Function Table
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
    
    >>> lf.setParamRule('omega', init=1.0)

.. doctest::
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', init=2.0, is_stem=True, is_clade=True)
    >>> print lf
    Likelihood Function Table
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

The likelihood function ``setParamRule`` method also has the arguments of ``edge`` and ``edges``. These allow specific naming of the tree edge(s) to be affected by a rule. In general, however, the ``tip_names`` + ``outgroup_name`` combo is more robust.

Applications of scoped parameters
---------------------------------

The general use-cases for which a tree scope can be applied are:

1. constraining all edges identified by a rule to have a specific value which is constant and not modifiable

    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True, is_const=True)

2. all edges identified by a rule have the same but different value to the rest of the tree
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True)

3. allowing all edges identified by a rule to have different values of the parameter with the remaining tree edges having the same value
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True, is_independent=True)

4. allowing all edges to have a different value

    >>> lf.setParamRule('omega', is_independent=True)

I'll demonstrate these cases sequentially as they involve gradually increasing the degrees of freedom in the model. First we'll constrain ``omega`` to equal 1 on the primate edges. I'll then optimise the model.

.. note:: here I'm specifying a constant value for the parameter and so I **must** use the argument ``value`` to set it. This not to be confused with the argument ``init`` that is used for providing initial (starting) values for fitting.

.. doctest::
    :hide:
    
    >>> lf.setParamRule('omega', init=1.0)

.. doctest::
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True, value=1.0, is_const=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
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
    ==============
    motif   mprobs
    --------------
        T     0.23
        C     0.19
        A     0.37
        G     0.21
    --------------
    >>> print lf.getLogLikelihood()
    -8640.9...
    >>> print lf.getNumFreeParams()
    9

I'll now free up ``omega`` on the primate clade, but making it a single value shared by all primate lineages.

.. doctest::
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True, is_const=False)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
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
    ==============
    motif   mprobs
    --------------
        T     0.23
        C     0.19
        A     0.37
        G     0.21
    --------------
    >>> print lf.getLogLikelihood()
    -8639.7...
    >>> print lf.getNumFreeParams()
    10

Finally I'll allow all primate edges to have different values of ``omega``.

.. doctest::
    
    >>> lf.setParamRule('omega', tip_names=['Human', 'HowlerMon'],
    ...      outgroup_name='Mouse', is_clade=True, is_independent=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
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
    ==============
    motif   mprobs
    --------------
        T     0.23
        C     0.19
        A     0.37
        G     0.21
    --------------
    >>> print lf.getLogLikelihood()
    -8638.9...
    >>> print lf.getNumFreeParams()
    11

We now allow ``omega`` to be different on all edges.

.. doctest::
    
    >>> lf.setParamRule('omega', is_independent=True)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
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
    ==============
    motif   mprobs
    --------------
        T     0.23
        C     0.19
        A     0.37
        G     0.21
    --------------
    >>> print lf.getLogLikelihood()
    -8636.1...
    >>> print lf.getNumFreeParams()
    15

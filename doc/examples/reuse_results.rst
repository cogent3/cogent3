Reusing results to speed up optimisation
========================================

An example of how to use the maximum-likelihood parameter estimates from one model as starting values for another model. In this file we do something silly, by saving a result and then reloading it. This is silly because the analyses are run consecutively. A better approach when running consecutively is to simply use the annotated tree directly.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import Y98

We'll create a simple model, optimise it and save it for later reuse

.. doctest::

    >>> al = LoadSeqs("data/test.paml")
    >>> t = LoadTree("data/test.tree")
    >>> sm = Y98()
    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setTablesFormat(digits=2,space=2)
    >>> lf.setAlignment(al)
    >>> lf.optimise(local=True, show_progress=False)
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

The essential object for reuse is an annotated tree these capture the parameter estimates from the above optimisation we can either use this directly in the same run, or we can save the tree to file in ``xml`` format and reload the tree at a later time for use. In this example I'll illustrate the latter scenario.

.. doctest::

    >>> at=lf.getAnnotatedTree()
    >>> at.writeToFile('tree.xml')

We load the tree as per usual

.. doctest::

    >>> nt = LoadTree('tree.xml')

Now create a more parameter rich model, in this case by allowing the ``Human`` edge to have a different value of ``omega``. By providing the annotated tree, the parameter estimates from the above run will be used as starting values for the new model.

.. doctest::

    >>> new_lf = sm.makeLikelihoodFunction(nt)
    >>> new_lf.setTablesFormat(digits=2,space=2)
    >>> new_lf.setParamRule('omega', edge='Human',
    ... is_independent=True)
    >>> new_lf.setAlignment(al)
    >>> new_lf.optimise(local=True, show_progress=False)
    >>> print new_lf
    Likelihood Function Table
    =====
    kappa
    -----
     9.07
    -----
    ====================================
         edge  parent  length      omega
    ------------------------------------
        Human  edge.0    0.10  999999.97
    HowlerMon  edge.0    0.05       1.57
       edge.0  edge.1    0.06       1.57
        Mouse  edge.1    0.90       1.57
       edge.1    root    0.00       1.57
    NineBande    root    0.11       1.57
     DogFaced    root    0.18       1.57...

.. note:: A parameter rich model applied to a small data set is unreliable.

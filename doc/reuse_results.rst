Reusing results to speed up optimisation
========================================

An example of how to use the maximum-likelihood parameter estimates from one model as starting values for another model. In this file we do something silly, by saving a result and then reloading it. This is silly because the analyses are run consecutively. A better approach when running consecutively is to simply use the annotated tree directly.

.. pycode::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import Y98

We'll create a simple model, optimise it and save it for later reuse

.. pycode::

    >>> al = LoadSeqs("data/test.paml")
    >>> t = LoadTree("data/test.tree")
    >>> sm = Y98()
    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setAlignment(al)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
    ================
     kappa     omega
    ----------------
    9.2759    1.8713
    ----------------
    =============================
         edge    parent    length
    -----------------------------
        Human    edge.0    0.0968
    HowlerMon    edge.0    0.0540
       edge.0    edge.1    0.0654
        Mouse    edge.1    0.9116
       edge.1      root    0.0000
    NineBande      root    0.1073
     DogFaced      root    0.1801...

The essential object for reuse is an annotated tree these capture the parameter estimates from the above optimisation we can either use this directly in the same run, or we can save the tree to file in ``xml`` format and reload the tree at a later time for use. In this example I'll illustrate the latter scenario.

.. pycode::

    >>> at=lf.getAnnotatedTree()
    >>> at.writeToFile('tree.xml')

We load the tree as per usual

.. pycode::

    >>> nt = LoadTree('tree.xml')

Now create a more parameter rich model, in this case by allowing the ``Human`` edge to have a different value of ``omega``. By providing the annotated tree, the parameter estimates from the above run will be used as starting values for the new model.

.. pycode::

    >>> new_lf = sm.makeLikelihoodFunction(nt)
    >>> new_lf.setParamRule('omega', edge='Human',
    ... is_independent=True)
    >>> new_lf.setAlignment(al)
    >>> new_lf.optimise(local=True, show_progress=False)
    >>> print new_lf
    Likelihood Function Table
    ======
     kappa
    ------
    9.0706
    ------
    ============================================
         edge    parent    length          omega
    --------------------------------------------
        Human    edge.0    0.1001    999999.9965
    HowlerMon    edge.0    0.0510         1.5666
       edge.0    edge.1    0.0649         1.5666
        Mouse    edge.1    0.8984         1.5666
       edge.1      root    0.0000         1.5666
    NineBande      root    0.1064         1.5666
     DogFaced      root    0.1793         1.5666...

:Note: A parameter rich model applied to a small data set is unreliable.

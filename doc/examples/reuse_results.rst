Reusing results to speed up optimisation
========================================

.. sectionauthor:: Gavin Huttley

An example of how to use the maximum-likelihood parameter estimates from one model as starting values for another model. In this file we do something silly, by saving a result and then reloading it. This is silly because the analyses are run consecutively. A better approach when running consecutively is to simply use the annotated tree directly.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import MG94HKY

We'll create a simple model, optimise it and save it for later reuse

.. doctest::

    >>> aln = LoadSeqs("data/long_testseqs.fasta")
    >>> t = LoadTree("data/test.tree")
    >>> sm = MG94HKY()
    >>> lf = sm.makeLikelihoodFunction(t, digits=2, space=2)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
    ============
    kappa  omega
    ------------
     3.85   0.90
    ------------
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
    =============
    motif  mprobs
    -------------
        T    0.23
        C    0.19
        A    0.37
        G    0.21
    -------------

The essential object for reuse is an annotated tree these capture the parameter estimates from the above optimisation we can either use this directly in the same run, or we can save the tree to file in ``xml`` format and reload the tree at a later time for use. In this example I'll illustrate the latter scenario.

.. doctest::

    >>> at=lf.getAnnotatedTree()
    >>> at.writeToFile('tree.xml')

We load the tree as per usual

.. doctest::

    >>> nt = LoadTree('tree.xml')

Now create a more parameter rich model, in this case by allowing the ``Human`` edge to have a different value of ``omega``. By providing the annotated tree, the parameter estimates from the above run will be used as starting values for the new model.

.. doctest::

    >>> new_lf = sm.makeLikelihoodFunction(nt, digits=2, space=2)
    >>> new_lf.setParamRule('omega', edge='Human',
    ...                     is_independent=True)
    >>> new_lf.setAlignment(aln)
    >>> new_lf.optimise(local=True, show_progress=False)
    >>> print new_lf
    Likelihood Function Table
    =====
    kappa
    -----
     3.85
    -----
    ================================
         edge  parent  length  omega
    --------------------------------
        Human  edge.0    0.09   0.59
    HowlerMon  edge.0    0.12   0.92
       edge.0  edge.1    0.12   0.92
        Mouse  edge.1    0.84   0.92
       edge.1    root    0.06   0.92
    NineBande    root    0.28   0.92
     DogFaced    root    0.34   0.92
    --------------------------------
    =============
    motif  mprobs
    -------------
        T    0.23
        C    0.19
        A    0.37
        G    0.21
    -------------

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('tree.xml')

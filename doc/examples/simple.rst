The simplest script
===================

.. sectionauthor:: Gavin Huttley

This is just about the simplest possible Cogent script for evolutionary modelling. We use a canned nucleotide substitution model (the ``HKY85`` model) on just three primate species. As there is only unrooted tree possible, the sequence names are all that's required to make the tree.

.. doctest::

    >>> from cogent.evolve.models import HKY85
    >>> from cogent import LoadSeqs, LoadTree
    >>> model = HKY85()
    >>> aln = LoadSeqs("data/primate_cdx2_promoter.fasta")
    >>> tree = LoadTree(tip_names=aln.Names)
    >>> lf = model.makeLikelihoodFunction(tree)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(show_progress = False)
    >>> print lf
    Likelihood Function Table
    ======
     kappa
    ------
    5.9589
    ------
    ===========================
       edge    parent    length
    ---------------------------
      human      root    0.0040
    macaque      root    0.0384
      chimp      root    0.0061
    ---------------------------
    ===============
    motif    mprobs
    ---------------
        T    0.2552
        C    0.2581
        A    0.2439
        G    0.2428
    ---------------

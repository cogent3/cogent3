The simplest script
===================

.. sectionauthor:: Gavin Huttley

This is just about the simplest possible ``cogent3`` script for evolutionary modelling. We use a canned nucleotide substitution model (the ``HKY85`` model) on just three primate species. As there is only one unrooted tree possible, the sequence names are all that's required to make the tree.

.. doctest::

    >>> from cogent3.evolve.models import HKY85
    >>> from cogent3 import LoadSeqs, LoadTree
    >>> model = HKY85()
    >>> aln = LoadSeqs("data/primate_cdx2_promoter.fasta")
    >>> tree = LoadTree(tip_names=aln.names)
    >>> lf = model.make_likelihood_function(tree)
    >>> lf.set_alignment(aln)
    >>> lf.optimise(show_progress=False)
    >>> print(lf)
    Likelihood function statistics
    log-likelihood = -2494.9537
    number of free parameters = 4
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

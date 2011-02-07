********************
Building phylogenies
********************

.. Anuj Pahwa, Gavin Huttley

Built-in Phylogenetic reconstruction
====================================

By distance method
------------------

Given an alignment, a phylogenetic tree can be generated based on the pair-wise distance matrix computed from the alignment.

Fast pairwise distance estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a limited number of evolutionary models a fast implementation is available. Here we use the Tamura and Nei 1993 model.

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.evolve.pairwise_distance import TN93Pair
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> dist_calc = TN93Pair(DNA, alignment=aln)
    >>> dist_calc.run()

We can obtain the distances as a ``dict`` for direct usage in phylogenetic reconstruction

.. doctest::
    
    >>> dists = dist_calc.getPairwiseDistances()

or as a table for display / saving

.. doctest::
    
    >>> print dist_calc.Dists[:4,:4] # truncated to fit screens
    Pairwise Distances
    ============================================
    Seq1 \ Seq2    Galago    HowlerMon    Rhesus
    --------------------------------------------
         Galago         *       0.2157    0.1962
      HowlerMon    0.2157            *    0.0736
         Rhesus    0.1962       0.0736         *
      Orangutan    0.1944       0.0719    0.0411
    --------------------------------------------

Other statistics are also available, such the as the standard errors of the estimates.

.. doctest::
    
    >>> print dist_calc.StdErr[:4,:4] # truncated to fit screens
    Standard Error of Pairwise Distances
    ============================================
    Seq1 \ Seq2    Galago    HowlerMon    Rhesus
    --------------------------------------------
         Galago         *       0.0103    0.0096
      HowlerMon    0.0103            *    0.0054
         Rhesus    0.0096       0.0054         *
      Orangutan    0.0095       0.0053    0.0039
    --------------------------------------------


More general estimation of pairwise distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The standard cogent likelihood function can also be used to estimate distances. Because these require numerical optimisation they can be significantly slower than the fast estimation approach above.

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.phylo import distance
    >>> from cogent.evolve.models import F81
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> d = distance.EstimateDistances(aln, submodel=F81())
    >>> d.run()

The example above will use the F81 nucleotide substitution model and run the ``distance.EstimateDistances()`` method with the default options for the optimiser. To configure the optimiser a dictionary of optimisation options can be passed onto the ``run`` command. The example below configures the ``Powell`` optimiser to run a maximum of 10000 evaluations, with a maximum of 5 restarts (a total of 5 x 10000 = 50000 evaluations).

.. doctest::
    
    >>> dist_opt_args = dict(max_restarts=5, max_evaluations=10000)
    >>> d.run(dist_opt_args=dist_opt_args)
    >>> print d
    ============================================================================================
    Seq1 \ Seq2    Galago    HowlerMon    Rhesus    Orangutan    Gorilla     Human    Chimpanzee
    --------------------------------------------------------------------------------------------
         Galago         *       0.2112    0.1930       0.1915     0.1891    0.1934        0.1892
      HowlerMon    0.2112            *    0.0729       0.0713     0.0693    0.0729        0.0697
         Rhesus    0.1930       0.0729         *       0.0410     0.0391    0.0421        0.0395
      Orangutan    0.1915       0.0713    0.0410            *     0.0136    0.0173        0.0140
        Gorilla    0.1891       0.0693    0.0391       0.0136          *    0.0086        0.0054
          Human    0.1934       0.0729    0.0421       0.0173     0.0086         *        0.0089
     Chimpanzee    0.1892       0.0697    0.0395       0.0140     0.0054    0.0089             *
    --------------------------------------------------------------------------------------------

Building A Phylogenetic Tree From Pairwise Distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phylogenetic Trees can be built by using the neighbour joining algorithm by providing a dictionary of pairwise distances. This dictionary can be obtained either from the output of ``distance.EstimateDistances()``

.. doctest::
    
    >>> from cogent.phylo import nj
    >>> njtree = nj.nj(d.getPairwiseDistances())
    >>> njtree = njtree.balanced()
    >>> print njtree.asciiArt()
                        /-Rhesus
              /edge.1--|
             |         |          /-HowlerMon
             |          \edge.0--|
             |                    \-Galago
    -root----|
             |--Orangutan
             |
             |          /-Human
              \edge.2--|
                       |          /-Gorilla
                        \edge.3--|
                                  \-Chimpanzee

Or created manually as shown below.

.. doctest::
    
    >>> dists = {('a', 'b'): 2.7, ('c', 'b'): 2.33, ('c', 'a'): 0.73}
    >>> njtree2 = nj.nj(dists)
    >>> print njtree2.asciiArt()
              /-a
             |
    -root----|--b
             |
              \-c

By least-squares
----------------

We illustrate the phylogeny reconstruction by least-squares using the F81 substitution model. We use the advanced-stepwise addition algorithm to search tree space. Here ``a`` is the number of taxa to exhaustively evaluate all possible phylogenies for. Successive taxa will are added to the top ``k`` trees (measured by the least-squares metric) and ``k`` trees are kept at each iteration.

.. doctest::
    
    >>> import cPickle
    >>> from cogent.phylo.least_squares import WLS
    >>> dists = cPickle.load(open('data/dists_for_phylo.pickle'))
    >>> ls = WLS(dists)
    >>> stat, tree = ls.trex(a = 5, k = 5, show_progress = False)

Other optional arguments that can be passed to the ``trex`` method are: ``return_all``, whether the ``k`` best trees at the final step are returned as a ``ScoredTreeCollection`` object; ``order``, a series of tip names whose order defines the sequence in which tips will be added during tree building (this allows the user to randomise the input order).

By ML
-----

We illustrate the phylogeny reconstruction using maximum-likelihood using the F81 substitution model. We use the advanced-stepwise addition algorithm to search tree space, setting 

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.phylo.maximum_likelihood import ML
    >>> from cogent.evolve.models import F81
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> ml = ML(F81(), aln)

The ``ML`` object also has the ``trex`` method and this can be used in the same way as for above, i.e. ``ml.trex()``. We don't do that here because this is a very slow method for phylogenetic reconstruction.

Building phylogenies with 3rd-party apps such as FastTree or RAxML
==================================================================

A thorough description is :ref:`appcontroller-phylogeny`.

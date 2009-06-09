.. _calculating-pairwise-distances:

Calculate pairwise distances between sequences
==============================================

.. sectionauthor:: Gavin Huttley

An example of how to calculate the pairwise distances for a set of sequences.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance

Import a substitution model (or create your own)

.. doctest::

    >>> from cogent.evolve.models import HKY85

Load my alignment

.. doctest::

    >>> al = LoadSeqs("data/long_testseqs.fasta")

Create a pairwise distances object with your alignment and substitution model

.. doctest::

    >>> d = distance.EstimateDistances(al, submodel= HKY85())

Printing ``d`` before execution shows its status.

.. doctest::

    >>> print d
    =========================================================================
    Seq1 \ Seq2       Human    HowlerMon       Mouse    NineBande    DogFaced
    -------------------------------------------------------------------------
          Human           *     Not Done    Not Done     Not Done    Not Done
      HowlerMon    Not Done            *    Not Done     Not Done    Not Done
          Mouse    Not Done     Not Done           *     Not Done    Not Done
      NineBande    Not Done     Not Done    Not Done            *    Not Done
       DogFaced    Not Done     Not Done    Not Done     Not Done           *
    -------------------------------------------------------------------------

Which in this case is to simply indicate nothing has been done.

.. doctest::

    >>> d.run(show_progress=False)
    >>> print d
    =====================================================================
    Seq1 \ Seq2     Human    HowlerMon     Mouse    NineBande    DogFaced
    ---------------------------------------------------------------------
          Human         *       0.0730    0.3363       0.1804      0.1972
      HowlerMon    0.0730            *    0.3487       0.1865      0.2078
          Mouse    0.3363       0.3487         *       0.3813      0.4022
      NineBande    0.1804       0.1865    0.3813            *      0.2019
       DogFaced    0.1972       0.2078    0.4022       0.2019           *
    ---------------------------------------------------------------------

Note that pairwise distances can be distributed for computation across multiple CPU's. In this case, when statistics (like distances) are requested only the master CPU returns data.

We'll write a phylip formatted distance matrix.

.. doctest::

    >>> d.writeToFile('dists_for_phylo.phylip', format="phylip")

We'll also save the distances to file in Python's pickle format.

.. doctest::

    >>> import cPickle
    >>> f = open('dists_for_phylo.pickle', "w")
    >>> cPickle.dump(d.getPairwiseDistances(), f)
    >>> f.close()

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> for file_name in 'dists_for_phylo.phylip', 'dists_for_phylo.pickle':
    ...     os.remove(file_name)

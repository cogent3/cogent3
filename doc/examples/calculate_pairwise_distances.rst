Calculate pairwise distances between sequences
==============================================

An example of how to calculate the pairwise distances for a set of sequences.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance

Import a substitution model (or create your own)

.. doctest::

    >>> from cogent.evolve.models import HKY85

Load my alignment

.. doctest::

    >>> al = LoadSeqs("data/test.paml")

Create a pairwise distances object with your alignment and substitution model

.. doctest::

    >>> d = distance.EstimateDistances(al, submodel= HKY85())

Printing ``d`` before execution shows its status.

.. doctest::

    >>> print d
    =========================================================================
    Seq1 \ Seq2    NineBande       Mouse       Human    HowlerMon    DogFaced
    -------------------------------------------------------------------------
      NineBande            *    Not Done    Not Done     Not Done    Not Done
          Mouse     Not Done           *    Not Done     Not Done    Not Done
          Human     Not Done    Not Done           *     Not Done    Not Done
      HowlerMon     Not Done    Not Done    Not Done            *    Not Done
       DogFaced     Not Done    Not Done    Not Done     Not Done           *
    -------------------------------------------------------------------------

Which in this case is to simply indicate nothing has been done.

.. doctest::

    >>> d.run(show_progress=False)
    >>> print d
    =====================================================================
    Seq1 \ Seq2    NineBande     Mouse     Human    HowlerMon    DogFaced
    ---------------------------------------------------------------------
      NineBande            *    0.2196    0.0890       0.0700      0.0891
          Mouse       0.2196         *    0.2737       0.2736      0.2467
          Human       0.0890    0.2737         *       0.0530      0.1092
      HowlerMon       0.0700    0.2736    0.0530            *      0.0894
       DogFaced       0.0891    0.2467    0.1092       0.0894           *
    ---------------------------------------------------------------------

Note that pairwise distances can be distributed for computation across multiple CPU's. In this case, when statistics (like distances) are requested only the master CPU returns data.

We'll write a phylip formatted distance matrix.

.. doctest::

    >>> d.writeToFile('junk.phylip', format="phylip")

We'll also save the distances to file in Python's pickle format.

.. doctest::

    >>> import cPickle
    >>> f = open('dists_for_phylo.pickle', "w")
    >>> cPickle.dump(d.getPairwiseDistances(), f)
    >>> f.close()

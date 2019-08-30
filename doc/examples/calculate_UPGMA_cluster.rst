Make a UPGMA cluster
====================

.. sectionauthor:: Catherine Lozupone

An example of how to calculate the pairwise distances for a set of sequences.

**NOTE:** UPGMA should not be used for phylogenetic reconstruction.

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> from cogent3.evolve import distance
    >>> from cogent3.cluster.UPGMA import upgma

Import a substitution model (or create your own)

.. doctest::

    >>> from cogent3.evolve.models import HKY85

Load the alignment.

.. doctest::

    >>> al = load_aligned_seqs("data/test.paml")

Create a pairwise distances object calculator for the alignment, providing a substitution model instance.

.. doctest::

    >>> d = distance.EstimateDistances(al, submodel=HKY85())
    >>> d.run(show_progress=False)

Now use this matrix to build a UPGMA cluster.

.. doctest::

    >>> mycluster = upgma(d.get_pairwise_distances())
    >>> print(mycluster.ascii_art())  # doctest: +SKIP
                                  /-NineBande
                        /edge.1--|
                       |         |          /-HowlerMon
              /edge.0--|          \edge.2--|
             |         |                    \-Human
    -root----|         |
             |          \-DogFaced
             |
              \-Mouse

We demonstrate saving this UPGMA cluster to a file.

.. doctest::

    >>> mycluster.write('test_upgma.tree')

..
    We don't actually want to keep that file now, so I'm importing the ``os`` module to delete it.

.. doctest::
    :hide:

    >>> import os
    >>> os.remove('test_upgma.tree')

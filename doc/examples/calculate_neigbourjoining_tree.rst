Make a neighbor joining tree
============================

.. sectionauthor:: Gavin Huttley

An example of how to calculate the pairwise distances for a set of sequences.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance, nj

Import a substitution model (or create your own)

.. doctest::

    >>> from cogent.evolve.models import HKY85

Load the alignment.

.. doctest::

    >>> al = LoadSeqs("data/long_testseqs.fasta")

Create a pairwise distances object calculator for the alignment, providing a substitution model instance.

.. doctest::

    >>> d = distance.EstimateDistances(al, submodel= HKY85())
    >>> d.run(show_progress=False)

Now use this matrix to build a neighbour joining tree.

.. doctest::

    >>> mytree = nj.nj(d.getPairwiseDistances())

We can visualise this tree by ``print mytree.asciiArt()``, which generates the equivalent of:

.. code-block:: python
    
                        /-Human
              /edge.0--|
             |          \-HowlerMon
             |
    -root----|          /-NineBande
             |-edge.1--|
             |          \-DogFaced
             |
              \-Mouse

We can save this tree to file.

.. doctest::

    >>> mytree.writeToFile('test_nj.tree')

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('test_nj.tree')

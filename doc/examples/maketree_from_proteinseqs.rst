Making a phylogenetic tree from a protein sequence alignment
============================================================

.. sectionauthor:: Gavin Huttley

In this example we pull together the distance calculation and tree building with the additional twist of using an empirical protein substitution matrix. We will therefore be computing the tree from a protein sequence alignment. We will first do the standard cogent import for ``LoadSeqs``.

.. doctest::

    >>> from cogent import LoadSeqs, PROTEIN

We will use an empirical protein substitution matrix.

.. doctest::

    >>> from cogent.evolve.models import JTT92

The next components we need are for computing the matrix of pairwise sequence distances and then for estimating a neighbour joining tree from those distances.

.. doctest::

    >>> from cogent.phylo import nj, distance

Now load our sequence alignment, explicitly setting the alphabet to be protein.

.. doctest::

    >>> aln = LoadSeqs('data/abglobin_aa.phylip', interleaved=True,
    ...                 moltype=PROTEIN)

Create an Empirical Protein Matrix Substitution model object. This will take the unscaled empirical matrix and use it and the motif frequencies to create a scaled Q matrix.

.. doctest::

    >>> sm = JTT92()

We now use this and the alignment to construct a distance calculator.

.. doctest::

    >>> d = distance.EstimateDistances(aln, submodel = sm)
    >>> d.run(show_progress=False)

The resulting distances are passed to the nj function.

.. doctest::

    >>> mytree = nj.nj(d.getPairwiseDistances())

The shape of the resulting tree can be readily view by printing ``mytree.asciiArt()``. The result will be equivalent to.

.. code-block:: python

              /-human
             |
             |          /-rabbit
    -root----|-edge.1--|
             |          \-rat
             |
             |          /-goat-cow
              \edge.0--|
                        \-marsupial

This tree can be saved to file, the ``with_distances`` argument specifies that branch lengths are to be included in the newick formatted output.

.. doctest::

    >>> mytree.writeToFile('test_nj.tree', with_distances=True)

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('test_nj.tree')

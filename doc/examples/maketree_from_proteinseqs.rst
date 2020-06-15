Making a phylogenetic tree from a protein sequence alignment
============================================================

.. sectionauthor:: Gavin Huttley

In this example we pull together the distance calculation and tree building with the additional twist of using an empirical protein substitution matrix. We will therefore be computing the tree from a protein sequence alignment. We will first do the standard ``cogent3`` import for ``load_aligned_seqs``.

.. jupyter-execute::
    :linenos:

    from cogent3 import load_aligned_seqs

We will use an empirical protein substitution matrix.

.. jupyter-execute::
    :linenos:

    from cogent3.evolve.models import JTT92

The next components we need are for computing the matrix of pairwise sequence distances and then for estimating a neighbour joining tree from those distances.

.. jupyter-execute::
    :linenos:

    from cogent3.phylo import nj
    from cogent3.evolve import distance

Now load our sequence alignment, explicitly setting the alphabet to be protein.

.. jupyter-execute::
    :linenos:

    aln = load_aligned_seqs("data/abglobin_aa.phylip", moltype="protein")

Create an Empirical Protein Matrix Substitution model object. This will take the unscaled empirical matrix and use it and the motif frequencies to create a scaled Q matrix.

.. jupyter-execute::
    :linenos:

    sm = JTT92()

We now use this and the alignment to construct a distance calculator.

.. jupyter-execute::
    :linenos:

    d = distance.EstimateDistances(aln, submodel=sm)
    d.run(show_progress=False)

The resulting distances are passed to the nj function.

.. jupyter-execute::
    :linenos:

    mytree = nj.nj(d.get_pairwise_distances())

The shape of the resulting tree can be readily view by printing ``mytree.ascii_art()``. The result will be equivalent to.

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

.. jupyter-execute::
    :linenos:

    mytree.write("test_nj.tree", with_distances=True)

.. clean up

.. jupyter-execute::
    :hide-code:

    import os

    os.remove("test_nj.tree")

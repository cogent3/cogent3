.. jupyter-execute::
    :hide-code:

    import set_working_directory

********************
Building phylogenies
********************

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

Building A Phylogenetic Tree From Pairwise Distances
====================================================

Directly via ``alignment.quick_tree()``
=======================================

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    tree = aln.quick_tree(calc="TN93")
    tree = tree.balanced()  # purely for display
    print(tree.ascii_art())

Directly via ``alignment.quick_tree()`` with a third-party hook
===============================================================

You can use the IQ-TREE rapid-NJ algorithm for ``quick_tree()`` by installing `piqtree <https://piqtree.readthedocs.io/en/latest/>`_ and setting ``use_hook="piqtree"``.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    tree = aln.quick_tree(calc="TN93", use_hook="piqtree")
    # dropping branch lengths to simplify display
    dnd = tree.get_figure(contemporaneous=True, width=600, height=600)
    dnd.show()

Using the ``DistanceMatrix`` object
-----------------------------------

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dists = aln.distance_matrix(calc="TN93")
    tree = dists.quick_tree()
    tree = tree.balanced()  # purely for display
    print(tree.ascii_art())

Explicitly via ``DistanceMatrix`` and ``cogent3.phylo.nj.nj()```
----------------------------------------------------------------

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.phylo import nj

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dists = aln.distance_matrix(calc="TN93")
    tree = nj.nj(dists, show_progress=False)
    tree = tree.balanced()  # purely for display
    print(tree.ascii_art())

Directly from a pairwise distance ``dict``
------------------------------------------

.. jupyter-execute::

    from cogent3.phylo import nj

    dists = {("a", "b"): 2.7, ("c", "b"): 2.33, ("c", "a"): 0.73}
    tree = nj.nj(dists, show_progress=False)
    print(tree.ascii_art())

By Least-squares
================

We illustrate the phylogeny reconstruction by least-squares using the F81 substitution model. We use the advanced-stepwise addition algorithm to search tree space. Here ``a`` is the number of taxa to exhaustively evaluate all possible phylogenies for. Successive taxa are added to the top ``k`` trees (measured by the least-squares metric) and ``k`` trees are kept at each iteration.

.. jupyter-execute::

    from cogent3.phylo.least_squares import WLS
    from cogent3.util.deserialise import deserialise_object

    dists = deserialise_object("data/dists_for_phylo.json")
    ls = WLS(dists)
    stat, tree = ls.trex(a=5, k=5, show_progress=False)

Other optional arguments that can be passed to the ``trex`` method are: ``return_all``, whether the ``k`` best trees at the final step are returned as a ``ScoredTreeCollection`` object; ``order``, a series of tip names whose order defines the sequence in which tips will be added during tree building (this allows the user to randomise the input order).

By ML
=====

We illustrate the phylogeny reconstruction using maximum-likelihood using the F81 substitution model. We use the advanced-stepwise addition algorithm to search tree space.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.evolve.models import F81
    from cogent3.phylo.maximum_likelihood import ML

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    ml = ML(F81(), aln)

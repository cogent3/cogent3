.. jupyter-execute::
    :hide-code:

    import set_working_directory

Make a neighbor joining tree
============================

.. sectionauthor:: Gavin Huttley

An example of how to calculate the pairwise distances for a set of sequences.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.evolve import distance
    from cogent3.phylo import nj

Import a substitution model (or create your own)

.. jupyter-execute::

    from cogent3.evolve.models import get_model

Load the alignment.

.. jupyter-execute::

    al = load_aligned_seqs("data/long_testseqs.fasta")

Create a pairwise distances object calculator for the alignment, providing a substitution model instance.

.. jupyter-execute::

    d = distance.EstimateDistances(al, submodel=get_model("HKY85"))
    d.run(show_progress=False)

Now use this matrix to build a neighbour joining tree.

.. jupyter-execute::

    mytree = nj.nj(d.get_pairwise_distances(), show_progress=False)
    print(mytree.ascii_art())

We can save this tree to file.

.. jupyter-execute::

    mytree.write("test_nj.tree")

.. clean up

.. jupyter-execute::
    :hide-code:

    import os

    os.remove("test_nj.tree")

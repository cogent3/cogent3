.. jupyter-execute::
    :hide-code:

    import set_working_directory

Make a UPGMA cluster
====================

.. sectionauthor:: Catherine Lozupone

An example of how to calculate the pairwise distances for a set of sequences.

.. note:: UPGMA should not be used for phylogenetic reconstruction.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.cluster.UPGMA import upgma
    from cogent3.evolve import distance

Import a substitution model (or create your own)

.. jupyter-execute::

    from cogent3.evolve.models import HKY85

Load the alignment.

.. jupyter-execute::

    al = load_aligned_seqs("data/test.paml")

Create a pairwise distances object calculator for the alignment, providing a substitution model instance.

.. jupyter-execute::

    d = distance.EstimateDistances(al, submodel=HKY85())
    d.run(show_progress=False)

Now use this matrix to build a UPGMA cluster.

.. jupyter-execute::

    mycluster = upgma(d.get_pairwise_distances())
    print(mycluster.ascii_art())

We demonstrate saving this UPGMA cluster to a file.

.. jupyter-execute::

    mycluster.write("test_upgma.tree")

..  We don't actually want to keep that file now, so I'm importing the ``os`` module to delete it.

.. jupyter-execute::
    :hide-code:

    import os

    os.remove("test_upgma.tree")

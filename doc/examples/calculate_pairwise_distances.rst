.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _calculating-pairwise-distances:

Calculate pairwise distances between sequences
==============================================

.. sectionauthor:: Gavin Huttley

An example of how to calculate the pairwise distances for a set of sequences.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.evolve import distance

Import a substitution model (or create your own)

.. jupyter-execute::

    from cogent3.evolve.models import HKY85

Load my alignment

.. jupyter-execute::

    al = load_aligned_seqs("data/long_testseqs.fasta")

Create a pairwise distances object with your alignment and substitution model and run it.

.. jupyter-execute::

    d = distance.EstimateDistances(al, submodel=HKY85())
    d.run(show_progress=False)
    d.get_pairwise_distances()

Note that pairwise distances can be distributed for computation across multiple CPU's. In this case, when statistics (like distances) are requested only the master CPU returns data.

We'll write a phylip formatted distance matrix.

.. jupyter-execute::

    d.write("dists_for_phylo.phylip", format="phylip")

.. todo:: write out in json format

We'll also save the distances to file in Python's pickle format.

.. jupyter-execute::

    import pickle

    with open("dists_for_phylo.pickle", "wb") as f:
        pickle.dump(d.get_pairwise_distances(), f)

.. clean up

.. jupyter-execute::
    :hide-code:

    import os

    for file_name in "dists_for_phylo.phylip", "dists_for_phylo.pickle":
        os.remove(file_name)

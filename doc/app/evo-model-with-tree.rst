.. jupyter-execute::
    :hide-code:

    import set_working_directory

Apply a non-stationary nucleotide model to an alignment with a tree
===================================================================

We analyse an alignment with sequences from 6 primates.

.. jupyter-execute::

    from cogent3.app import io

    reader = io.load_aligned(format="fasta", moltype="dna")
    aln = reader("data/primate_brca1.fasta")
    aln.names

Specify the tree via a tree instance
------------------------------------

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.app import evo

    tree = load_tree("data/primate_brca1.tree")
    gn = evo.model("GN", tree=tree)
    gn

Specify the tree via a path.
----------------------------

.. jupyter-execute::

    gn = evo.model("GN", tree="data/primate_brca1.tree")
    gn

Apply the model to an alignment
-------------------------------

.. jupyter-execute::

    fitted = gn(aln)
    fitted

In the above, no value is shown for ``unique_Q``. This can happen because of numerical precision issues.

.. note:: in the display of the ``lf`` below, the “length” parameter is not the ENS. It is, instead, just a scalar.

.. jupyter-execute::

    fitted.lf

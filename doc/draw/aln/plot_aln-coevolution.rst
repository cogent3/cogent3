.. jupyter-execute::
    :hide-code:

    import set_working_directory

Coevolution analysis
====================

A method on the alignment provides an interface to the simple (and yet robust and fast) methods for estimating coevolution. The default measure is normalised mutual information (NMI).

.. todo:: add citation for RMI

Display coevolution as a heatmap
--------------------------------

Using the ``drawable`` argument causes the returned object to have a ``drawable`` attribute (type ``Drawable`` which has ``show()`` and ``write()`` methods), for the corresponding plot types -- a heatmap in this case.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
    aln = aln.no_degenerates(motif_length=3)
    aln = aln.get_translation()
    aln = aln[:100]  # sliced to simplify the visual display
    coevo = aln.coevolution(stat="nmi", show_progress=False, drawable="heatmap")
    coevo.drawable.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-coevolution.png"

    coevo.drawable.write(outpath)

Display coevolution scores as a Violin plot
-------------------------------------------

.. jupyter-execute::

    coevo = aln.coevolution(stat="nmi", show_progress=False, drawable="violin")
    coevo.drawable.show(width=300)

Display coevolution scores as a Boxplot
---------------------------------------

.. jupyter-execute::

    coevo = aln.coevolution(stat="nmi", show_progress=False, drawable="box")
    coevo.drawable.show(width=300)

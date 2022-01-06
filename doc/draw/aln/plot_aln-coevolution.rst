.. jupyter-execute::
    :hide-code:

    import set_working_directory

Coevolution analysis
====================

A method on the alignment provides an interface to the simpler (and yet robust and fast) methods for estimating coevolution. The default measure is normalised mutual information (NMI).

.. todo:: add citation for NMI

Display coevolution as a heatmap
--------------------------------

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
    aln = aln.no_degenerates(motif_length=3)
    aln = aln.get_translation()
    aln = aln[:100]  # for compute speed in testing the documentation
    coevo = aln.coevolution(show_progress=False, drawable="heatmap")
    coevo.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-coevolution.png"

    coevo.drawable.write(outpath)

Display coevolution scores as a Violin plot
-------------------------------------------

.. jupyter-execute::

    coevo = aln.coevolution(show_progress=False, drawable="violin")
    coevo.show(width=300)

Display coevolution scores as a Boxplot
---------------------------------------

.. jupyter-execute::

    coevo = aln.coevolution(show_progress=False, drawable="box")
    coevo.show(width=300)

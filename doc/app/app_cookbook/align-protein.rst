.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using a protein model
=====================

We use apps to load unaligned DNA sequences and to translate them into amino acids.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    loader = get_app("load_unaligned", format="fasta")
    to_aa = get_app("translate_seqs")
    process = loader + to_aa
    seqs = process("data/SCA1-cds.fasta")

Protein alignment with default settings
---------------------------------------

The default setting for "protein" is a WG01 model.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    aa_aligner = get_app("progressive_align", "protein")
    aligned = aa_aligner(seqs)
    aligned

Specify a different distance measure for estimating the guide tree
------------------------------------------------------------------

The distance measures available are percent or paralinear.

.. note:: An estimated guide tree has its branch lengths scaled so they are consistent with usage in a codon model.

.. jupyter-execute::
    :raises:

    aa_aligner = get_app("progressive_align", "protein", distance="paralinear")
    aligned = aa_aligner(seqs)
    aligned

Alignment settings and file provenance are recorded in the ``info`` attribute
-----------------------------------------------------------------------------

.. jupyter-execute::
    :raises:

    aligned.info

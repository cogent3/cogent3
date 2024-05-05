.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using a codon model
===================

We load the unaligned sequences we will use in our examples.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    loader = get_app("load_unaligned", format="fasta")
    seqs = loader("data/SCA1-cds.fasta")

.. note:: We use an app loader, but since this is just a single file we could have used the ``cogent3.load_unaligned_seqs()`` function.

Codon alignment with default settings
-------------------------------------

The default settings will result in estimation of a guide tree (using percent identity between the sequences). The default "codon" model is MG94HKY.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    codon_aligner = get_app("progressive_align", "codon")
    aligned = codon_aligner(seqs)
    aligned

.. note:: If you specify ``unique_guides=True``, a guide tree will be estimated for every alignment.

Specify a different distance measure for estimating the guide tree
------------------------------------------------------------------

The distance measures available are the same as for the nucleotide case (percent, TN93 or paralinear).

.. note:: An estimated guide tree has its branch lengths scaled so they are consistent with usage in a codon model.

.. jupyter-execute::
    :raises:

    nt_aligner = get_app("progressive_align", "codon", distance="paralinear")
    aligned = nt_aligner(seqs)
    aligned

Providing a guide tree
----------------------

.. jupyter-execute::
    :raises:

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    codon_aligner = get_app("progressive_align", "codon", guide_tree=tree)
    aligned = codon_aligner(seqs)
    aligned

.. warning:: The guide tree must have branch lengths, otherwise a ``ValueError`` is raised.

Specifying the gap parameters
-----------------------------

.. jupyter-execute::
    :raises:

    codon_aligner = get_app("progressive_align",
        "codon", guide_tree=tree, indel_rate=0.001, indel_length=0.01
    )
    aligned = codon_aligner(seqs)
    aligned

Specifying the substitution model and parameters
------------------------------------------------

Any ``cogent3`` codon substitution model can be used. (See ``cogent3.available_models()``.)

.. jupyter-execute::
    :raises:

    codon_aligner = get_app("progressive_align",
        "CNFHKY", guide_tree=tree, param_vals=dict(omega=0.1, kappa=3)
    )
    aligned = codon_aligner(seqs)
    aligned

.. note:: If you provide parameter values, those must be consistent with the model definition.

Alignment settings and file provenance are recorded in the ``info`` attribute
-----------------------------------------------------------------------------

The parameters used to construct the alignment, including the guide tree and substitution model, are record in the alignment ``info`` attribute.

.. jupyter-execute::
    :raises:

    aligned.info

.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using a nucleotide model
------------------------

We load the unaligned sequences we will use in our examples.

.. jupyter-execute::

    from cogent3.app import io

    reader = io.load_unaligned(format="fasta")
    seqs = reader("data/SCA1-cds.fasta")

Nucleotide alignment with default settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default setting for "nucleotide" is a HKY85 model.

.. jupyter-execute::

    from cogent3.app.align import progressive_align

    nt_aligner = progressive_align("nucleotide")
    aligned = nt_aligner(seqs)
    aligned

Specify a different distance measure for estimating the guide tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the nucleotide case, you can use TN93 or paralinear.

.. jupyter-execute::

    nt_aligner = progressive_align("nucleotide", distance="TN93")
    aligned = nt_aligner(seqs)
    aligned

Providing a guide tree
^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    nt_aligner = progressive_align("nucleotide", guide_tree=tree)
    aligned = nt_aligner(seqs)
    aligned

.. note:: You can also specify ``unique_guides=True``, which means a guide tree will be estimated for every alignment.

Specifying the substitution model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use any nucleotide substitution model. For a list of all available, see ``cogent3.available_models()``.

.. jupyter-execute::

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    nt_aligner = progressive_align("F81", guide_tree=tree)
    aligned = nt_aligner(seqs)
    aligned

Alignment settings and file provenance are recorded in the ``info`` attribute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    aligned.info

.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using a nucleotide model
------------------------

We load the unaligned sequences we will use in our examples.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_unaligned", format="fasta")
    seqs = loader("data/SCA1-cds.fasta")

.. note:: We use an app loader, but since this is just a single file we could have used the ``cogent3.load_unaligned_seqs()`` function.

Nucleotide alignment with default settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default setting for "nucleotide" is a HKY85 model.

.. jupyter-execute::

    from cogent3 import get_app

    nt_aligner = get_app("progressive_align", "nucleotide")
    aligned = nt_aligner(seqs)
    aligned

.. note:: If you specify ``unique_guides=True``, a guide tree will be estimated for every alignment.

Specify a different distance measure for estimating the guide tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the nucleotide case, you can use TN93 or paralinear.

.. jupyter-execute::

    nt_aligner = get_app("progressive_align", "nucleotide", distance="TN93")
    aligned = nt_aligner(seqs)
    aligned

Providing a guide tree
^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    nt_aligner = get_app("progressive_align", "nucleotide", guide_tree=tree)
    aligned = nt_aligner(seqs)
    aligned

.. warning:: The guide tree must have branch lengths, otherwise a ``ValueError`` is raised.

Specifying the substitution model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use any ``cogent3`` nucleotide substitution model. For a list of all available, see ``cogent3.available_models()``.

.. jupyter-execute::

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    nt_aligner = get_app("progressive_align", "F81", guide_tree=tree)
    aligned = nt_aligner(seqs)
    aligned

Alignment settings and file provenance are recorded in the ``info`` attribute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    aligned.info

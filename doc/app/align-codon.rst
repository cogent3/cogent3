.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using a codon model
-------------------

We load the unaligned sequences we will use in our examples.

.. jupyter-execute::

    from cogent3.app import io

    reader = io.load_unaligned(format="fasta")
    seqs = reader("data/SCA1-cds.fasta")

Codon alignment with default settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default settings will result in estimation of a guide tree (using percent identity between the sequences). The default "codon" model is MG94HKY.

.. jupyter-execute::

    from cogent3.app.align import progressive_align

    codon_aligner = progressive_align("codon")
    aligned = codon_aligner(seqs)
    aligned

The parameters used to construct the alignment, including the guide tree and substitution model, are record in the ``info`` attribute.

.. jupyter-execute::

    aligned.info

.. note:: If can also specify ``unique_guides=True``, which means a guide tree will be estimated for every alignment.

Specify a different distance measure for estimating the guide tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distance measures available are the same as for the nucleotide case (percent, TN93 or paralinear). 

.. note:: An estimated guide tree has its branch lengths scaled so they are consistent with usage in a codon model.

.. jupyter-execute::

    nt_aligner = progressive_align("codon", distance="paralinear")
    aligned = nt_aligner(seqs)
    aligned

Providing a guide tree
^^^^^^^^^^^^^^^^^^^^^^

.. note:: The guide tree needs to have branch lengths, otherwise a ``ValueError`` is raised.

.. jupyter-execute::

    tree = "((Chimp:0.001,Human:0.001):0.0076,Macaque:0.01,((Rat:0.01,Mouse:0.01):0.02,Mouse_Lemur:0.02):0.01)"
    codon_aligner = progressive_align("codon", guide_tree=tree)
    aligned = codon_aligner(seqs)
    aligned

Specifying the gap parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    codon_aligner = progressive_align(
        "codon", guide_tree=tree, indel_rate=0.001, indel_length=0.01
    )
    aligned = codon_aligner(seqs)
    aligned

Specifying the substitution model and parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Any codon substitution model can be used. (See ``cogent3.available_models()``.) If you provide parameter values, those must be consistent with the model definition.

.. jupyter-execute::

    codon_aligner = progressive_align(
        "CNFHKY", guide_tree=tree, param_vals=dict(omega=0.1, kappa=3)
    )
    aligned = codon_aligner(seqs)
    aligned

Alignment settings and file provenance are recorded in the ``info`` attribute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    aligned.info

Using a protein model
---------------------

We load the unaligned sequences we will use in our examples and translate them.

.. doctest::

    >>> from cogent3.app import io, translate
    >>> reader = io.load_unaligned(format="fasta")
    >>> to_aa = translate.translate_seqs()
    >>> process = reader + to_aa
    >>> seqs = process("data/SCA1-cds.fasta")

Protein alignment with default settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default setting for "protein" is a WG01 model.

.. doctest::

    >>> from cogent3.app.align import progressive_align
    >>> aa_aligner = progressive_align("protein")
    >>> aligned = aa_aligner(seqs)
    >>> aligned
    6 x 825 protein alignment...

Specify a different distance measure for estimating the guide tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distance measures available are percent or paralinear.

.. note:: An estimated guide tree has its branch lengths scaled so they are consistent with usage in a codon model.

.. doctest::
    
    >>> aa_aligner = progressive_align("protein", distance="paralinear")
    >>> aligned = aa_aligner(seqs)
    >>> aligned
    6 x 825 protein alignment...

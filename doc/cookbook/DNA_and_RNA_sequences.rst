.. _dna-rna-seqs:

Sequences
---------

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

The ``Sequence`` object provides generic biological sequence manipulation functions, plus functions that are critical for the ``evolve`` module calculations.

Generic molecular types
^^^^^^^^^^^^^^^^^^^^^^^

Sequence properties are affected by the moltype you specify. The default type for a sequence is ``"text"``.

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq.moltype.label

.. jupyter-execute::

    my_seq

In some circumstances you can also have a ``"bytes"`` moltype, which I'll explicitly construct here.

.. jupyter-execute::

    my_seq = make_seq("AGTACACTGGT", moltype="bytes")
    my_seq.moltype.label

.. jupyter-execute::

    my_seq

DNA and RNA sequences
^^^^^^^^^^^^^^^^^^^^^

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Tony Walters, Meg Pirrung

Creating a DNA sequence from a string
"""""""""""""""""""""""""""""""""""""

Sequence properties are affected by the moltype you specify. Here we specify the ``DNA`` molecular type.

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq

Creating a RNA sequence from a string
"""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    rnaseq = make_seq("ACGUACGUACGUACGU", moltype="rna")

Converting to FASTA format
""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq

Convert a RNA sequence to FASTA format
""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    rnaseq = make_seq("ACGUACGUACGUACGU", moltype="rna")
    rnaseq

Creating a named sequence
"""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", "my_gene", moltype="dna")
    my_seq
    type(my_seq)

Setting or changing the name of a sequence
""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq.name = "my_gene"
    my_seq

Complementing a DNA sequence
""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq.complement()

Reverse complementing a DNA sequence
""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    my_seq.rc()

.. _translation:

Translate a sequence to protein
"""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("GCTTGGGAAAGTCAAATGGAA", name="s1", moltype="dna")
    pep = my_seq.get_translation()
    type(pep)

.. jupyter-execute::

    pep

The default is to trim a terminating stop if it exists. If you set ``trim_stop=False`` and there is a terminating stop, an ``AlphabetError`` is raised.

.. jupyter-execute::
    :hide-code:

    from cogent3.core.new_alphabet import AlphabetError

.. jupyter-execute::
    :raises: AlphabetError

    from cogent3 import make_seq

    my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    my_seq.get_translation(trim_stop=False)

You can also specify the :ref:`genetic code <genetic-codes>`.

.. jupyter-execute::

    my_seq.get_translation(gc="Vertebrate Mitochondrial") # or gc=2

Translating a DNA sequence containing stop codons
"""""""""""""""""""""""""""""""""""""""""""""""""

By default, ``get_translation()`` will fail if there are any stop codons in frame in the sequence. You can allow translation in these cases by setting the optional argument ``include_stop=True``.

.. jupyter-execute::

    from cogent3 import make_seq

    seq = make_seq("ATGTGATGGTAA", name="s1", moltype="dna")
    pep = seq.get_translation(include_stop=True)
    pep

Converting a DNA sequence to RNA
""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("ACGTACGTACGTACGT", moltype="dna")
    rnaseq = my_seq.to_rna()
    rnaseq

Convert an RNA sequence to DNA
""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    rnaseq = make_seq("ACGUACGUACGUACGU", moltype="rna")
    dnaseq = rnaseq.to_dna()
    dnaseq

Testing complementarity
"""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    a = make_seq("AGTACACTGGT", moltype="dna")
    a.can_pair(a.complement())

.. jupyter-execute::

    a.can_pair(a.rc())

Joining two DNA sequences
"""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    extra_seq = make_seq("CTGAC", moltype="dna")
    long_seq = my_seq + extra_seq
    long_seq

Getting all *k*-mers from a sequence
""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    list(my_seq.iter_kmers(k=2))

.. note:: By default, any *k*-mer that contains an ambiguity code is excluded from the output.

You can include ALL *k*-mers by setting ``strict=False``.

.. jupyter-execute::

    my_seq = make_seq("AGTANACTGGT", moltype="dna")
    list(my_seq.iter_kmers(k=2, strict=False))

Slicing DNA sequences
"""""""""""""""""""""

.. jupyter-execute::

    my_seq[1:6]

Obtaining the codons from a ``DnaSequence`` object
""""""""""""""""""""""""""""""""""""""""""""""""""

Use the method ``get_in_motif_size``

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    codons = my_seq.get_in_motif_size(3)
    codons

Getting 3rd positions from codons
"""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    seq = make_seq("ATGATGATGATG", moltype="dna")
    pos3 = seq[2::3]
    assert str(pos3) == "GGGG"

Getting 1st and 2nd positions from codons
"""""""""""""""""""""""""""""""""""""""""

In this instance we can use features.

.. jupyter-execute::

    from cogent3 import make_seq

    seq = make_seq("ATGATGATGATG", moltype="dna")
    indices = [(i, i + 2) for i in range(len(seq))[::3]]
    pos12 = seq.add_feature(biotype="pos12", name="pos12", spans=indices)
    pos12 = pos12.get_slice()
    assert str(pos12) == "ATATATAT"

Return a randomised version of the sequence
"""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    rnaseq.shuffle()

Remove gaps from a sequence
"""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s = make_seq("--AUUAUGCUAU-UAU--", moltype="rna")
    s.degap()


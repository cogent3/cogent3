.. _dna-rna-seqs:

Sequences
---------

The ``Sequence`` object provides generic biological sequence manipulation functions, plus functions that are critical for the ``evolve`` module calculations.

Generic molecular types
^^^^^^^^^^^^^^^^^^^^^^^

Sequence properties are affected by the moltype you specify. The default type for a sequence is ``"text"``.

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT")
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

Sequence properties are affected by the moltype you specify. Here we specify the ``DNA`` ``MolType``.

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

You can also use a convenience ``make_seq()`` function, providing the moltype as a string.

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

Translate a ``DnaSequence`` to protein
""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("GCTTGGGAAAGTCAAATGGAA", name="s1", moltype="dna")
    pep = my_seq.get_translation()
    type(pep)

.. jupyter-execute::

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

Slicing DNA sequences
"""""""""""""""""""""

.. jupyter-execute::

    my_seq[1:6]

Getting 3rd positions from codons
"""""""""""""""""""""""""""""""""

The easiest approach is to work off the ``cogent3`` ``ArrayAlignment`` object.

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

Return a randomized version of the sequence
"""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    rnaseq.shuffle()

Remove gaps from a sequence
"""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s = make_seq("--AUUAUGCUAU-UAu--", moltype="rna")
    s.degap()

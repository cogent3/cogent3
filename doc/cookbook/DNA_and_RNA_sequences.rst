.. _dna-rna-seqs:

``Sequence``
============

The ``Sequence`` object contains classes that represent biological sequence data. These provide generic biological sequence manipulation functions, plus functions that are critical for the ``evolve`` module calculations.

.. warning:: Do not import sequence classes directly! It is expected that you will access them through ``MolType`` objects. The molecular types can be accessed via the ``cogent3.get_moltype()`` function. Sequence classes depend on information from the ``MolType`` that is **only** available after ``MolType`` has been imported. Sequences are intended to be immutable. This is not enforced by the code for performance reasons, but don't alter the ``MolType`` or the sequence data after creation.

DNA and RNA sequences
---------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Tony Walters, Meg Pirrung

Creating a DNA sequence from a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All sequence and alignment objects have a molecular type, or ``MolType`` which provides key properties for validating sequence characters. Here we use the ``DNA`` ``MolType`` to create a DNA sequence.

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("AGTACACTGGT")
    my_seq
    print(my_seq)
    str(my_seq)

Creating a RNA sequence from a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import RNA

    rnaseq = RNA.make_seq("ACGUACGUACGUACGU")

Converting to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("AGTACACTGGT")
    print(my_seq.to_fasta())

Convert a RNA sequence to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import RNA

    rnaseq = RNA.make_seq("ACGUACGUACGUACGU")
    rnaseq.to_fasta()

Creating a named sequence
^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use a convenience ``make_seq()`` function, providing the moltype as a string.

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", "my_gene", moltype="dna")
    my_seq
    type(my_seq)

Setting or changing the name of a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGT", moltype="dna")
    my_seq.name = "my_gene"
    print(my_seq.to_fasta())

Complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("AGTACACTGGT")
    print(my_seq.complement())

Reverse complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    print(my_seq.rc())

The ``rc`` method name is easier to type

.. jupyter-execute::

    print(my_seq.rc())

.. _translation:

Translate a ``DnaSequence`` to protein
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("GCTTGGGAAAGTCAAATGGAA", "protein-X")
    pep = my_seq.get_translation()
    type(pep)
    print(pep.to_fasta())

Converting a DNA sequence to RNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("ACGTACGTACGTACGT")
    print(my_seq.to_rna())

Convert an RNA sequence to DNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import RNA

    rnaseq = RNA.make_seq("ACGUACGUACGUACGU")
    print(rnaseq.to_dna())

Testing complementarity
^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    a = DNA.make_seq("AGTACACTGGT")
    a.can_pair(a.complement())
    a.can_pair(a.rc())

Joining two DNA sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import DNA

    my_seq = DNA.make_seq("AGTACACTGGT")
    extra_seq = DNA.make_seq("CTGAC")
    long_seq = my_seq + extra_seq
    long_seq
    str(long_seq)

Slicing DNA sequences
^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    my_seq[1:6]

Getting 3rd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest approach is to work off the ``cogent3`` ``ArrayAlignment`` object.

.. jupyter-execute::

    from cogent3 import DNA

    seq = DNA.make_seq("ATGATGATGATG")
    pos3 = seq[2::3]
    assert str(pos3) == "GGGG"

Getting 1st and 2nd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this instance we can use the annotatable sequence classes.

.. jupyter-execute::

    from cogent3 import DNA

    seq = DNA.make_seq("ATGATGATGATG")
    indices = [(i, i + 2) for i in range(len(seq))[::3]]
    pos12 = seq.add_feature(biotype="pos12", name="pos12", spans=indices)
    pos12 = pos12.get_slice()
    assert str(pos12) == "ATATATAT"

Return a randomized version of the sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   print rnaseq.shuffle()
   ACAACUGGCUCUGAUG

Remove gaps from a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import RNA

    s = RNA.make_seq("--AUUAUGCUAU-UAu--")
    print(s.degap())

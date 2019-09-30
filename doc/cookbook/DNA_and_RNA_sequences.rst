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

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq("AGTACACTGGT")
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> print(my_seq)
    AGTACACTGGT
    >>> str(my_seq)
    'AGTACACTGGT'

Creating a RNA sequence from a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import RNA
    >>> rnaseq = RNA.make_seq('ACGUACGUACGUACGU')

Converting to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq('AGTACACTGGT')
    >>> print(my_seq.to_fasta())
    >0
    AGTACACTGGT
    <BLANKLINE>

Convert a RNA sequence to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import RNA
    >>> rnaseq = RNA.make_seq('ACGUACGUACGUACGU')
    >>> rnaseq.to_fasta()
    '>0\nACGUACGUACGUACGU\n'

Creating a named sequence
^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use a convenience ``make_seq()`` function, providing the moltype as a string.

.. doctest::

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq('AGTACACTGGT','my_gene', moltype="dna")
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> type(my_seq)
    <class 'cogent3.core.sequence.DnaSequence'>

Setting or changing the name of a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq('AGTACACTGGT', moltype="dna")
    >>> my_seq.name = 'my_gene'
    >>> print(my_seq.to_fasta())
    >my_gene
    AGTACACTGGT
    <BLANKLINE>

Complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq("AGTACACTGGT")
    >>> print(my_seq.complement())
    TCATGTGACCA

Reverse complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print(my_seq.rc())
    ACCAGTGTACT

The ``rc`` method name is easier to type

.. doctest::

    >>> print(my_seq.rc())
    ACCAGTGTACT

.. _translation:

Translate a ``DnaSequence`` to protein
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq('GCTTGGGAAAGTCAAATGGAA','protein-X')
    >>> pep = my_seq.get_translation()
    >>> type(pep)
    <class 'cogent3.core.sequence.ProteinSequence'>
    >>> print(pep.to_fasta())
    >protein-X
    AWESQME
    <BLANKLINE>

Converting a DNA sequence to RNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq('ACGTACGTACGTACGT')
    >>> print(my_seq.to_rna())
    ACGUACGUACGUACGU

Convert an RNA sequence to DNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import RNA
   >>> rnaseq = RNA.make_seq('ACGUACGUACGUACGU')
   >>> print(rnaseq.to_dna())
   ACGTACGTACGTACGT

Testing complementarity
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> a = DNA.make_seq("AGTACACTGGT")
    >>> a.can_pair(a.complement())
    False
    >>> a.can_pair(a.rc())
    True

Joining two DNA sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import DNA
    >>> my_seq = DNA.make_seq("AGTACACTGGT")
    >>> extra_seq = DNA.make_seq("CTGAC")
    >>> long_seq = my_seq + extra_seq
    >>> long_seq
    DnaSequence(AGTACAC... 16)
    >>> str(long_seq)
    'AGTACACTGGTCTGAC'

Slicing DNA sequences
^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> my_seq[1:6]
    DnaSequence(GTACA)

Getting 3rd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest approach is to work off the ``cogent3`` ``ArrayAlignment`` object.

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent3 import DNA
    >>> seq = DNA.make_array_seq('ATGATGATGATG')
    >>> pos3 = seq[2::3]
    >>> assert str(pos3) == 'GGGG'

Getting 1st and 2nd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this instance we can use the annotatable sequence classes.

.. doctest::

    >>> from cogent3 import DNA
    >>> seq = DNA.make_seq('ATGATGATGATG')
    >>> indices = [(i, i+2) for i in range(len(seq))[::3]]
    >>> pos12 = seq.add_feature('pos12', 'pos12', indices)
    >>> pos12 = pos12.get_slice()
    >>> assert str(pos12) == 'ATATATAT'

Return a randomized version of the sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   print rnaseq.shuffle()
   ACAACUGGCUCUGAUG

Remove gaps from a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import RNA
   >>> s = RNA.make_seq('--AUUAUGCUAU-UAu--')
   >>> print(s.degap())
   AUUAUGCUAUUAU
  
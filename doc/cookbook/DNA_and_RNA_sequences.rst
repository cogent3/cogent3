.. _dna-rna-seqs:

DNA and RNA sequences
---------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Tony Walters, Meg Pirrung

Creating a DNA sequence from a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All sequence and alignment objects have a molecular type, or ``MolType`` which provides key properties for validating sequence characters. Here we use the ``DNA`` ``MolType`` to create a DNA sequence.

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> print my_seq
    AGTACACTGGT
    >>> str(my_seq)
    'AGTACACTGGT'

Creating a RNA sequence from a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import RNA
    >>> rnaseq = RNA.makeSequence('ACGUACGUACGUACGU')

Converting to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT')
    >>> print my_seq.toFasta()
    >0
    AGTACACTGGT

Convert a RNA sequence to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import RNA
    >>> rnaseq = RNA.makeSequence('ACGUACGUACGUACGU')
    >>> rnaseq.toFasta()
    '>0\nACGUACGUACGUACGU'

Creating a named sequence
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT','my_gene')
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> type(my_seq)
    <class 'cogent.core.sequence.DnaSequence'>

Setting or changing the name of a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT')
    >>> my_seq.Name = 'my_gene'
    >>> print my_seq.toFasta()
    >my_gene
    AGTACACTGGT

Complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> print my_seq.complement()
    TCATGTGACCA

Reverse complementing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print my_seq.reversecomplement()
    ACCAGTGTACT

The ``rc`` method name is easier to type

.. doctest::

    >>> print my_seq.rc()
    ACCAGTGTACT

.. _translation:

Translate a ``DnaSequence`` to protein
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('GCTTGGGAAAGTCAAATGGAA','protein-X')
    >>> pep = my_seq.getTranslation()
    >>> type(pep)
    <class 'cogent.core.sequence.ProteinSequence'>
    >>> print pep.toFasta()
    >protein-X
    AWESQME

Converting a DNA sequence to RNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('ACGTACGTACGTACGT')
    >>> print my_seq.toRna()
    ACGUACGUACGUACGU

Convert an RNA sequence to DNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import RNA
   >>> rnaseq = RNA.makeSequence('ACGUACGUACGUACGU')
   >>> print rnaseq.toDna()
   ACGTACGTACGTACGT

Testing complementarity
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> a = DNA.makeSequence("AGTACACTGGT")
    >>> a.canPair(a.complement())
    False
    >>> a.canPair(a.reversecomplement())
    True

Joining two DNA sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> extra_seq = DNA.makeSequence("CTGAC")
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

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')

Creating the position indices, note that we start at the 2nd index (the 'first' codon's 3rd position) indicate each position as a *span* (``i -- i+1``).

.. doctest::

    >>> indices = [(i, i+1) for i in range(len(seq))[2::3]]

Create the sequence feature and use it to slice the sequence.

.. doctest::

    >>> pos3 = seq.addFeature('pos3', 'pos3', indices)
    >>> pos3 = pos3.getSlice()
    >>> assert str(pos3) == 'GGGG'

Getting 1st and 2nd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The only difference here to above is that our spans cover 2 positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')
    >>> indices = [(i, i+2) for i in range(len(seq))[::3]]
    >>> pos12 = seq.addFeature('pos12', 'pos12', indices)
    >>> pos12 = pos12.getSlice()
    >>> assert str(pos12) == 'ATATATAT'

Return a randomized version of the sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   print rnaseq.shuffle()
   ACAACUGGCUCUGAUG

Remove gaps from a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import RNA
   >>> s = RNA.makeSequence('--AUUAUGCUAU-UAu--')
   >>> print s.degap()
   AUUAUGCUAUUAU

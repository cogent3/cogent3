Alphabets
---------

.. authors Gavin Huttley

``Alphabet`` and ``MolType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MolType`` instances have an ``Alphabet``.

.. doctest::
    
    >>> from cogent import DNA, PROTEIN
    >>> print DNA.Alphabet
    ('T', 'C', 'A', 'G')
    >>> print PROTEIN.Alphabet
    ('A', 'C', 'D', 'E', ...

``Alphabet`` instances have a ``MolType``.

.. doctest::
    
    >>> PROTEIN.Alphabet.MolType == PROTEIN
    True

Creating tuple alphabets
^^^^^^^^^^^^^^^^^^^^^^^^

You can create a tuple alphabet of, for example, dinucleotides or trinucleotides.

.. doctest::
    
    >>> dinuc_alphabet = DNA.Alphabet.getWordAlphabet(2)
    >>> print dinuc_alphabet
    ('TT', 'CT', 'AT', 'GT', ...
    >>> trinuc_alphabet = DNA.Alphabet.getWordAlphabet(3)
    >>> print trinuc_alphabet
    ('TTT', 'CTT', 'ATT', ...

Convert a sequence into integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> seq = 'TAGT'
    >>> indices = DNA.Alphabet.toIndices(seq)
    >>> indices
    [0, 2, 3, 0]

Convert integers to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> seq = DNA.Alphabet.fromIndices([0,2,3,0])
    >>> seq
    ['T', 'A', 'G', 'T']

or

.. doctest::
    
    >>> seq = DNA.Alphabet.fromOrdinalsToSequence([0,2,3,0])
    >>> seq
    DnaSequence(TAGT)

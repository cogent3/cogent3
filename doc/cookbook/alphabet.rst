Alphabets
---------

.. authors Gavin Huttley

``Alphabet`` and ``MolType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MolType`` instances have an ``Alphabet``.

.. doctest::

    >>> from cogent3 import DNA, PROTEIN
    >>> print(DNA.alphabet)
    ('T', 'C', 'A', 'G')
    >>> print(PROTEIN.alphabet)
    ('A', 'C', 'D', 'E', ...

``Alphabet`` instances have a ``MolType``.

.. doctest::

    >>> PROTEIN.alphabet.moltype == PROTEIN
    True

Creating tuple alphabets
^^^^^^^^^^^^^^^^^^^^^^^^

You can create a tuple alphabet of, for example, dinucleotides or trinucleotides.

.. doctest::

    >>> dinuc_alphabet = DNA.alphabet.get_word_alphabet(2)
    >>> print(dinuc_alphabet)
    ('TT', 'CT', 'AT', 'GT', ...
    >>> trinuc_alphabet = DNA.alphabet.get_word_alphabet(3)
    >>> print(trinuc_alphabet)
    ('TTT', 'CTT', 'ATT', ...

Convert a sequence into integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> seq = 'TAGT'
    >>> indices = DNA.alphabet.to_indices(seq)
    >>> indices
    [0, 2, 3, 0]

Convert integers to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> seq = DNA.alphabet.from_indices([0,2,3,0])
    >>> seq
    ['T', 'A', 'G', 'T']

or

.. doctest::

    >>> seq = DNA.alphabet.from_ordinals_to_seq([0,2,3,0])
    >>> seq
    DnaSequence(TAGT)

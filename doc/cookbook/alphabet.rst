Alphabets
---------

.. authors Gavin Huttley

``Alphabet`` and ``MolType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MolType`` instances have an ``Alphabet``.

.. jupyter-execute::
    :linenos:

    from cogent3 import DNA, PROTEIN

    print(DNA.alphabet)
    print(PROTEIN.alphabet)

``Alphabet`` instances have a ``MolType``.

.. jupyter-execute::
    :linenos:

    PROTEIN.alphabet.moltype == PROTEIN

Creating tuple alphabets
^^^^^^^^^^^^^^^^^^^^^^^^

You can create a tuple alphabet of, for example, dinucleotides or trinucleotides.

.. jupyter-execute::
    :linenos:

    dinuc_alphabet = DNA.alphabet.get_word_alphabet(2)
    print(dinuc_alphabet)
    trinuc_alphabet = DNA.alphabet.get_word_alphabet(3)
    print(trinuc_alphabet)

Convert a sequence into integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :linenos:

    seq = "TAGT"
    indices = DNA.alphabet.to_indices(seq)
    indices

Convert integers to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :linenos:

    seq = DNA.alphabet.from_indices([0, 2, 3, 0])
    seq

or

.. jupyter-execute::
    :linenos:

    seq = DNA.alphabet.from_ordinals_to_seq([0, 2, 3, 0])
    seq

Alphabets
---------

.. authors Gavin Huttley

``Alphabet`` and ``MolType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MolType`` instances have an ``Alphabet``.

.. jupyter-execute::

    from cogent3 import get_moltype

    dna = get_moltype("dna")
    print(dna.alphabet)
    protein = get_moltype("protein")
    print(protein.alphabet)

``Alphabet`` instances reference the ``MolType`` that created them.

.. jupyter-execute::

    dna.alphabet.moltype is dna

Creating tuple alphabets
^^^^^^^^^^^^^^^^^^^^^^^^

You can create a tuple alphabet of, for example, dinucleotides or trinucleotides.

.. jupyter-execute::

    dinuc_alphabet = dna.alphabet.get_word_alphabet(2)
    print(dinuc_alphabet)
    trinuc_alphabet = dna.alphabet.get_word_alphabet(3)
    print(trinuc_alphabet)

Convert a sequence into integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    seq = "TAGT"
    indices = dna.alphabet.to_indices(seq)
    indices

Convert integers to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    seq = dna.alphabet.from_indices([0, 2, 3, 0])
    seq


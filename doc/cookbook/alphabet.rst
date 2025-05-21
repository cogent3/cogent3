.. jupyter-execute::
    :hide-code:

    import os

    # new type only from now
    os.environ["COGENT3_NEW_TYPE"] = "1"

Alphabets
---------

.. authors Gavin Huttley

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

``CharAlphabet`` and ``MolType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``MolType`` instances have an ``CharAlphabet``.

.. jupyter-execute::

    from cogent3 import get_moltype

    dna = get_moltype("dna")
    print(dna.alphabet)
    protein = get_moltype("protein")
    print(protein.alphabet)

``CharAlphabet`` instances reference the ``MolType`` that created them.

.. jupyter-execute::

    dna.alphabet.moltype is dna

Creating tuple alphabets
^^^^^^^^^^^^^^^^^^^^^^^^

You can create a tuple alphabet of, for example, dinucleotides or trinucleotides.

.. jupyter-execute::

    dinuc_alphabet = dna.alphabet.get_kmer_alphabet(2)
    trinuc_alphabet = dna.alphabet.get_kmer_alphabet(3)
    print(dinuc_alphabet, trinuc_alphabet, sep="\n")

Convert a sequence into integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    seq = "TAGT"
    indices = dna.alphabet.to_indices(seq)
    indices

Convert integers to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    import numpy
    
    data = numpy.array([0, 2, 3, 0], dtype=numpy.uint8)
    seq = dna.alphabet.from_indices(data)
    seq

Converting a sequence into k-mer indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use a ``KmerAlphabet`` to convert a standard sequence into a ``numpy`` array of integers. In this case, each integer is the encoding of the dinucleotide string into the index of that dinucleotide. Because the ``CharAlphabet`` and ``KmerAlphabet`` both inherit from  ``tuple``, they have the built-in ``.index()`` method.

.. jupyter-execute::

    import numpy

    bases = dna.alphabet
    dinucs = bases.get_kmer_alphabet(2)
    dinucs.index("TG")    

The ``to_indices()`` method is faster and provides more flexibility. We use that on the single dinucleotide

.. jupyter-execute::

    dinucs.to_indices("TG")

and on a longer sequence where we want the independent *k*-mers.

.. jupyter-execute::

    seq = "TGTGGCACAAATACTCATGCCAGCTCATTA"
    dinuc_indices = dinucs.to_indices(seq, independent_kmer=True)
    dinuc_indices

We can also convert the sequence into all possible k-mers.

.. jupyter-execute::

    dinucs.to_indices(seq, independent_kmer=False)

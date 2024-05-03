.. jupyter-execute::
    :hide-code:

    import set_working_directory

Select `n` sequences from a collection
--------------------------------------

Let's load alignment of primates to use in examples. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    loader = get_app("load_aligned", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    aln

Select the first `n` seuquences from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Initilising ``take_n_seqs`` with the argument ``number=3`` creates an app which returns the first 3 sequences from an alignment 

.. note::  "first n" refers to the ordering in the fasta file. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    first_3 = get_app("take_n_seqs", number=3)
    first_3(aln)

Randomly selecting `n` seuquences from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using ``random=True`` and ``number=3`` returns 3 random sequences. An optional argument for a ``seed`` can be provided to ensure the same sequences are returned each time the app is called.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    random_n = get_app("take_n_seqs", random=True, number=3, seed=1)
    random_n(aln)

Selecting the same sequences from multiple alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Providing the argument ``fixed_choice=True`` ensures the same sequences are returned when (randomly) sampling sequences across several alignments.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    loader = get_app("load_aligned", moltype="dna")
    aln1 = loader("data/primate_brca1.fasta")
    aln2 = loader("data/brca1.fasta")

    aln1.names

.. jupyter-execute::
    :raises:

    aln2.names

.. jupyter-execute::
    :raises:

    fixed_choice = get_app("take_n_seqs", number=2, random=True, fixed_choice=True)
    result1 = fixed_choice(aln1).names
    result2 = fixed_choice(aln2).names
    result1 == result2
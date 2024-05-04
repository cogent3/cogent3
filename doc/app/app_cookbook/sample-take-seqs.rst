.. jupyter-execute::
    :hide-code:

    import set_working_directory

Select named sequences from a collection
----------------------------------------

Let's load alignment of primates to use in examples. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    loader = get_app("load_aligned", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    aln

Selecting sequences from a collection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the ``take_named_seqs`` app we can select sequences from a collection that match provided names. For instance, we could be interested in an alignment of sequences from Humans with those from our two closest relatives. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    greater_ape_selector = get_app("take_named_seqs", "Human", "Chimpanzee", "Gorilla")
    ape_aln = greater_ape_selector(aln)
    ape_aln

Using ``take_named_seqs`` in a composed app process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    select_seqs = get_app("take_named_seqs", "Human", "Rhesus", "Galago")
    process = loader + select_seqs
    hrg_aln = process("data/primate_brca1.fasta")
    hrg_aln.names

Removing sequences from a collection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating the ``take_named_seqs`` app with the argument ``negate=True`` will exclude sequences that match the provided names. For instance, perhaps the Galago has got to go. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    no_galagos = get_app("take_named_seqs", "Galago", negate=True)
    no_galago_aln = no_galagos(hrg_aln)
    no_galago_aln.names

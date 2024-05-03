.. jupyter-execute::
    :hide-code:

    import set_working_directory

Sample nucleotides from a given codon position
----------------------------------------------

The ``take_codon_positions`` app allows you to extract all nucleotides at a given codon position from an alignment. 

Let's create a sample alignment for our example. 

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs({"s1": "ACGACGACG", "s2": "GATGATGAT"}, moltype="dna")
    aln

Extract the third codon position from an alignment 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can achieve this by creating the ``take_codon_positions`` app with ``3`` as a positional argument.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    take_pos3 = get_app("take_codon_positions", 3, moltype="dna")
    result = take_pos3(aln)
    result

Extract the first and second codon positions from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can achieve this by creating the ``take_codon_positions`` app with ``1`` and ``2`` as a positional argument. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    take_pos12 = get_app("take_codon_positions", 1, 2, moltype="dna")
    result = take_pos12(aln)
    result

Extract only the third codon positions from four-fold degenerate codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can achieve this by creating the ``take_codon_positions`` app with the argument ``fourfold_degenerate=True``. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, make_aligned_seqs

    aln_ff = make_aligned_seqs({"s1": "GCAAGCGTTTAT", "s2": "GCTTTTGTCAAT"})
    take_fourfold = get_app("take_codon_positions", fourfold_degenerate=True, moltype="dna")
    result = take_fourfold(aln_ff)
    result

Create a composed process which samples only the third codon position
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :hide-code:

    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

Let's set up a data store containing all the files with the ".fasta" suffix in the data directory, limiting the data store to two members as a minimum example.

.. jupyter-execute::
    :raises:

    from cogent3 import open_data_store

    fasta_seq_dstore = open_data_store("data", suffix="fasta", mode="r", limit=2)

Now let's set up a process composing the following apps: ``load_aligned`` (loads the sequences ), ``take_codon_positions`` (extracts the third codon position), and ``write_seqs`` (writes the filtered sequences to a data store). 

.. note:: Learn the basics of turning apps into composed processes :ref:`here! <apps>` 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, suffix="fa", mode="w")

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    cpos3 = get_app("take_codon_positions", 3)
    writer = get_app("write_seqs", out_dstore, format="fasta")

    process = loader + cpos3 + writer

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.

Now let's apply ``process`` to our data store! ``result`` is a data store containing the filtered alignments, which we can index to see individual data members. We could take a closer look using the ``.read()`` method on data members. 

.. jupyter-execute::
    :raises:

    result = process.apply_to(fasta_seq_dstore)
    result.describe

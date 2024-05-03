.. jupyter-execute::
    :hide-code:

    import set_working_directory

Loading aligned sequence data
-----------------------------

We can load aligned sequence data using the ``load_aligned`` app. When making the app, you can optionally provide arguments for the molecular type of the sequence and the format of the data. 

Loading aligned DNA sequences from a single fasta file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we load the brca1 gene in bats, providing the molecular type (``moltype="dna"``) and file format (``format="fasta"``). 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app
    
    load_aligned_app = get_app("load_aligned", moltype="dna", format="fasta")
    aln = load_aligned_app("data/brca1-bats.fasta")
    aln

Loading aligned protien sequences from a single phylip file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we load a globin alignment, providing the molecular type (``moltype="protein"``) and file format (``format="phylip"``). 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app
    
    load_aligned_app = get_app("load_aligned", moltype="protein", format="phylip")
    aln = load_aligned_app("data/abglobin_aa.phylip")
    aln
    
Loading aligned DNA sequences from multiple fasta files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above examples, the result is a single alignment, which could have been achieved using standard cogent3 (``load_aligned_seqs()``). The real power of apps is for batch processing of a large number of files.

To apply apps to multiple files we need to set two things up:

.. jupyter-execute::
    :hide-code:

    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

1. A data store that identifies the files we are interested in 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Here, we create a data store containing all the files with the ".fasta" suffix in the data directory, limiting the data store to two members as a minimum example.

.. jupyter-execute::
    :raises:

    from cogent3 import open_data_store

    fasta_seq_dstore = open_data_store("data", suffix="fasta", mode="r", limit=2)

2. A composed process that defines our workflow 
"""""""""""""""""""""""""""""""""""""""""""""""

In this example, our process loads the sequences, filters the sequences to keep only the third codon position, and then writes the filtered sequences to a data store. 

.. note:: Apps that are "writers" typically require a data store to write to, learn more about writers :ref:`here! <writers>` 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, suffix="fa", mode="w")

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    cpos3 = get_app("take_codon_positions", 3)
    writer = get_app("write_seqs", out_dstore, format="fasta")

    process = loader + cpos3 + writer

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.

Now we're good to go, we can apply ``process`` to our data store!
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``result`` is a data store, which you can index to see individual data members - which are our alignments. We can take a closer look using the ``.read()`` method on data members (truncating to 50 characters). 

.. jupyter-execute::
    :raises:

    result = process.apply_to(fasta_seq_dstore)
    print(result[0].read()[:50])

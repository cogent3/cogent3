.. jupyter-execute::
    :hide-code:

    import set_working_directory

Filter sequence collections and alignments by length
----------------------------------------------------

Let's load a collection of globin sequences. Note that we must have a the molecular type specified. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app

    loader = get_app("load_unaligned", moltype="dna", format="fasta")
    aln = loader("data/SCA1-cds.fasta")
    aln

Remove sequences shorter than a minimum length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating the ``min_length`` app and providing a positional argument specifying the minimum length allows us to filter an alignment, removing sequences which do not satisfy our threshold. 

For instance, we can remove sequences from the globin alignment which are shorter than 240 amino acids long. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app

    over_240s = get_app("min_length", 240)
    aln_fltrd = over_240s(aln)
    aln_fltrd

Using the ``min_length`` length app to filter multiple alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :hide-code:

    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

``min_length`` is a useful app to use in data processing pipelines where downstream analysis requires sequences that exceed a given length.

In the following example, we compose a process that loads alignments and removes sequences less than 300 nucleotides in length, before writing them to a data store. We apply this process to a data store of the fasta files in the data directory. We restrict this data store to two members as a minimum example. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    dstore = open_data_store("data", suffix="fasta", limit=2)

    reader = get_app("load_aligned", format="fasta", moltype="dna")
    min_length = get_app("min_length", 300)
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="fasta")
    writer = get_app("write_seqs", out_dstore)

    process = reader + min_length + writer
    result = process.apply_to(dstore)
    result.describe
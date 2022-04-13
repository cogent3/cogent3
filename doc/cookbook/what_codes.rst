.. jupyter-execute::
    :hide-code:

    import set_working_directory

***********************
Available genetic codes
***********************

.. jupyter-execute::

    from cogent3 import available_codes

    available_codes()

In cases where a ``cogent3`` object method has a ``gc`` argument, you can just use the number under "Code ID" column.

For example:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    nt_seqs = load_aligned_seqs("data/brca1-bats.fasta", moltype="dna")
    nt_seqs[:21]

We specify the genetic code, and that codons that are incomplete as they contain a gap, are converted to ``?``.

.. jupyter-execute::

    aa_seqs = nt_seqs.get_translation(gc=1, incomplete_ok=True)
    aa_seqs[:20]

Getting a genetic code with ``get_code()``
==========================================

This function can be used directly to get a genetic code. We will get the code with ID 4.

.. jupyter-execute::

    from cogent3 import get_code

    gc = get_code(4)
    gc

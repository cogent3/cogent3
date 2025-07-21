.. jupyter-execute::
    :hide-code:

    import set_working_directory

Protein sequences
-----------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Creating a ProteinSequence with a name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_seq

    p = make_seq("THISISAPRQTEIN", "myProtein", moltype="protein")
    type(p)

.. jupyter-execute::

    p

Converting a DNA sequence string to protein sequence string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code.translate("TTTGCAAAC")

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Converting a nucleic acid sequence object to protein
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_seq
    
    nuc = make_seq("TTTGCAAAC", moltype="dna")
    pep = nuc.get_translation()
    pep

Loading protein sequences from a Phylip file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    seq = load_aligned_seqs("data/abglobin_aa.phylip", moltype="protein")

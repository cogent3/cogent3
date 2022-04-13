.. jupyter-execute::
    :hide-code:

    import set_working_directory

Protein sequences
-----------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Creating a ProteinSequence with a name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import PROTEIN

    p = PROTEIN.make_seq("THISISAPRQTEIN", "myProtein")
    type(p)
    str(p)

Converting a DNA sequence string to protein sequence string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3.core.genetic_code import DEFAULT as standard_code

    standard_code.translate("TTTGCAAAC")

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Loading protein sequences from a Phylip file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    seq = load_aligned_seqs("data/abglobin_aa.phylip", moltype="protein")

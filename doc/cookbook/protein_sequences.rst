.. jupyter-execute::
    :hide-code:

    import set_working_directory

Protein sequences
-----------------

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

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

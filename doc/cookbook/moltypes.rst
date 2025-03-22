.. jupyter-execute::
    :hide-code:

    import set_working_directory

***************
Molecular types
***************

.. note:: **Alpha Release of the New MolType API**

   We are pleased to announce an alpha release of our new ``MolType`` API. This version can be accessed by specifying the argument ``new_type=True`` in the ``get_moltype()`` function. 

   Please be aware that this alpha release has not been fully integrated with the library. Users are encouraged to explore its capabilities but should proceed with caution!

The ``MolType`` object provides services for resolving ambiguities, or providing the correct ambiguity for recoding. It also maintains the mappings between different kinds of alphabets, sequences and alignments.

If your analysis involves handling ambiguous states, or translation via a genetic code, it's critical to specify the appropriate moltype.

Available molecular types
=========================

.. jupyter-execute::

    from cogent3 import available_moltypes

    available_moltypes()

For statements that have a ``moltype`` argument, use the entry under the "Abbreviation" column. For example:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    seqs = load_aligned_seqs("data/brca1-bats.fasta", moltype="dna")

Getting a ``MolType``
=====================

.. jupyter-execute::

    from cogent3 import get_moltype

    dna = get_moltype("dna")
    dna

Using a ``MolType`` to get ambiguity codes
==========================================

Just using ``dna`` from above.

.. jupyter-execute::

    dna.ambiguities

Nucleic acid ``MolType`` and complementing
==========================================

.. jupyter-execute::

    dna.complement("AGG")

Making sequences
================

Use the either the top level ``cogent3.make_seq`` function, or the method on the ``MolType`` instance.

.. jupyter-execute::

    seq = dna.make_seq(seq="AGGCTT", name="seq1")
    seq

Verify sequences
================

.. jupyter-execute::

    rna = get_moltype("rna")
    rna.is_valid("ACGUACGUACGUACGU")

Making a custom ``MolType``
===========================

We demonstrate this by customising DNA so it allows ``.`` as gaps

.. jupyter-execute::

    from cogent3.core import moltype as mt

    DNAgapped = mt.MolType(
        seq_constructor=mt.DnaSequence,
        motifset=mt.IUPAC_DNA_chars,
        ambiguities=mt.IUPAC_DNA_ambiguities,
        complements=mt.IUPAC_DNA_ambiguities_complements,
        pairs=mt.DnaStandardPairs,
        gaps=".",
    )
    seq = DNAgapped.make_seq("ACG.")
    seq

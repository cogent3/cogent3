.. jupyter-execute::
    :hide-code:

    import set_working_directory

***************
Molecular types
***************

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

We demonstrate this by customising DNA so it allows a ``.`` as a gap character.

.. jupyter-execute::

    from cogent3.core import moltype, sequence

    mt = moltype.MolType(
            monomers="".join(moltype.IUPAC_DNA_chars),
            ambiguities=moltype.IUPAC_DNA_ambiguities,
            name="dna.gap",
            complements=moltype.IUPAC_DNA_ambiguities_complements,
            make_seq=sequence.DnaSequence,
            pairing_rules=moltype.DNA_STANDARD_PAIRS,
            mw_calculator=moltype.DnaMW,
            coerce_to=moltype.coerce_to_dna,
            gap=".",
        )
    seq = mt.make_seq(seq="ACG.")
    seq

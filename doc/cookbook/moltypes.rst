.. jupyter-execute::
    :hide-code:

    import set_working_directory

***************
Molecular types
***************

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

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

    from cogent3.core import new_moltype, new_sequence

    mt = new_moltype.MolType(
            monomers="".join(new_moltype.IUPAC_DNA_chars),
            ambiguities=new_moltype.IUPAC_DNA_ambiguities,
            name="dna.gap",
            complements=new_moltype.IUPAC_DNA_ambiguities_complements,
            make_seq=new_sequence.DnaSequence,
            pairing_rules=new_moltype.DNA_STANDARD_PAIRS,
            mw_calculator=new_moltype.DnaMW,
            coerce_to=new_moltype.coerce_to_dna,
            gap=".",
        )
    seq = mt.make_seq(seq="ACG.")
    seq

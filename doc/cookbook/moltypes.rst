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

``MolType`` definition of degenerate codes
==========================================

.. jupyter-execute::

    dna.degenerates

Nucleic acid ``MolType`` and complementing
==========================================

.. jupyter-execute::

    dna.complement("AGG")

Making sequences
================

Use the either the top level ``cogent3.make_seq`` function, or the method on the ``MolType`` instance.

.. jupyter-execute::

    seq = dna.make_seq("AGGCTT", name="seq1")
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

.. warning:: At present, constructing a custom ``MolType`` that overrides a builtin one affects the original (in this instance, the ``DnaSequence`` class). All subsequent calls to the original class in the running process that made the change are affected. The below code is resetting this attribute now to allow the rest of the documentation to be executed.

.. jupyter-execute::

    from cogent3 import DNA
    from cogent3.core.sequence import DnaSequence

    DnaSequence.moltype = DNA

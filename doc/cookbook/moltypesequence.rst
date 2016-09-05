**********************************************
Using the ``MolType`` and ``Sequence`` objects
**********************************************

.. authors Meg Pirrung

MolType
=======

``MolType`` provides services for resolving ambiguities, or providing the correct ambiguity for recoding. It also maintains the mappings between different kinds of alphabets, sequences and alignments.

One issue with ``MolType``'s is that they need to know about ``Sequence``, ``Alphabet``, and other objects, but, at the same time, those objects need to know about the ``MolType``. It is thus essential that the connection between these other types and the ``MolType`` can be made after the objects are created.

Setting up a ``MolType`` object with an RNA sequence
----------------------------------------------------

.. doctest::

   >>> from cogent3.core.moltype import MolType, IUPAC_RNA_chars,\
   ...   IUPAC_RNA_ambiguities, RnaStandardPairs, RnaMW,\
   ...   IUPAC_RNA_ambiguities_complements
   >>> from cogent3.core.sequence import NucleicAcidSequence
   >>> testrnaseq = 'ACGUACGUACGUACGU'
   >>> RnaMolType = MolType(
   ...     seq_constructor=NucleicAcidSequence(testrnaseq),
   ...     motifset=IUPAC_RNA_chars,
   ...     ambiguities=IUPAC_RNA_ambiguities,
   ...     label="rna_with_lowercase",
   ...     mw_calculator=RnaMW,
   ...     complements=IUPAC_RNA_ambiguities_complements,
   ...     pairs= RnaStandardPairs,
   ...     add_lower=True,
   ...     preserve_existing_moltypes=True,
   ...     make_alphabet_group=True,
   ...     )

Setting up a ``MolType`` object with a DNA sequence
---------------------------------------------------

.. doctest::

    >>> from cogent3.core.moltype import MolType, IUPAC_DNA_chars,\
    ...   IUPAC_DNA_ambiguities, DnaMW, IUPAC_DNA_ambiguities_complements,\
    ...   DnaStandardPairs
    >>> testdnaseq = 'ACGTACGTACGUACGT'
    >>> DnaMolType = MolType(
    ...     seq_constructor=NucleicAcidSequence(testdnaseq),
    ...     motifset=IUPAC_DNA_chars,
    ...     ambiguities=IUPAC_DNA_ambiguities,
    ...     label="dna_with_lowercase",
    ...     mw_calculator=DnaMW,
    ...     complements=IUPAC_DNA_ambiguities_complements,
    ...     pairs=DnaStandardPairs,
    ...     add_lower=True,
    ...     preserve_existing_moltypes=True,
    ...     make_alphabet_group=True,
    ...     )


Setting up a DNA ``MolType`` object allowing ``.`` as gaps
----------------------------------------------------------

.. doctest::

   >>> from cogent3.core import moltype as mt
   >>> DNAgapped = mt.MolType(seq_constructor=mt.DnaSequence,
   ...                        motifset=mt.IUPAC_DNA_chars,
   ...                        ambiguities=mt.IUPAC_DNA_ambiguities,
   ...                        complements=mt.IUPAC_DNA_ambiguities_complements,
   ...                        pairs=mt.DnaStandardPairs,
   ...                        gaps='.')
   >>> seq = DNAgapped.make_seq('ACG.')


Setting up a ``MolType`` object with a protein sequence
-------------------------------------------------------

.. doctest::

    >>> from cogent3.core.moltype import MolType, IUPAC_PROTEIN_chars,\
    ...   IUPAC_PROTEIN_ambiguities, ProteinMW
    >>> from cogent3.core.sequence import ProteinSequence, ArrayProteinSequence
    >>> protstr = 'TEST'
    >>> ProteinMolType = MolType(
    ...     seq_constructor=ProteinSequence(protstr),
    ...     motifset=IUPAC_PROTEIN_chars,
    ...     ambiguities=IUPAC_PROTEIN_ambiguities,
    ...     mw_calculator=ProteinMW,
    ...     make_alphabet_group=True,
    ...     array_seq_constructor=ArrayProteinSequence,
    ...     label="protein")
    >>> protseq = ProteinMolType.make_seq

Verify sequences
----------------

.. doctest::

   >>> rnastr = 'ACGUACGUACGUACGU'
   >>> dnastr = 'ACGTACGTACGTACGT'
   >>> RnaMolType.is_valid(rnastr)
   True
   >>> RnaMolType.is_valid(dnastr)
   False
   >>> RnaMolType.is_valid(NucleicAcidSequence(dnastr).to_rna())
   True

``Sequence``
============

The ``Sequence`` object contains classes that represent biological sequence data. These provide generic biological sequence manipulation functions, plus functions that are critical for the ``evolve`` module calculations.

.. warning:: Do not import sequence classes directly! It is expected that you will access them through ``MolType`` objects. The most common molecular types ``DNA``, ``RNA``, ``PROTEIN`` are provided as top level imports in cogent (e.g. ``cogent3.DNA``). Sequence classes depend on information from the ``MolType`` that is **only** available after ``MolType`` has been imported. Sequences are intended to be immutable. This is not enforced by the code for performance reasons, but don't alter the ``MolType`` or the sequence data after creation.

More detailed usage of sequence objects can be found in :ref:`dna-rna-seqs`.

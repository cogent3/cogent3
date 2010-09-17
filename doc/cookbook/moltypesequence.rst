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

   >>> from cogent.core.moltype import MolType, IUPAC_RNA_chars,\
   ...   IUPAC_RNA_ambiguities, RnaStandardPairs, RnaMW,\
   ...   IUPAC_RNA_ambiguities_complements
   >>> from cogent.core.sequence import NucleicAcidSequence
   >>> testrnaseq = 'ACGUACGUACGUACGU'
   >>> RnaMolType = MolType(
   ...     Sequence = NucleicAcidSequence(testrnaseq),
   ...     motifset = IUPAC_RNA_chars,
   ...     Ambiguities = IUPAC_RNA_ambiguities,
   ...     label = "rna_with_lowercase",
   ...     MWCalculator = RnaMW,
   ...     Complements = IUPAC_RNA_ambiguities_complements,
   ...     Pairs = RnaStandardPairs,
   ...     add_lower=True,
   ...     preserve_existing_moltypes=True,
   ...     make_alphabet_group=True,
   ...     )

Setting up a ``MolType`` object with a DNA sequence
---------------------------------------------------

.. doctest::

    >>> from cogent.core.moltype import MolType, IUPAC_DNA_chars,\
    ...   IUPAC_DNA_ambiguities, DnaMW, IUPAC_DNA_ambiguities_complements,\
    ...   DnaStandardPairs
   >>> testdnaseq = 'ACGTACGTACGUACGT'
   >>> DnaMolType = MolType(
   ...     Sequence = NucleicAcidSequence(testdnaseq),
   ...     motifset = IUPAC_DNA_chars,
   ...     Ambiguities = IUPAC_DNA_ambiguities,
   ...     label = "dna_with_lowercase",
   ...     MWCalculator = DnaMW,
   ...     Complements = IUPAC_DNA_ambiguities_complements,
   ...     Pairs = DnaStandardPairs,
   ...     add_lower=True,
   ...     preserve_existing_moltypes=True,
   ...     make_alphabet_group=True,
   ...     )


Setting up a ``MolType`` object with a protein sequence
-------------------------------------------------------

.. doctest::

    >>> from cogent.core.moltype import MolType, IUPAC_PROTEIN_chars,\
    ...   IUPAC_PROTEIN_ambiguities, ProteinMW
   >>> from cogent.core.sequence import ProteinSequence, ModelProteinSequence
   >>> protstr = 'TEST'
   >>> ProteinMolType = MolType(
   ...     Sequence = ProteinSequence(protstr),
   ...     motifset = IUPAC_PROTEIN_chars,
   ...     Ambiguities = IUPAC_PROTEIN_ambiguities,
   ...     MWCalculator = ProteinMW,
   ...     make_alphabet_group=True,
   ...     ModelSeq = ModelProteinSequence,
   ...     label = "protein")
   >>> protseq = ProteinMolType.Sequence

Verify sequences
----------------

.. doctest::

   >>> rnastr = 'ACGUACGUACGUACGU'
   >>> dnastr = 'ACGTACGTACGTACGT'
   >>> RnaMolType.isValid(rnastr)
   True
   >>> RnaMolType.isValid(dnastr)
   False
   >>> RnaMolType.isValid(NucleicAcidSequence(dnastr).toRna())
   True

``Sequence``
============

The ``Sequence`` object contains classes that represent biological sequence data. These provide generic biological sequence manipulation functions, plus functions that are critical for the ``evolve`` module calculations.

.. warning:: Do not import sequence classes directly! It is expected that you will access them through ``MolType`` objects. The most common molecular types ``DNA``, ``RNA``, ``PROTEIN`` are provided as top level imports in cogent (e.g. ``cogent.DNA``). Sequence classes depend on information from the ``MolType`` that is **only** available after ``MolType`` has been imported. Sequences are intended to be immutable. This is not enforced by the code for performance reasons, but don't alter the ``MolType`` or the sequence data after creation.

More detailed usage of sequence objects can be found in :ref:`dna-rna-seqs`.
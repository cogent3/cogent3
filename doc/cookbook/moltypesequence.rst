**************************************
Using the MolType and Sequence objects
**************************************

.. Meg Pirrung

MolType
=======
MolType provides services for resolving ambiguities, or providing the
correct ambiguity for recoding. It also maintains the mappings between
different kinds of alphabets, sequences and alignments.

One issue with MolTypes is that they need to know about Sequence, Alphabet,
and other objects, but, at the same time, those objects need to know about
the MolType. It is thus essential that the connection between these other
types and the MolType can be made after the objects are created.

Setting up a MolType object with an RNA sequence
------------------------------------------------

.. doctest::

   >>> from cogent.core.moltype import *
   >>> from cogent.core.sequence import *
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

Setting up a MolType object with a DNA sequence
-----------------------------------------------

.. doctest::

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


Setting up a MolType object with a protein sequence
---------------------------------------------------

.. doctest::

   >>> from cogent.core.moltype import *
   >>> from cogent.core.sequence import *
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

Sequence
========
The sequence object contains classes that represent biological sequence data. These
provide generic biological sequence manipulation functions, plus functions
that are critical for the EVOLVE calculations.

.. Warning::
   Do not import sequence classes directly! It is expected that you will
   access them through the moltype module. Sequence classes depend on information
   from the MolType that is **only** available after MolType has been imported.

   Sequences are intended to be immutable. This is not enforced by the code for
   performance reasons, but don't alter the MolType or the sequence data after
   creation.

Convert an RNA sequence to DNA
---------------------------------

.. doctest::

   >>> rnaseq = NucleicAcidSequence('ACGUACGUACGUACGU')
   >>> print rnaseq.toDna()
   ACGTACGTACGTACGT

Convert a DNA sequence to RNA
-----------------------------

.. doctest::

   >>> dnaseq = NucleicAcidSequence('ACGTACGTACGTACGT')
   >>> print dnaseq.toRna()
   ACGUACGUACGUACGU

Translate DNA into protein
--------------------------

.. doctest::

   >>> s = NucleicAcidSequence('GCTTGGGAAAGTCAAATGGAA')
   >>> print s.getTranslation()
   AWESQME

Convert a sequence to FASTA format
----------------------------------

.. doctest::

   >>> rnaseq.toFasta()
   '>0\nACGUACGUACGUACGU'


Return a randomized version of the sequence
-------------------------------------------

.. doctest::

   >>> print rnaseq.shuffle()
   ACAACUGGCUCUGAUG



Remove gaps from a sequence
---------------------------

.. doctest::

   >>> s = Sequence('--AUUAUGCUAU-UAu--')
   >>> print s.degap()
   AUUAUGCUAUUAU


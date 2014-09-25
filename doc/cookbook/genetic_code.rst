Getting a genetic code
----------------------

The standard genetic code.

.. doctest::
    
    >>> from cogent.core.genetic_code import GeneticCodes
    >>> standard_code = GeneticCodes[1]

The vertebrate mt genetic code.

.. doctest::
    
    >>> from cogent.core.genetic_code import GeneticCodes
    >>> mt_gc = GeneticCodes[2]
    >>> print mt_gc.Name
    Vertebrate Mitochondrial

To see the key -> genetic code mapping, use a loop.

.. doctest::
    
    >>> for key, code in GeneticCodes.items():
    ...     print key, code.Name
    1 Standard Nuclear
    2 Vertebrate Mitochondrial
    3 Yeast Mitochondrial...

Translate DNA sequences
-----------------------

.. doctest::

    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Translate all six frames
------------------------

.. doctest::
    
    >>> from cogent import DNA
    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> seq = DNA.makeSequence('ATGCTAACATAAA')
    >>> translations = standard_code.sixframes(seq)
    >>> print translations
    ['MLT*', 'C*HK', 'ANI', 'FMLA', 'LC*H', 'YVS']

Find out how many stops in a frame
----------------------------------

.. doctest::
    
    >>> from cogent import DNA
    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> seq = DNA.makeSequence('ATGCTAACATAAA')
    >>> stops_frame1 = standard_code.getStopIndices(seq, start=0)
    >>> stops_frame1
    [9]
    >>> stop_index = stops_frame1[0]
    >>> seq[stop_index:stop_index+3]
    DnaSequence(TAA)


Translate a codon
-----------------

.. doctest::

    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code['TTT']
    'F'

or get the codons for a single amino acid

.. doctest::

    >>> standard_code['A']
    ['GCT', 'GCC', 'GCA', 'GCG']

Look up the amino acid corresponding to a single codon
------------------------------------------------------

.. doctest::

    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code['TTT']
    'F'

Or get all the codons for one amino acid
----------------------------------------

.. doctest::

    >>> standard_code['A']
    ['GCT', 'GCC', 'GCA', 'GCG']

For a group of amino acids
--------------------------

.. doctest::

    >>> targets = ['A','C']
    >>> codons = [standard_code[aa] for aa in targets]
    >>> codons
    [['GCT', 'GCC', 'GCA', 'GCG'], ['TGT', 'TGC']]
    >>> flat_list = sum(codons,[])
    >>> flat_list
    ['GCT', 'GCC', 'GCA', 'GCG', 'TGT', 'TGC']

Converting the ``CodonAlphabet`` to codon series
------------------------------------------------

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGTT")
    >>> sorted(my_seq.CodonAlphabet())
    ['AAA', 'AAC', 'AAG', 'AAT'...
    >>> len(my_seq.CodonAlphabet())
    61

Obtaining the codons from a ``DnaSequence`` object
--------------------------------------------------

Use the method ``getInMotifSize``

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('ATGCACTGGTAA','my_gene')
    >>> codons = my_seq.getInMotifSize(3)
    >>> print codons
    ['ATG', 'CAC', 'TGG', 'TAA']

You can't translate a sequence that contains a stop codon.

.. doctest::
    
    >>> pep = my_seq.getTranslation()
    Traceback (most recent call last):
    AlphabetError: TAA

Remove the stop codon first
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('ATGCACTGGTAA','my_gene')
    >>> seq = my_seq.withoutTerminalStopCodon()
    >>> pep = seq.getTranslation()
    >>> print pep.toFasta()
    >my_gene
    MHW
    >>> print type(pep)
    <class 'cogent.core.sequence.ProteinSequence'>

Or we can just grab the correct slice from the ``DnaSequence`` object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('CAAATGTATTAA','my_gene')
    >>> pep = my_seq[:-3].getTranslation().toFasta()
    >>> print pep
    >my_gene
    QMY


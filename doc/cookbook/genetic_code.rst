Translate DNA sequences
-----------------------

.. doctest::

    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown here :ref:`translation`.

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

Use the method ``getInMotifSize()``

.. doctest::

    >>> from cogent import LoadSeqs,DNA
    >>> from cogent.core.alphabet import AlphabetError
    >>> my_seq = DNA.makeSequence('ATGCACTGGTAA','my_gene')
    >>> codons = my_seq.getInMotifSize(3)
    >>> print codons
    ['ATG', 'CAC', 'TGG', 'TAA']
    >>> try:
    ...     pep = my_seq.getTranslation()
    ... except AlphabetError, e:
    ...     print 'AlphabetError', e
    ...
    AlphabetError TAA

Remove the stop codon first
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs,DNA
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

    >>> from cogent import LoadSeqs,DNA
    >>> from cogent.core.alphabet import AlphabetError
    >>> my_seq = DNA.makeSequence('ATGCACTGGTAA','my_gene')
    >>> pep = my_seq[:-3].getTranslation().toFasta()
    >>> print pep
    >my_gene
    MHW


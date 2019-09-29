Translate DNA sequences
-----------------------

.. doctest::

    >>> from cogent3 import get_code
    >>> standard_code = get_code(1)
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Translate all six frames
------------------------

.. doctest::

    >>> from cogent3 import make_seq, get_code
    >>> standard_code = get_code(1)
    >>> seq = make_seq('ATGCTAACATAAA', moltype="dna")
    >>> translations = standard_code.sixframes(seq)
    >>> print(translations)
    ['MLT*', 'C*HK', 'ANI', 'FMLA', 'LC*H', 'YVS']

Find out how many stops in a frame
----------------------------------

.. doctest::

    >>> from cogent3 import make_seq, get_code
    >>> standard_code = get_code(1)
    >>> seq = make_seq('ATGCTAACATAAA', moltype="dna")
    >>> stops_frame1 = standard_code.get_stop_indices(seq, start=0)
    >>> stops_frame1
    [9]
    >>> stop_index = stops_frame1[0]
    >>> seq[stop_index:stop_index+3]
    DnaSequence(TAA)


Translate a codon
-----------------

.. doctest::

    >>> from cogent3 import make_seq, get_code
    >>> standard_code = get_code(1)
    >>> standard_code['TTT']
    'F'

or get the codons for a single amino acid

.. doctest::

    >>> standard_code['A']
    ['GCT', 'GCC', 'GCA', 'GCG']

Look up the amino acid corresponding to a single codon
------------------------------------------------------

.. doctest::

    >>> from cogent3 import get_code
    >>> standard_code = get_code(1)
    >>> standard_code['TTT']
    'F'

Or get all the codons for one amino acid
----------------------------------------

.. doctest::

    >>> from cogent3 import get_code
    >>> standard_code = get_code(1)
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

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq("AGTACACTGGTT", moltype="dna")
    >>> sorted(my_seq.codon_alphabet())
    ['AAA', 'AAC', 'AAG', 'AAT'...
    >>> len(my_seq.codon_alphabet())
    61

Obtaining the codons from a ``DnaSequence`` object
--------------------------------------------------

Use the method ``get_in_motif_size``

.. doctest::

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    >>> codons = my_seq.get_in_motif_size(3)
    >>> print(codons)
    ['ATG', 'CAC', 'TGG', 'TAA']

You can't translate a sequence that contains a stop codon.

.. doctest::

    >>> pep = my_seq.get_translation()
    Traceback (most recent call last):
    cogent3.core.alphabet.AlphabetError: TAA

Remove the stop codon first
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq('ATGCACTGGTAA',name='my_gene', moltype="dna")
    >>> seq = my_seq.trim_stop_codon()
    >>> pep = seq.get_translation()
    >>> print(pep.to_fasta())
    >my_gene
    MHW
    <BLANKLINE>
    >>> print(type(pep))
    <class 'cogent3.core.sequence.ProteinSequence'>

Or we can just grab the correct slice from the ``DnaSequence`` object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_seq
    >>> my_seq = make_seq('CAAATGTATTAA',name='my_gene', moltype="dna")
    >>> pep = my_seq[:-3].get_translation()
    >>> print(pep.to_fasta())
    >my_gene
    QMY
    <BLANKLINE>

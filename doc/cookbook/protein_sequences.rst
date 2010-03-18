Protein sequences
-----------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Creating a ProteinSequence with a name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import PROTEIN
    >>> p = PROTEIN.makeSequence('THISISAPRQTEIN','myProtein')
    >>> type(p)
    <class 'cogent.core.sequence.ProteinSequence'>
    >>> str(p)
    'THISISAPRQTEIN'

Converting a DNA sequence string to protein sequence string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Loading protein sequences from a Phylip file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs, PROTEIN
    >>> seq = LoadSeqs('data/abglobin_aa.phylip', moltype=PROTEIN,
    ...              aligned=True)

Loading other formats, or collections of sequences is shown in :ref:`load-seqs`.
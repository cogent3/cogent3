Protein sequences
-----------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Creating a ProteinSequence with a name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import PROTEIN
    >>> p = PROTEIN.make_seq('THISISAPRQTEIN','myProtein')
    >>> type(p)
    <class 'cogent3.core.sequence.ProteinSequence'>
    >>> str(p)
    'THISISAPRQTEIN'

Converting a DNA sequence string to protein sequence string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3.core.genetic_code import DEFAULT as standard_code
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Loading protein sequences from a Phylip file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> seq = load_aligned_seqs('data/abglobin_aa.phylip', moltype="protein")

Loading other formats, or collections of sequences is shown in :ref:`load-seqs`.
Perform a coevolutionary analysis on biological sequence alignments
===================================================================

.. sectionauthor:: Greg Caporaso

This document describes how to perform a coevolutionary analysis on a ``DenseAlignment`` object. Coevolutionary analyses identify correlated substitution patterns between ``DenseAlignment`` positions (columns). Several coevolution detection methods are currently provided via the PyCogent coevolution module. ``DenseAlignment`` objects must always be used as input to these functions. 

Before using an alignment in a coevolutionary analysis, you should be confident in the alignment. Poorly aligned sequences can yield very misleading results. There can be no ambiguous residue/base codes (e.g., B/Z/X in protein alignments) -- while some of the algorithms could tolerate them (e.g. Mutual Information), others which rely on information such as background residue frequencies (e.g. Statistical Coupling Analysis) cannot handle them. Some recoded amino acid alphabets will also not handle ambiguous residues. The best strategy is just to exclude ambiguous codes all together. To test for invalid characters before starting an analysis you can do the following:

.. doctest::

   >>> from cogent import LoadSeqs, PROTEIN, DNA, RNA
   >>> from cogent.core.alignment import DenseAlignment
   >>> from cogent.evolve.coevolution import validate_alignment
   >>> aln = LoadSeqs(data={'1':'GAA','2':'CTA', '3':'CTC','4':'-TC'},moltype=PROTEIN,aligned=DenseAlignment)
   >>> validate_alignment(aln)

To run a coevolutionary analysis, first create a ``DenseAlignment``:

.. doctest::

   >>> from cogent import LoadSeqs, PROTEIN, DNA, RNA
   >>> from cogent.core.alignment import DenseAlignment
   >>> aln = LoadSeqs(data={'1':'AAA','2':'CTA', '3':'CTC','4':'-TC'},moltype=PROTEIN,aligned=DenseAlignment)

Perform a coevolutionary analysis on a pair of positions in the alignment using mutual information (``mi``):

.. doctest::
    
    >>> from cogent.evolve.coevolution import coevolve_pair_functions, coevolve_pair
    >>> coevolve_pair(coevolve_pair_functions['mi'],aln,pos1=1,pos2=2)
    0.31127...

Perform a coevolutionary analysis on a pair of positions in the alignment using statistical coupling analysis (``sca``):

.. doctest::
    
    >>> from cogent.evolve.coevolution import coevolve_pair_functions, coevolve_pair
    >>> coevolve_pair(coevolve_pair_functions['sca'],aln,pos1=1,pos2=2,cutoff=0.5)
    0.98053...

Perform a coevolutionary analysis on one position and all other positions in the alignment using mutual information (``mi``):

.. doctest::
    
    >>> from cogent.evolve.coevolution import coevolve_position_functions, coevolve_position
    >>> coevolve_position(coevolve_position_functions['mi'],aln,position=1)
    array([        NaN,  0.81127812,  0.31127812])

Perform a coevolutionary analysis on all pairs of positions in the alignment using mutual information (``mi``):

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> from cogent.evolve.coevolution import coevolve_alignment_functions, coevolve_alignment
    >>> coevolve_alignment(coevolve_alignment_functions['mi'],aln)
    array([[        NaN,         NaN,         NaN],
           [        NaN,  0.81127812,  0.31127812],
           [        NaN,  0.31127812, 1.        ]])

View the available algorithms for computing coevolution values:

.. doctest::
    
    >>> print coevolve_pair_functions.keys()
    ['mi', 'sca', 'an', 'gctmpca', 'rmi', 'nmi']

Perform an intermolecular coevolutionary analysis using mutual information (``mi``). Note that there are strict requirements on the sequence identifiers for intermolecular analyses, and some important considerations involved in preparing alignments for these analyses. See the coevolve_alignments docstring (i.e., ``help(coevolve_alignments)`` from the python interpreter) for information. Briefly, sequence identifiers are split on ``+`` symbols. The ids before the + must match perfectly between the two alignments as these are used to match the sequences between alignments. In the following example, these are common species names: human, chicken, echidna, and pig. The text after the ``+`` can be anything, and should probably be the original database identifiers of the sequences.

.. doctest::

    >>> from cogent.evolve.coevolution import coevolve_alignment_functions,\
    ...   coevolve_alignments
    >>> aln1 = LoadSeqs(data={'human+protein1':'AAA','pig+protein1':'CTA',
    ...  'chicken+protein1':'CTC','echidna+weird_db_identifier':'-TC'},
    ...   moltype=PROTEIN,aligned=DenseAlignment)
    >>> aln2 = LoadSeqs(data={'pig+protein2':'AAAY','chicken+protein2':'CTAY',
    ...  'echidna+protein2':'CTCF','human+protein2':'-TCF'},
    ...   moltype=PROTEIN,aligned=DenseAlignment)
    >>> coevolve_alignments(coevolve_alignment_functions['mi'],aln1,aln2)
    array([[        NaN,         NaN,         NaN],
           [        NaN,  0.12255625,  0.31127812],
           [        NaN,  0.31127812,  0.        ],
           [        NaN,  0.31127812,  0.        ]])

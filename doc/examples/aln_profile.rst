Creating and manipulating alignment profiles
============================================

.. sectionauthor:: Sandra Smit

This is an example of how to create a profile from an alignment and how to do particular tricks with it. First, import the necessary stuff.

.. doctest::

    >>> from cogent.core.profile import Profile
    >>> from cogent import LoadSeqs, RNA

Then load an example alignment of 20 phe-tRNA sequences which we will use to create the profile

.. doctest::

    >>> aln = LoadSeqs("data/trna_profile.fasta", moltype=RNA)

Examine the number of sequences in the alignment and the alignment length
    
.. doctest:: 

    >>> print len(aln.Seqs)
    20
    >>> print len(aln)
    77

Create a profile containing the counts of each base at each alignment position

.. doctest::
    
    >>> pf = aln.getPosFreqs()
    >>> print pf.prettyPrint(include_header=True, column_limit=6, col_sep='   ')
     U    C    A    G    -   B
     0    0    0   20    0   0
     0   12    0    8    0   0
     1   18    0    1    0   0
     7    9    0    4    0   0...

Normalize the positions to get the relative frequencies at each position

.. doctest::
    
    >>> pf.normalizePositions()
    >>> print pf.prettyPrint(include_header=True, column_limit=6, col_sep='   ')
         U        C        A        G        -        B
    0.0000   0.0000   0.0000   1.0000   0.0000   0.0000
    0.0000   0.6000   0.0000   0.4000   0.0000   0.0000
    0.0500   0.9000   0.0000   0.0500   0.0000   0.0000
    0.3500   0.4500   0.0000   0.2000   0.0000   0.0000...

Make sure the data in the profile is valid. The method isValid checks whether all rows add up to one and whether the profile has a valid Alphabet and CharacterOrder.

.. doctest::

    >>> print pf.isValid()
    True

A profile can be used to calculate consensus sequences from the alignment. To illustrate the different options for consensus calculation, let's examine the frequency data at the fifth position of the alignment (index=4)

.. doctest::

    >>> print '\n'.join(['%s: %.3f'%(c,f) for (c,f) in\
    ...     zip(pf.CharOrder, pf.dataAt(4)) if f!=0])
    U: 0.050
    C: 0.400
    A: 0.250
    G: 0.300

The easiest consensus calculation will simply take the most frequent character at each position.

.. doctest::

    >>> print pf.toConsensus(fully_degenerate=False)
    GCCCCGGUAGCUCAGU--GGUAGAGCAGGGGACUGAAAAUCCCCGUGUCGGCGGUUCGAUUCCGUCCCGGGGCACCA

You can also specify to use the degenerate character needed to cover all symbols occurring at a certain alignment position (fully_degenerate=True). At index 4 in the alignment U, C, A, and G occur, thus the fully degenerate symbol needed is 'N'. Alternatively, using the cutoff value, you can ask for the degenerate symbol needed to cover a certain frequency. At a cutoff of 0.8, we need both C, G, and A at index 4 to cover this value, which results in the degenerate character 'V'. For the lower cutoff of 0.6, C and G suffice, and thus the character in the consensus sequence is 'S'.

.. doctest::

    >>> pf.Alphabet=RNA
    >>> print pf.toConsensus(fully_degenerate=True)
    GSBBNNDUAGCUCAGH??GGKAGAGCRBNVGRYUGAARAYCBNVNKGUCVBBDGWUCRAWHCHSNBHNNNVSC?CHM
    >>> print pf.toConsensus(cutoff=0.8)
    GSCYVBRUAGCUCAGU??GGUAGAGCASVSGAYUGAAAAUCYBSRUGUCSSYGGUUCGAUUCCGBSYSBRGSCACCA
    >>> print pf.toConsensus(cutoff=0.6)
    GCCYSGRUAGCUCAGU??GGUAGAGCAGRGGACUGAAAAUCCYCGUGUCGGYGGUUCGAUUCCGYCYCKRGGCACCA

A profile could also function as the description of a certain motif. As an example, let's create a profile description for the T-pseudouridine-C-loop which starts at index 54 and ends at index 59 (based on the reference structure matching the alignment).

.. doctest::
    
    >>> loop_profile = Profile(pf.Data[54:60,:], Alphabet=RNA, CharOrder=pf.CharOrder)
    >>> print loop_profile.prettyPrint(include_header=True, column_limit=6, col_sep='   ')
         U        C        A        G        -        B
    0.9500   0.0000   0.0500   0.0000   0.0000   0.0000
    1.0000   0.0000   0.0000   0.0000   0.0000   0.0000
    0.0000   1.0000   0.0000   0.0000   0.0000   0.0000
    0.0000   0.0000   0.0500   0.9500   0.0000   0.0000
    0.0000   0.0000   1.0000   0.0000   0.0000   0.0000
    0.8500   0.0000   0.1500   0.0000   0.0000   0.0000

We can calculate how well this profile matches in a certain sequence (or profile) by using the score method. As an example we see where the loop profile best fits into the yeast phe-tRNA sequence. As expected, we find the best hit at index 54 (with a score of 5.75).

.. doctest::

    >>> yeast = RNA.Sequence(\
    ...     'GCGGAUUUAGCUCAGUU-GGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA')
    >>> scores = loop_profile.score(yeast)
    >>> print scores
    [ 2.8   0.9   0.85  0.15  2.05  2.    3.75  0.95  1.2   1.    2.9   2.75
      0.    0.05  1.    2.9   2.05  1.95  0.2   1.95  0.05  1.    0.    2.
      0.15  2.    1.2   1.95  0.9   0.05  1.15  2.15  2.05  1.15  2.8   0.1
      0.9   0.    2.05  2.05  2.95  1.    1.8   0.95  0.05  0.85  2.    2.8
      0.95  1.85  2.75  1.    0.95  1.15  5.75  1.    0.    0.15  3.05  2.15
      1.    1.2   2.15  1.9   0.95  0.    0.05  1.05  4.05  1.95  1.05  0.15]
    >>> print max(scores)
    5.75
    >>> print scores.argmax()
    54


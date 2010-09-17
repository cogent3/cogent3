********************************
Analysis of sequence composition
********************************

.. sectionauthor:: Jesse Zaneveld

PyCogent provides several tools for analyzing the composition of DNA, RNA, or
protein sequences.

Loading your sequence
=====================

Let us say that we wish to study the sequence composition of the *Y. pseudotuberculosis* PB1 DNA Polymerase III beta subunit.

First we input the sequence as a string.

.. doctest::

    >>> y_pseudo_seq = \
    ...        """ atgaaatttatcattgaacgtgagcatctgctaaaaccactgcaacaggtcagtagcccg
    ...        ctgggtggacgccctacgttgcctattttgggtaacttgttgctgcaagtcacggaaggc
    ...        tctttgcggctgaccggtaccgacttggagatggagatggtggcttgtgttgccttgtct
    ...        cagtcccatgagccgggtgctaccacagtacccgcacggaagttttttgatatctggcgt
    ...        ggtttacccgaaggggcggaaattacggtagcgttggatggtgatcgcctgctagtgcgc
    ...        tctggtcgcagccgtttctcgctgtctaccttgcctgcgattgacttccctaatctggat
    ...        gactggcagagtgaggttgaattcactttaccgcaggctacgttaaagcgtctgattgag
    ...        tccactcagttttcgatggcccatcaggatgtccgttattatttgaacggcatgctgttt
    ...        gagaccgaaggcgaagagttacgtactgtggcgaccgatgggcatcgcttggctgtatgc
    ...        tcaatgcctattggccagacgttaccctcacattcggtgatcgtgccgcgtaaaggtgtg
    ...        atggagctggttcggttgctggatggtggtgatacccccttgcggctgcaaattggcagt
    ...        aataatattcgtgctcatgtgggcgattttattttcacatctaagctggttgatggccgt
    ...        ttcccggattatcgccgcgtattgccgaagaatcctgataaaatgctggaagccggttgc
    ...        gatttactgaaacaggcattttcgcgtgcggcaattctgtcaaatgagaagttccgtggt
    ...        gttcggctctatgtcagccacaatcaactcaaaatcactgctaataatcctgaacaggaa
    ...        gaagcagaagagatcctcgatgttagctacgaggggacagaaatggagatcggtttcaac
    ...        gtcagctatgtgcttgatgtgctaaatgcactgaagtgcgaagatgtgcgcctgttattg
    ...        actgactctgtatccagtgtgcagattgaagacagcgccagccaagctgcagcctatgtc
    ...        gtcatgccaatgcgtttgtag"""

To check that our results are reasonable, we can also load a small example string.

.. doctest::

    >>> example_seq = "GCGTTT"

In order to calculate compositional statistics, we need to import one of the ``Usage`` objects from ``cogent.core.usage``, create an object from our string, and normalize the counts contained in the string into frequencies. ``Usage`` objects include ``BaseUsage``, ``PositionalBaseUsage``, ``CodonUsage``, and ``AminoAcidUsage``.

Let us start with the ``BaseUsage`` object. The first few steps will be the same for the other Usage objects, however (as we will see below).

GC content
==========

Total GC  content
-----------------

GC content is one commonly used compositional statistic. To calculate the total GC content of our gene, we will need to initiate and normalize a ``BaseUsage`` object.

.. doctest::
    
    >>> from cogent.core.usage import BaseUsage
    >>> example_bu = BaseUsage(example_seq)
    >>> # Print raw counts
    >>> print example_bu.content("GC")
    3.0
    >>> example_bu.normalize()
    >>> print example_bu.content("GC")
    0.5

We can now visually verify that the reported GC contents are correct, and use the same technique on our full sequence.

.. doctest::

    >>> y_pseudo_bu = BaseUsage(y_pseudo_seq)
    >>> # Print raw counts
    >>> y_pseudo_bu.content("GC")
    555.0
    >>> y_pseudo_bu.normalize()
    >>> print y_pseudo_bu.content("GC")
    0.50408719346

Positional GC content of Codons
-------------------------------

When analyzing protein coding genes, it is often useful to subdivide the GC content by codon position. In particular, the 3rd codon position ``CodonUsage`` objects allow us to calculate the GC content at each codon position.

First, let us calculate the GC content for the codons in the example sequence as follows.

.. doctest::

    >>> # Import CodonUsage object
    >>> from cogent.core.usage import CodonUsage
    >>> # Initiate & normalize CodonUsage object
    >>> example_seq_cu = CodonUsage(example_seq)
    >>> example_seq_cu.normalize() 
    >>> GC,P1,P2,P3 = example_seq_cu.positionalGC()

Here, GC is the overall GC content for the sequence, while P1, P2, and P3 are the GC content at the first, second, and third codon positions, respectively.

Printing the results for the example gives the following results.

.. doctest::
    
    >>> print "GC:", GC
    GC: 0.5
    >>> print "P1:", P1
    P1: 0.5
    >>> print "P2:", P2
    P2: 0.5
    >>> print "P3:", P3
    P3: 0.5

We can then do the same for our biological sequence.
    
.. doctest::

    >>> y_pseudo_cu = CodonUsage(y_pseudo_seq)
    >>> y_pseudo_cu.normalize()
    >>> y_pseudo_GC = y_pseudo_cu.positionalGC()
    >>> print y_pseudo_GC
    [0.51874999999999993, 0.58437499999999987, 0.47500000000000009, 0.49687499999999996]

These results could then be fed into downstream analyses.

One important note is that ``CodonUsage`` objects calculate the GC content of codons within nucleotide sequences, rather than the full GC content.  Therefore, ``BaseUsage`` rather than ``CodonUsage`` objects should be used for calculating the GC content of non-coding sequences.

Total Base Usage
================

A more detailed view of composition incorporates the relative counts or frequencies of all bases. We can calculate total base usage as follows.

.. doctest::
    
    >>> from cogent.core.usage import BaseUsage
    >>> example_bu = BaseUsage(example_seq)
    >>> # Print raw counts
    >>> for k in example_bu.RequiredKeys:
    ...    print k, example_bu[k]
    A 0.0
    C 1.0
    U 3.0
    G 2.0
    >>> example_bu.normalize()
    >>> for k in example_bu.RequiredKeys:
    ...    print k, example_bu[k]
    A 0.0
    C 0.166666666667
    U 0.5
    G 0.333333333333

Dinucleotide Content
====================

The ``DinucUsage`` object allows us to calculate Dinucleotide usage for our sequence.

Dinucleotide usage can be calculated using overlapping, non-overlapping, or '3-1' dinucleotides.

Given the sequence "AATTAAGCC", each method will count dinucleotide usage differently. Overlapping dinucleotide usage will count "AA", "AT", "TT", "TA", "AA", "AG", "GC", "CC". Non-overlapping dinucleotide usage will count "AA", "TT", "AA", "GC" 3-1 dinucleotide usage will count "TT", "AC".

Calculating the GC content at the third and first codon positions ("3-1" usage) is useful for some applications, such as gene transfer detection, because changes at these positions tend to produce the most conservative amino acid substitutions, and thus are thought to better reflect mutational (rather than selective) pressure.

Overlapping dinucleotide content
--------------------------------

To calculate overlapping dinucleotide usage for our *Y. pseudotuberculosis* PB1 sequence.

.. doctest::

    >>> from cogent.core.usage import DinucUsage
    >>> du  = DinucUsage(y_pseudo_seq, Overlapping=True)
    >>> du.normalize()

We can inspect individual dinucleotide usages and confirm that the results add to 100% as follows

.. doctest::
    
    >>> total = 0.0
    >>> for k in du.RequiredKeys:
    ...    print k, du[k]
    ...    total += du[k]
    UU 0.0757855822551
    UC 0.0517560073937
    UA 0.043438077634
    UG 0.103512014787
    CU 0.0619223659889
    CC 0.0517560073937
    CA 0.0517560073937
    CG 0.0573012939002
    AU 0.0674676524954
    AC 0.043438077634
    AA 0.0573012939002
    AG 0.054528650647
    GU 0.0711645101664
    GC 0.0794824399261
    GA 0.0674676524954
    GG 0.0619223659889
    >>> print "Total:",total
    Total: 1.0

Non-overlapping Dinucleotide Content
------------------------------------

To calculate non-overlapping dinucleotide usage we simply change the ``Overlapping`` parameter to ``False`` when initiating the ``DinucUsage`` object.

.. doctest::

    >>> from cogent.core.usage import DinucUsage
    >>> du_no  = DinucUsage(y_pseudo_seq, Overlapping=False)
    >>> du_no.normalize()
    >>> total = 0
    >>> for k in du_no.RequiredKeys:
    ...    print k, du_no[k]
    ...    total += du_no[k]
    UU 0.0733082706767
    UC 0.0507518796992
    UA 0.0375939849624
    UG 0.105263157895
    CU 0.0733082706767
    CC 0.046992481203
    CA 0.0394736842105
    CG 0.0601503759398
    AU 0.0751879699248
    AC 0.046992481203
    AA 0.062030075188
    AG 0.0545112781955
    GU 0.0601503759398
    GC 0.0845864661654
    GA 0.0676691729323
    GG 0.062030075188
    >>> print "Total:",total
    Total: 1.0

'3-1' Dinucleotide Content
--------------------------

To calculate dinucleotide usage considering only adjacent first and third codon positions, we set the Overlapping parameter to '3-1' when constructing our ``DinucUsage`` object

.. doctest::
    
    >>> from cogent.core.usage import DinucUsage
    >>> du_3_1  = DinucUsage(y_pseudo_seq, Overlapping='3-1')
    >>> du_3_1.normalize()
    >>> total = 0
    >>> for k in du_3_1.RequiredKeys:
    ...    print k, du_3_1[k]
    ...    total += du_3_1[k]
    UU 0.0720221606648
    UC 0.0664819944598
    UA 0.0360110803324
    UG 0.0914127423823
    CU 0.0387811634349
    CC 0.0415512465374
    CA 0.0554016620499
    CG 0.0554016620499
    AU 0.0498614958449
    AC 0.0470914127424
    AA 0.0664819944598
    AG 0.0747922437673
    GU 0.0886426592798
    GC 0.0886426592798
    GA 0.0609418282548
    GG 0.0664819944598
    >>> print "Total:",total
    Total: 1.0

Comparing dinucleotide usages
-----------------------------

Above, we noted that there are several ways to calculate dinucleotide usages on a single sequence, and that the choice of methods changes the reported frequencies somewhat. How could we quantify the effect this choice make on the result?

One way to test this is to calculate the Euclidean distance between the resulting frequencies. We can do this using the dinucleotide usage's

.. doctest::
    
    >>> du_vs_du_3_1_dist = du.distance(du_3_1)

As required of a true distance, the results are independent of the direction of the calculation.

.. doctest::
    
    >>> du_3_1_vs_du_dist = du_3_1.distance(du)
    >>> print du_3_1_vs_du_dist == du_vs_du_3_1_dist
    True

Caution regarding unnormalized distances  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that in this case we have already called ``du.normalize()`` on each ``DinucUsage`` object. You MUST call ``du.normalize()`` before calculating distances. Otherwise the distance calculated will be for the dinucleotide counts, rather than frequencies. Distances of counts can be non-zero even for sequences with identical dinucleotide usage, if those sequences are of different lengths.

k-words
-------

*To be written.*

Codon usage analyses
====================

In addition to allowing a more detailed examination of GC content in coding sequences, ``CodonUsage`` objects (as the name implies) let us examine the codon usage of our sequence.

.. doctest::
   
    >>> from cogent.core.usage import CodonUsage
    >>> y_pseudo_cu  = CodonUsage(y_pseudo_seq)
    >>> # Print raw counts
    >>> for k in y_pseudo_cu.RequiredKeys:
    ...    print k, y_pseudo_cu[k]
    UUU 8.0
    UUC 4.0
    UUA 5.0
    UUG 14.0
    UCU 4.0
    UCC 3.0
    UCA 5.0
    UCG 3.0
    UAU 8.0...

Note that before normalization the ``CodonUsage`` object holds raw counts of results. However, for most purposes, we will want frequencies, so we normalize the counts.

.. doctest::
    
    >>> y_pseudo_cu.normalize()
    >>> # Print normalized frequencies 
    >>> for k in y_pseudo_cu.RequiredKeys:
    ...    print k, y_pseudo_cu[k]
    UUU 0.0225988700565
    UUC 0.0112994350282
    UUA 0.0141242937853
    UUG 0.0395480225989
    UCU 0.0112994350282
    UCC 0.00847457627119
    UCA 0.0141242937853
    UCG 0.00847457627119
    UAU 0.0225988700565...

Relative Synonymous Codon Usage
-------------------------------

The RSCU or relative synonymous codon usage metric divides the frequency of each codon by the total frequency of all codons encoding the same amino acid.

.. doctest::
    
    >>> y_pseudo_cu.normalize()
    >>> y_pseudo_rscu = y_pseudo_cu.rscu()
    >>> # Print rscu frequencies 
    >>> for k in y_pseudo_rscu.keys():
    ...    print k, y_pseudo_rscu[k]
    ACC 0.263157894737
    GUC 0.238095238095
    ACA 0.210526315789
    ACG 0.263157894737
    AAC 0.4
    CCU 0.315789473684
    UGG 1.0
    AUC 0.266666666667
    GUA 0.190476190476...

PR2 bias
--------

*To be written*

Fingerprint analysis
--------------------

*To be written*

Amino Acid Usage
================

*To be written.*

Profiles
========

*To be written.*

Visualisation
=============

*To be written.*

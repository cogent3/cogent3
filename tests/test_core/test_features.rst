Complete version of manipulating sequence annotations
=====================================================

.. sectionauthor:: Peter Maxwell, Gavin Huttley

A Sequence with a couple of exons on it.

.. doctest::
    
    >>> from cogent import DNA
    >>> from cogent.core.annotation import Feature
    >>> s = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA",
    ... Name="Orig")
    >>> exon1 = s.addAnnotation(Feature, 'exon', 'fred', [(10,15)])
    >>> exon2 = s.addAnnotation(Feature, 'exon', 'trev', [(30,40)])

The corresponding sequence can be extracted either with slice notation or by asking the feature to do it, since the feature knows what sequence it belongs to.

.. doctest::
    
    >>> s[exon1]
    DnaSequence(CCCCC)
    >>> exon1.getSlice()
    DnaSequence(CCCCC)

Usually the only way to get a ``Feature`` object like ``exon1`` is to ask the sequence for it. There is one method for querying annotations by type and optionally by name:

.. doctest::
    
    >>> exons = s.getAnnotationsMatching('exon')
    >>> print exons
    [exon "fred" at [10:15]/48, exon "trev" at [30:40]/48]

To construct a pseudo-feature covering (or excluding) multiple features, use ``getRegionCoveringAll``:

.. doctest::
    
    >>> print s.getRegionCoveringAll(exons)
    region "exon" at [10:15, 30:40]/48
    >>> print s.getRegionCoveringAll(exons).getShadow()
    region "not exon" at [0:10, 15:30, 40:48]/48

eg: all the exon sequence:

.. doctest::
    
    >>> s.getRegionCoveringAll(exons).getSlice()
    DnaSequence(CCCCCTT... 15)

or with slice notation:
    
.. doctest::
    
    >>> s[exon1, exon2]
    DnaSequence(CCCCCTT... 15)

Though ``.getRegionCoveringAll`` also guarantees no overlaps within the result, slicing does not:

.. doctest::
    
    >>> print s.getRegionCoveringAll(exons+exons)
    region "exon" at [10:15, 30:40]/48
    >>> s[exon1, exon1, exon1, exon1, exon1]
    Traceback (most recent call last):
    ValueError: Uninvertable. Overlap: 10 < 15

You can use features, maps, slices or integers, but non-monotonic slices are not allowed:

.. doctest::
    
    >>> s[15:20, 5:16]
    Traceback (most recent call last):
    ValueError: Uninvertable. Overlap: 15 < 16

Features are themselves sliceable:

.. doctest::
    
    >>> exon1[0:3].getSlice()
    DnaSequence(CCC)

When sequences are concatenated they keep their (non-overlapping) annotations:
    
.. doctest::
    
    >>> c = s[exon1[4:]]+s
    >>> print len(c)
    49
    >>> for feat in  c.annotations:
    ...     print feat
    ...
    exon "fred" at [-4-, 0:1]/49
    exon "fred" at [11:16]/49
    exon "trev" at [31:41]/49

Since features know their parents you can't use a feature from one sequence to slice another:
    
.. doctest::
    
    >>> print c[exon1]
    Traceback (most recent call last):
    ValueError: Can't map exon "fred" at [10:15]/48 onto ...

Features are generally attached to the thing they annotate, but in those cases where a free-floating feature is created it can later be attached:

.. doctest::
    
    >>> len(s.annotations)
    2
    >>> region = s.getRegionCoveringAll(exons)
    >>> len(s.annotations)
    2
    >>> region.attach()
    >>> len(s.annotations)
    3
    >>> region.detach()
    >>> len(s.annotations)
    2

When dealing with sequences that can be reverse complemented (e.g. ``DnaSequence``) features are also reversed. Feature types like CDS, however, have strand specific meaning and thus they're preserved in that orientation. We create a sequence with a CDS that spans multiple exons, and show that after getting the reverse complement we have exactly the same result from getting the CDS annotation.

.. doctest::
    
    >>> plus = DNA.makeSequence("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA",
    ... Name="plus")
    >>> plus_cds = plus.addAnnotation(Feature, 'CDS', 'gene',
    ...                           [(2,6),(10,15),(25,35)])
    >>> print plus_cds.getSlice()
    GGGGCCCCCTTTTTTTTTT
    >>> minus = plus.rc()
    >>> minus_cds = minus.getAnnotationsMatching('CDS')[0]
    >>> print minus_cds.getSlice()
    GGGGCCCCCTTTTTTTTTT


Sequence features can be accessed via a containing ``Alignment``:

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data=[['x','-AAAAAAAAA'], ['y','TTTT--TTTT']])
    >>> print aln
    >x
    -AAAAAAAAA
    >y
    TTTT--TTTT
    <BLANKLINE>
    >>> exon = aln.getSeq('x').addAnnotation(Feature, 'exon', 'fred', [(3,8)])
    >>> aln_exons = aln.getAnnotationsFromSequence('x', 'exon')
    >>> aln_exons = aln.getAnnotationsFromAnySequence('exon')

But these will be returned as **alignment** features with locations in alignment coordinates.

.. doctest::
    
    >>> print exon
    exon "fred" at [3:8]/9
    >>> print aln_exons[0]
    exon "fred" at [4:9]/10
    >>> print aln_exons[0].getSlice()
    >x
    AAAAA
    >y
    --TTT
    <BLANKLINE>
    >>> aln_exons[0].attach()
    >>> len(aln.annotations)
    1

Similarly alignment features can be projected onto the aligned sequences, where they may end up falling across gaps:

.. doctest::
    
    >>> exons = aln.getProjectedAnnotations('y', 'exon') 
    >>> print exons 
    [exon "fred" at [-2-, 4:7]/8]
    >>> print aln.getSeq('y')[exons[0].map.withoutGaps()]
    TTT

We copy the annotations from another sequence.

.. doctest::
    
    >>> aln = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']])
    >>> s = DNA.makeSequence("AAAAAAAAA", Name="x")
    >>> exon = s.addAnnotation(Feature, 'exon', 'fred', [(3,8)])
    >>> exon = aln.getSeq('x').copyAnnotations(s)
    >>> aln_exons = list(aln.getAnnotationsFromSequence('x', 'exon'))
    >>> print aln_exons
    [exon "fred" at [4:9]/10]

We consider cases where there are terminal gaps.

.. doctest::
    
    >>> aln = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', '------TTTT']])
    >>> exon = aln.getSeq('x').addFeature('exon', 'fred', [(3,8)])
    >>> aln_exons = list(aln.getAnnotationsFromSequence('x', 'exon'))
    >>> print aln_exons
    [exon "fred" at [4:9]/10]
    >>> print aln_exons[0].getSlice()
    >x
    AAAAA
    >y
    --TTT
    <BLANKLINE>
    >>> aln = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--T---']])
    >>> exon = aln.getSeq('x').addFeature('exon', 'fred', [(3,8)])
    >>> aln_exons = list(aln.getAnnotationsFromSequence('x', 'exon'))
    >>> print aln_exons[0].getSlice()
    >x
    AAAAA
    >y
    --T--
    <BLANKLINE>

In this case, only those residues included within the feature are covered - note the omission of the T in ``y`` opposite the gap in ``x``.

.. doctest::
    
    >>> aln = LoadSeqs(data=[['x', 'C-CCCAAAAA'], ['y', '-T----TTTT']])
    >>> print aln
    >x
    C-CCCAAAAA
    >y
    -T----TTTT
    <BLANKLINE>
    >>> exon = aln.getSeq('x').addFeature('exon', 'ex1', [(0,4)])
    >>> print exon
    exon "ex1" at [0:4]/9
    >>> print exon.getSlice()
    CCCC
    >>> aln_exons = list(aln.getAnnotationsFromSequence('x', 'exon'))
    >>> print aln_exons
    [exon "ex1" at [0:1, 2:5]/10]
    >>> print aln_exons[0].getSlice()
    >x
    CCCC
    >y
    ----
    <BLANKLINE>


``Feature.asOneSpan()``, is applied to the exon that straddles the gap in ``x``. The result is we preserve that feature.

.. doctest::
    
    >>> print aln_exons[0].asOneSpan().getSlice()
    >x
    C-CCC
    >y
    -T---
    <BLANKLINE>

Features can provide their coordinates, useful for custom analyses.
    
.. doctest::
    
    >>> all_exons = aln.getRegionCoveringAll(aln_exons)
    >>> coords = all_exons.getCoordinates()
    >>> assert coords == [(0,1),(2,5)]

Annotated regions can be masked (observed sequence characters replaced by another), either through the sequence on which they reside or by projection from the alignment. Note that ``mask_char`` must be a valid character for the sequence ``MolType``. Either the features (multiple can be named), or their shadow, can be masked.

We create an alignment with a sequence that has two different annotation types.

.. doctest::
    
    >>> aln = LoadSeqs(data=[['x', 'C-CCCAAAAAGGGAA'], ['y', '-T----TTTTG-GTT']])
    >>> print aln
    >x
    C-CCCAAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>
    >>> exon = aln.getSeq('x').addFeature('exon', 'norwegian', [(0,4)])
    >>> print exon.getSlice()
    CCCC
    >>> repeat = aln.getSeq('x').addFeature('repeat', 'blue', [(9,12)])
    >>> print repeat.getSlice()
    GGG
    >>> repeat = aln.getSeq('y').addFeature('repeat', 'frog', [(5,7)])
    >>> print repeat.getSlice()
    GG

Each sequence should correctly mask either the single feature, it's shadow, or the multiple features, or shadow.

.. doctest::
    
    >>> print aln.getSeq('x').withMaskedAnnotations('exon', mask_char='?')
    ????AAAAAGGGAA
    >>> print aln.getSeq('x').withMaskedAnnotations('exon', mask_char='?',
    ...                                         shadow=True)
    CCCC??????????
    >>> print aln.getSeq('x').withMaskedAnnotations(['exon', 'repeat'],
    ...                                           mask_char='?')
    ????AAAAA???AA
    >>> print aln.getSeq('x').withMaskedAnnotations(['exon', 'repeat'],
    ...                                           mask_char='?', shadow=True)
    CCCC?????GGG??
    >>> print aln.getSeq('y').withMaskedAnnotations('exon', mask_char='?')
    TTTTTGGTT
    >>> print aln.getSeq('y').withMaskedAnnotations('repeat', mask_char='?')
    TTTTT??TT
    >>> print aln.getSeq('y').withMaskedAnnotations('repeat', mask_char='?',
    ...                                          shadow=True)
    ?????GG??

The same methods can be applied to annotated Alignment's.

.. doctest::
    
    >>> print aln.withMaskedAnnotations('exon', mask_char='?')
    >x
    ?-???AAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>
    >>> print aln.withMaskedAnnotations('exon', mask_char='?', shadow=True)
    >x
    C-CCC??????????
    >y
    -?----?????-???
    <BLANKLINE>
    >>> print aln.withMaskedAnnotations('repeat', mask_char='?')
    >x
    C-CCCAAAAA???AA
    >y
    -T----TTTT?-?TT
    <BLANKLINE>
    >>> print aln.withMaskedAnnotations('repeat', mask_char='?', shadow=True)
    >x
    ?-????????GGG??
    >y
    -?----????G-G??
    <BLANKLINE>
    >>> print aln.withMaskedAnnotations(['repeat', 'exon'], mask_char='?')
    >x
    ?-???AAAAA???AA
    >y
    -T----TTTT?-?TT
    <BLANKLINE>
    >>> print aln.withMaskedAnnotations(['repeat', 'exon'],shadow=True)
    >x
    C-CCC?????GGG??
    >y
    -?----????G-G??
    <BLANKLINE>

It shouldn't matter whether annotated coordinates are entered separately, or as a series.

.. doctest::
    
    >>> data = [['human', 'CGAAACGTTT'], ['mouse', 'CTAAACGTCG']]
    >>> as_series = LoadSeqs(data = data)
    >>> as_items = LoadSeqs(data = data)

We add annotations to the sequences as a series.

.. doctest::
    
    >>> as_series.getSeq('human').addFeature('cpgsite', 'cpg', [(0,2), (5,7)])
    cpgsite "cpg" at [0:2, 5:7]/10
    >>> as_series.getSeq('mouse').addFeature('cpgsite', 'cpg', [(5,7), (8,10)])
    cpgsite "cpg" at [5:7, 8:10]/10

We add the annotations to the sequences one segment at a time.

.. doctest::
    
    >>> as_items.getSeq('human').addFeature('cpgsite', 'cpg', [(0,2)])
    cpgsite "cpg" at [0:2]/10
    >>> as_items.getSeq('human').addFeature('cpgsite', 'cpg', [(5,7)])
    cpgsite "cpg" at [5:7]/10
    >>> as_items.getSeq('mouse').addFeature('cpgsite', 'cpg', [(5,7)])
    cpgsite "cpg" at [5:7]/10
    >>> as_items.getSeq('mouse').addFeature('cpgsite', 'cpg', [(8,10)])
    cpgsite "cpg" at [8:10]/10

These different constructions should generate the same output.

.. doctest::
    
    >>> serial = as_series.withMaskedAnnotations(['cpgsite'])
    >>> print serial
    >human
    ??AAA??TTT
    >mouse
    CTAAA??T??
    <BLANKLINE>
    >>> itemwise = as_items.withMaskedAnnotations(['cpgsite'])
    >>> print itemwise
    >human
    ??AAA??TTT
    >mouse
    CTAAA??T??
    <BLANKLINE>

Annotations should be correctly masked, whether the sequence has been reverse complemented or not. We use the plus/minus strand CDS containing sequences created above.

.. doctest::
    
    >>> print plus.withMaskedAnnotations("CDS")
    AA????AAAA?????AAAAAAAAAA??????????AAA
    >>> print minus.withMaskedAnnotations("CDS")
    TTT??????????TTTTTTTTTT?????TTTT????TT

.. todo::
    
    Not documented, Source features.

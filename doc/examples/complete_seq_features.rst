Complete version of manipulating sequence annotations
=====================================================

.. sectionauthor:: Peter Maxwell, Gavin Huttley

A Sequence with a couple of exons on it.

.. doctest::
    
    >>> from cogent3 import DNA
    >>> from cogent3.core.annotation import Feature
    >>> s = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA",
    ... name="Orig")
    >>> exon1 = s.add_annotation(Feature, 'exon', 'fred', [(10,15)])
    >>> exon2 = s.add_annotation(Feature, 'exon', 'trev', [(30,40)])

The corresponding sequence can be extracted either with slice notation or by asking the feature to do it, since the feature knows what sequence it belongs to.

.. doctest::
    
    >>> s[exon1]
    DnaSequence(CCCCC)
    >>> exon1.get_slice()
    DnaSequence(CCCCC)

Usually the only way to get a ``Feature`` object like ``exon1`` is to ask the sequence for it. There is one method for querying annotations by type and optionally by name:

.. doctest::
    
    >>> exons = s.get_annotations_matching('exon')
    >>> print(exons)
    [exon "fred" at [10:15]/48, exon "trev" at [30:40]/48]

If the sequence does not have a matching feature you get back an empty list, and slicing the sequence with that returns a sequence of length 0.

.. doctest::
    
    >>> dont_exist = s.get_annotations_matching('dont_exist')
    >>> dont_exist
    []
    >>> s[dont_exist]
    DnaSequence()

To construct a pseudo-feature covering (or excluding) multiple features, use ``get_region_covering_all``:

.. doctest::
    
    >>> print(s.get_region_covering_all(exons))
    region "exon" at [10:15, 30:40]/48
    >>> print(s.get_region_covering_all(exons).get_shadow())
    region "not exon" at [0:10, 15:30, 40:48]/48

eg: all the exon sequence:

.. doctest::
    
    >>> s.get_region_covering_all(exons).get_slice()
    DnaSequence(CCCCCTT... 15)

or with slice notation:
    
.. doctest::
    
    >>> s[exon1, exon2]
    DnaSequence(CCCCCTT... 15)

Though ``.get_region_covering_all`` also guarantees no overlaps within the result, slicing does not:

.. doctest::
    
    >>> print(s.get_region_covering_all(exons+exons))
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
    
    >>> exon1[0:3].get_slice()
    DnaSequence(CCC)

When sequences are concatenated they keep their (non-overlapping) annotations:
    
.. doctest::
    
    >>> c = s[exon1[4:]]+s
    >>> print(len(c))
    49
    >>> for feat in  c.annotations:
    ...     print(feat)
    ...
    exon "fred" at [-4-, 0:1]/49
    exon "fred" at [11:16]/49
    exon "trev" at [31:41]/49

Since features know their parents you can't use a feature from one sequence to slice another:
    
.. doctest::
    
    >>> print(c[exon1])
    Traceback (most recent call last):
    ValueError: Can't map exon "fred" at [10:15]/48 onto ...

Features are generally attached to the thing they annotate, but in those cases where a free-floating feature is created it can later be attached:

.. doctest::
    
    >>> len(s.annotations)
    2
    >>> region = s.get_region_covering_all(exons)
    >>> len(s.annotations)
    2
    >>> region.attach()
    >>> len(s.annotations)
    3
    >>> region.detach()
    >>> len(s.annotations)
    2

When dealing with sequences that can be reverse complemented (e.g. ``DnaSequence``) features are **not** reversed. Features are considered to have strand specific meaning (.e.g CDS, exons) and so stay on their original strands. We create a sequence with a CDS that spans multiple exons, and show that after getting the reverse complement we have exactly the same result from getting the CDS annotation.

.. doctest::
    
    >>> plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA",
    ... name="plus")
    >>> plus_cds = plus.add_annotation(Feature, 'CDS', 'gene',
    ...                           [(2,6),(10,15),(25,35)])
    >>> print(plus_cds.get_slice())
    GGGGCCCCCTTTTTTTTTT
    >>> minus = plus.rc()
    >>> minus_cds = minus.get_annotations_matching('CDS')[0]
    >>> print(minus_cds.get_slice())
    GGGGCCCCCTTTTTTTTTT


Sequence features can be accessed via a containing ``Alignment``:

.. doctest::
    
    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs([['x','-AAAAAAAAA'], ['y','TTTT--TTTT']], array_align=False)
    >>> print(aln)
    >x
    -AAAAAAAAA
    >y
    TTTT--TTTT
    <BLANKLINE>
    >>> exon = aln.get_seq('x').add_annotation(Feature, 'exon', 'fred', [(3,8)])
    >>> aln_exons = aln.get_annotations_from_seq('x', 'exon')
    >>> aln_exons = aln.get_annotations_from_any_seq('exon')

But these will be returned as **alignment** features with locations in alignment coordinates.

.. doctest::
    
    >>> print(exon)
    exon "fred" at [3:8]/9
    >>> print(aln_exons[0])
    exon "fred" at [4:9]/10
    >>> print(aln_exons[0].get_slice())
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
    
    >>> exons = aln.get_projected_annotations('y', 'exon') 
    >>> print(exons)
    [exon "fred" at [-2-, 4:7]/8]
    >>> print(aln.get_seq('y')[exons[0].map.without_gaps()])
    TTT

We copy the annotations from another sequence,

.. doctest::
    
    >>> aln = make_aligned_seqs([['x', '-AAAAAAAAA'], ['y', 'TTTT--CCCC']], array_align=False)
    >>> s = DNA.make_seq("AAAAAAAAA", name="x")
    >>> exon = s.add_annotation(Feature, 'exon', 'fred', [(3,8)])
    >>> exon = aln.get_seq('x').copy_annotations(s)
    >>> aln_exons = list(aln.get_annotations_from_seq('x', 'exon'))
    >>> print(aln_exons)
    [exon "fred" at [4:9]/10]

even if the name is different.

.. doctest::
    
    >>> exon = aln.get_seq('y').copy_annotations(s)
    >>> aln_exons = list(aln.get_annotations_from_seq('y', 'exon'))
    >>> print(aln_exons)
    [exon "fred" at [3:4, 6:10]/10]
    >>> print(aln[aln_exons])
    >x
    AAAAA
    >y
    TCCCC
    <BLANKLINE>

If the feature lies outside the sequence being copied to, you get a lost span

.. doctest::

    >>> aln = make_aligned_seqs([['x', '-AAAA'], ['y', 'TTTTT']], array_align=False)
    >>> seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> exon = seq.add_feature('exon', 'A', [(5,8)])
    >>> aln.get_seq('x').copy_annotations(seq)
    >>> copied = list(aln.get_annotations_from_seq('x', 'exon'))
    >>> copied
    [exon "A" at [5:5, -4-]/5]
    >>> copied[0].get_slice()
    2 x 4 text alignment: x[----], y[----]

You can copy to a sequence with a different name, in a different alignment if the feature lies within the length

.. doctest::

    >>> aln = make_aligned_seqs([['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']], array_align=False)
    >>> seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.add_feature('exon', 'A', [(5,8)])
    >>> aln.get_seq('y').copy_annotations(seq)
    >>> copied = list(aln.get_annotations_from_seq('y', 'exon'))
    >>> copied
    [exon "A" at [7:10]/10]

If the sequence is shorter, again you get a lost span.

.. doctest::

    >>> aln = make_aligned_seqs([['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']], array_align=False)
    >>> diff_len_seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCCCCCCCCCC', 'x')
    >>> nonmatch = diff_len_seq.add_feature('repeat', 'A', [(12,14)])
    >>> aln.get_seq('y').copy_annotations(diff_len_seq)
    >>> copied = list(aln.get_annotations_from_seq('y', 'repeat'))
    >>> copied
    [repeat "A" at [10:10, -6-]/10]

We consider cases where there are terminal gaps.

.. doctest::
    
    >>> aln = make_aligned_seqs([['x', '-AAAAAAAAA'], ['y', '------TTTT']], array_align=False)
    >>> exon = aln.get_seq('x').add_feature('exon', 'fred', [(3,8)])
    >>> aln_exons = list(aln.get_annotations_from_seq('x', 'exon'))
    >>> print(aln_exons)
    [exon "fred" at [4:9]/10]
    >>> print(aln_exons[0].get_slice())
    >x
    AAAAA
    >y
    --TTT
    <BLANKLINE>
    >>> aln = make_aligned_seqs([['x', '-AAAAAAAAA'], ['y', 'TTTT--T---']], array_align=False)
    >>> exon = aln.get_seq('x').add_feature('exon', 'fred', [(3,8)])
    >>> aln_exons = list(aln.get_annotations_from_seq('x', 'exon'))
    >>> print(aln_exons[0].get_slice())
    >x
    AAAAA
    >y
    --T--
    <BLANKLINE>

In this case, only those residues included within the feature are covered - note the omission of the T in ``y`` opposite the gap in ``x``.

.. doctest::
    
    >>> aln = make_aligned_seqs([['x', 'C-CCCAAAAA'], ['y', '-T----TTTT']],
    ...                      moltype="dna", array_align=False)
    >>> print(aln)
    >x
    C-CCCAAAAA
    >y
    -T----TTTT
    <BLANKLINE>
    >>> exon = aln.get_seq('x').add_feature('exon', 'ex1', [(0,4)])
    >>> print(exon)
    exon "ex1" at [0:4]/9
    >>> print(exon.get_slice())
    CCCC
    >>> aln_exons = list(aln.get_annotations_from_seq('x', 'exon'))
    >>> print(aln_exons)
    [exon "ex1" at [0:1, 2:5]/10]
    >>> print(aln_exons[0].get_slice())
    >x
    CCCC
    >y
    ----
    <BLANKLINE>


``Feature.as_one_span()``, is applied to the exon that straddles the gap in ``x``. The result is we preserve that feature.

.. doctest::
    
    >>> print(aln_exons[0].as_one_span().get_slice())
    >x
    C-CCC
    >y
    -T---
    <BLANKLINE>

These properties also are consistently replicated with reverse complemented sequences.

.. doctest::
    
    >>> aln_rc = aln.rc()
    >>> rc_exons = list(aln_rc.get_annotations_from_any_seq('exon'))
    >>> print(aln_rc[rc_exons]) # not using as_one_span, so gap removed from x
    >x
    CCCC
    >y
    ----
    <BLANKLINE>
    >>> print(aln_rc[rc_exons[0].as_one_span()])
    >x
    C-CCC
    >y
    -T---
    <BLANKLINE>

Features can provide their coordinates, useful for custom analyses.
    
.. doctest::
    
    >>> all_exons = aln.get_region_covering_all(aln_exons)
    >>> coords = all_exons.get_coordinates()
    >>> assert coords == [(0,1),(2,5)]

Annotated regions can be masked (observed sequence characters replaced by another), either through the sequence on which they reside or by projection from the alignment. Note that ``mask_char`` must be a valid character for the sequence ``MolType``. Either the features (multiple can be named), or their shadow, can be masked.

We create an alignment with a sequence that has two different annotation types.

.. doctest::
    
    >>> aln = make_aligned_seqs([['x', 'C-CCCAAAAAGGGAA'], ['y', '-T----TTTTG-GTT']],
    ...               array_align=False)
    >>> print(aln)
    >x
    C-CCCAAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>
    >>> exon = aln.get_seq('x').add_feature('exon', 'norwegian', [(0,4)])
    >>> print(exon.get_slice())
    CCCC
    >>> repeat = aln.get_seq('x').add_feature('repeat', 'blue', [(9,12)])
    >>> print(repeat.get_slice())
    GGG
    >>> repeat = aln.get_seq('y').add_feature('repeat', 'frog', [(5,7)])
    >>> print(repeat.get_slice())
    GG

Each sequence should correctly mask either the single feature, it's shadow, or the multiple features, or shadow.

.. doctest::
    
    >>> print(aln.get_seq('x').with_masked_annotations('exon', mask_char='?'))
    ????AAAAAGGGAA
    >>> print(aln.get_seq('x').with_masked_annotations('exon', mask_char='?',
    ...                                         shadow=True))
    CCCC??????????
    >>> print(aln.get_seq('x').with_masked_annotations(['exon', 'repeat'],
    ...                                           mask_char='?'))
    ????AAAAA???AA
    >>> print(aln.get_seq('x').with_masked_annotations(['exon', 'repeat'],
    ...                                           mask_char='?', shadow=True))
    CCCC?????GGG??
    >>> print(aln.get_seq('y').with_masked_annotations('exon', mask_char='?'))
    TTTTTGGTT
    >>> print(aln.get_seq('y').with_masked_annotations('repeat', mask_char='?'))
    TTTTT??TT
    >>> print(aln.get_seq('y').with_masked_annotations('repeat', mask_char='?',
    ...                                          shadow=True))
    ?????GG??

The same methods can be applied to annotated Alignment's.

.. doctest::
    
    >>> print(aln.with_masked_annotations('exon', mask_char='?'))
    >x
    ?-???AAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>
    >>> print(aln.with_masked_annotations('exon', mask_char='?', shadow=True))
    >x
    C-CCC??????????
    >y
    -?----?????-???
    <BLANKLINE>
    >>> print(aln.with_masked_annotations('repeat', mask_char='?'))
    >x
    C-CCCAAAAA???AA
    >y
    -T----TTTT?-?TT
    <BLANKLINE>
    >>> print(aln.with_masked_annotations('repeat', mask_char='?', shadow=True))
    >x
    ?-????????GGG??
    >y
    -?----????G-G??
    <BLANKLINE>
    >>> print(aln.with_masked_annotations(['repeat', 'exon'], mask_char='?'))
    >x
    ?-???AAAAA???AA
    >y
    -T----TTTT?-?TT
    <BLANKLINE>
    >>> print(aln.with_masked_annotations(['repeat', 'exon'],shadow=True))
    >x
    C-CCC?????GGG??
    >y
    -?----????G-G??
    <BLANKLINE>

It shouldn't matter whether annotated coordinates are entered separately, or as a series.

.. doctest::
    
    >>> data = [['human', 'CGAAACGTTT'], ['mouse', 'CTAAACGTCG']]
    >>> as_series = make_aligned_seqs(data, array_align=False)
    >>> as_items = make_aligned_seqs(data, array_align=False)

We add annotations to the sequences as a series.

.. doctest::
    
    >>> as_series.get_seq('human').add_feature('cpgsite', 'cpg', [(0,2), (5,7)])
    cpgsite "cpg" at [0:2, 5:7]/10
    >>> as_series.get_seq('mouse').add_feature('cpgsite', 'cpg', [(5,7), (8,10)])
    cpgsite "cpg" at [5:7, 8:10]/10

We add the annotations to the sequences one segment at a time.

.. doctest::
    
    >>> as_items.get_seq('human').add_feature('cpgsite', 'cpg', [(0,2)])
    cpgsite "cpg" at [0:2]/10
    >>> as_items.get_seq('human').add_feature('cpgsite', 'cpg', [(5,7)])
    cpgsite "cpg" at [5:7]/10
    >>> as_items.get_seq('mouse').add_feature('cpgsite', 'cpg', [(5,7)])
    cpgsite "cpg" at [5:7]/10
    >>> as_items.get_seq('mouse').add_feature('cpgsite', 'cpg', [(8,10)])
    cpgsite "cpg" at [8:10]/10

These different constructions should generate the same output.

.. doctest::
    
    >>> serial = as_series.with_masked_annotations(['cpgsite'])
    >>> print(serial)
    >human
    ??AAA??TTT
    >mouse
    CTAAA??T??
    <BLANKLINE>
    >>> itemwise = as_items.with_masked_annotations(['cpgsite'])
    >>> print(itemwise)
    >human
    ??AAA??TTT
    >mouse
    CTAAA??T??
    <BLANKLINE>

Annotations should be correctly masked, whether the sequence has been reverse complemented or not. We use the plus/minus strand CDS containing sequences created above.

.. doctest::
    
    >>> print(plus.with_masked_annotations("CDS"))
    AA????AAAA?????AAAAAAAAAA??????????AAA
    >>> print(minus.with_masked_annotations("CDS"))
    TTT??????????TTTTTTTTTT?????TTTT????TT

.. todo::
    
    Not documented, Source features.

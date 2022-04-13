Annotating alignments and sequences
===================================

.. sectionauthor:: Peter Maxwell, Gavin Huttley

A Sequence with a couple of exons on it.

.. jupyter-execute::

    from cogent3 import DNA
    from cogent3.core.annotation import Feature

    s = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig")
    exon1 = s.add_annotation(Feature, "exon", "fred", [(10, 15)])
    exon2 = s.add_annotation(Feature, "exon", "trev", [(30, 40)])

The corresponding sequence can be extracted either with slice notation or by asking the feature to do it, since the feature knows what sequence it belongs to.

.. jupyter-execute::

    s[exon1]
    exon1.get_slice()

Usually the only way to get a ``Feature`` object like ``exon1`` is to ask the sequence for it. There is one method for querying annotations by type and optionally by name:

.. jupyter-execute::

    exons = s.get_annotations_matching("exon")
    print(exons)

If the sequence does not have a matching feature you get back an empty list, and slicing the sequence with that returns a sequence of length 0.

.. jupyter-execute::

    dont_exist = s.get_annotations_matching("dont_exist")
    dont_exist
    s[dont_exist]

To construct a pseudo-feature covering (or excluding) multiple features, use ``get_region_covering_all``:

.. jupyter-execute::

    print(s.get_region_covering_all(exons))
    print(s.get_region_covering_all(exons).get_shadow())

eg: all the exon sequence:

.. jupyter-execute::

    s.get_region_covering_all(exons).get_slice()

or with slice notation:
    
.. jupyter-execute::

    s[exon1, exon2]

Though ``.get_region_covering_all`` also guarantees no overlaps within the result, slicing does not:

.. jupyter-execute::
    :raises: ValueError

    print(s.get_region_covering_all(exons + exons))
    s[exon1, exon1, exon1, exon1, exon1]

You can use features, maps, slices or integers, but non-monotonic slices are not allowed:

.. jupyter-execute::
    :raises: ValueError

    s[15:20, 5:16]

Features are themselves sliceable:

.. jupyter-execute::

    exon1[0:3].get_slice()

When sequences are concatenated they keep their (non-overlapping) annotations:
    
.. jupyter-execute::

    c = s[exon1[4:]] + s
    print(len(c))
    for feat in c.annotations:
        print(feat)

Since features know their parents you can't use a feature from one sequence to slice another:
    
.. jupyter-execute::
    :raises: ValueError

    print(c[exon1])

Features are generally attached to the thing they annotate, but in those cases where a free-floating feature is created it can later be attached:

.. jupyter-execute::

    len(s.annotations)
    region = s.get_region_covering_all(exons)
    len(s.annotations)
    region.attach()
    len(s.annotations)
    region.detach()
    len(s.annotations)

When dealing with sequences that can be reverse complemented (e.g. ``DnaSequence``) features are **not** reversed. Features are considered to have strand specific meaning (.e.g CDS, exons) and so stay on their original strands. We create a sequence with a CDS that spans multiple exons, and show that after getting the reverse complement we have exactly the same result from getting the CDS annotation.

.. jupyter-execute::

    plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
    plus_cds = plus.add_annotation(Feature, "CDS", "gene", [(2, 6), (10, 15), (25, 35)])
    print(plus_cds.get_slice())
    minus = plus.rc()
    minus_cds = minus.get_annotations_matching("CDS")[0]
    print(minus_cds.get_slice())

Sequence features can be accessed via a containing ``Alignment``:

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    print(aln)
    exon = aln.get_seq("x").add_annotation(Feature, "exon", "fred", [(3, 8)])
    aln_exons = aln.get_annotations_from_seq("x", "exon")
    aln_exons = aln.get_annotations_from_any_seq("exon")

But these will be returned as **alignment** features with locations in alignment coordinates.

.. jupyter-execute::

    print(exon)
    print(aln_exons[0])
    print(aln_exons[0].get_slice())
    aln_exons[0].attach()
    len(aln.annotations)

Similarly alignment features can be projected onto the aligned sequences, where they may end up falling across gaps:

.. jupyter-execute::

    exons = aln.get_projected_annotations("y", "exon")
    print(exons)
    print(aln.get_seq("y")[exons[0].map.without_gaps()])

We copy the annotations from another sequence,

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--CCCC"]], array_align=False
    )
    s = DNA.make_seq("AAAAAAAAA", name="x")
    exon = s.add_annotation(Feature, "exon", "fred", [(3, 8)])
    exon = aln.get_seq("x").copy_annotations(s)
    aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
    print(aln_exons)

even if the name is different.

.. jupyter-execute::

    exon = aln.get_seq("y").copy_annotations(s)
    aln_exons = list(aln.get_annotations_from_seq("y", "exon"))
    print(aln_exons)
    print(aln[aln_exons])

If the feature lies outside the sequence being copied to, you get a lost span

.. jupyter-execute::

    aln = make_aligned_seqs([["x", "-AAAA"], ["y", "TTTTT"]], array_align=False)
    seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
    exon = seq.add_feature("exon", "A", [(5, 8)])
    aln.get_seq("x").copy_annotations(seq)
    copied = list(aln.get_annotations_from_seq("x", "exon"))
    copied
    copied[0].get_slice()

You can copy to a sequence with a different name, in a different alignment if the feature lies within the length

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
    match_exon = seq.add_feature("exon", "A", [(5, 8)])
    aln.get_seq("y").copy_annotations(seq)
    copied = list(aln.get_annotations_from_seq("y", "exon"))
    copied

If the sequence is shorter, again you get a lost span.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    diff_len_seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "x")
    nonmatch = diff_len_seq.add_feature("repeat", "A", [(12, 14)])
    aln.get_seq("y").copy_annotations(diff_len_seq)
    copied = list(aln.get_annotations_from_seq("y", "repeat"))
    copied

We consider cases where there are terminal gaps.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "------TTTT"]], array_align=False
    )
    exon = aln.get_seq("x").add_feature("exon", "fred", [(3, 8)])
    aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
    print(aln_exons)
    print(aln_exons[0].get_slice())
    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--T---"]], array_align=False
    )
    exon = aln.get_seq("x").add_feature("exon", "fred", [(3, 8)])
    aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
    print(aln_exons[0].get_slice())

In this case, only those residues included within the feature are covered - note the omission of the T in ``y`` opposite the gap in ``x``.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]], moltype="dna", array_align=False
    )
    print(aln)
    exon = aln.get_seq("x").add_feature("exon", "ex1", [(0, 4)])
    print(exon)
    print(exon.get_slice())
    aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
    print(aln_exons)
    print(aln_exons[0].get_slice())

``Feature.as_one_span()``, is applied to the exon that straddles the gap in ``x``. The result is we preserve that feature.

.. jupyter-execute::

    print(aln_exons[0].as_one_span().get_slice())

These properties also are consistently replicated with reverse complemented sequences.

.. jupyter-execute::

    aln_rc = aln.rc()
    rc_exons = list(aln_rc.get_annotations_from_any_seq("exon"))
    print(aln_rc[rc_exons])  # not using as_one_span, so gap removed from x
    print(aln_rc[rc_exons[0].as_one_span()])

Features can provide their coordinates, useful for custom analyses.
    
.. jupyter-execute::

    all_exons = aln.get_region_covering_all(aln_exons)
    coords = all_exons.get_coordinates()
    assert coords == [(0, 1), (2, 5)]

Annotated regions can be masked (observed sequence characters replaced by another), either through the sequence on which they reside or by projection from the alignment. Note that ``mask_char`` must be a valid character for the sequence ``MolType``. Either the features (multiple can be named), or their shadow, can be masked.

We create an alignment with a sequence that has two different annotation types.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "C-CCCAAAAAGGGAA"], ["y", "-T----TTTTG-GTT"]], array_align=False
    )
    print(aln)
    exon = aln.get_seq("x").add_feature("exon", "norwegian", [(0, 4)])
    print(exon.get_slice())
    repeat = aln.get_seq("x").add_feature("repeat", "blue", [(9, 12)])
    print(repeat.get_slice())
    repeat = aln.get_seq("y").add_feature("repeat", "frog", [(5, 7)])
    print(repeat.get_slice())

Each sequence should correctly mask either the single feature, it's shadow, or the multiple features, or shadow.

.. jupyter-execute::

    print(aln.get_seq("x").with_masked_annotations("exon", mask_char="?"))
    print(aln.get_seq("x").with_masked_annotations("exon", mask_char="?", shadow=True))
    print(aln.get_seq("x").with_masked_annotations(["exon", "repeat"], mask_char="?"))
    print(
        aln.get_seq("x").with_masked_annotations(
            ["exon", "repeat"], mask_char="?", shadow=True
        )
    )
    print(aln.get_seq("y").with_masked_annotations("exon", mask_char="?"))
    print(aln.get_seq("y").with_masked_annotations("repeat", mask_char="?"))
    print(
        aln.get_seq("y").with_masked_annotations("repeat", mask_char="?", shadow=True)
    )

The same methods can be applied to annotated Alignment's.

.. jupyter-execute::

    print(aln.with_masked_annotations("exon", mask_char="?"))
    print(aln.with_masked_annotations("exon", mask_char="?", shadow=True))
    print(aln.with_masked_annotations("repeat", mask_char="?"))
    print(aln.with_masked_annotations("repeat", mask_char="?", shadow=True))
    print(aln.with_masked_annotations(["repeat", "exon"], mask_char="?"))
    print(aln.with_masked_annotations(["repeat", "exon"], shadow=True))

It shouldn't matter whether annotated coordinates are entered separately, or as a series.

.. jupyter-execute::

    data = [["human", "CGAAACGTTT"], ["mouse", "CTAAACGTCG"]]
    as_series = make_aligned_seqs(data, array_align=False)
    as_items = make_aligned_seqs(data, array_align=False)

We add annotations to the sequences as a series.

.. jupyter-execute::

    as_series.get_seq("human").add_feature("cpgsite", "cpg", [(0, 2), (5, 7)])
    as_series.get_seq("mouse").add_feature("cpgsite", "cpg", [(5, 7), (8, 10)])

We add the annotations to the sequences one segment at a time.

.. jupyter-execute::

    as_items.get_seq("human").add_feature("cpgsite", "cpg", [(0, 2)])
    as_items.get_seq("human").add_feature("cpgsite", "cpg", [(5, 7)])
    as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(5, 7)])
    as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(8, 10)])

These different constructions should generate the same output.

.. jupyter-execute::

    serial = as_series.with_masked_annotations(["cpgsite"])
    print(serial)
    itemwise = as_items.with_masked_annotations(["cpgsite"])
    print(itemwise)

Annotations should be correctly masked, whether the sequence has been reverse complemented or not. We use the plus/minus strand CDS containing sequences created above.

.. jupyter-execute::

    print(plus.with_masked_annotations("CDS"))
    print(minus.with_masked_annotations("CDS"))

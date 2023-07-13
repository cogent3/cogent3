Annotating alignments and sequences
===================================

.. sectionauthor:: Peter Maxwell, Gavin Huttley

A Sequence with a couple of exons on it.

.. jupyter-execute::

    from cogent3 import DNA

    s = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig")
    exon1 = s.add_feature(biotype="exon", name="fred", spans=[(10, 15)])
    exon2 = s.add_feature(biotype="exon", name="trev", spans=[(30, 40)])

The corresponding sequence can be extracted either with slice notation or by asking the feature to do it, since the feature knows what sequence it belongs to.

.. jupyter-execute::

    s[exon1]

is the same as

.. jupyter-execute::

    exon1.get_slice()

For sequences annotated from external data (e.g. gff file), the way to get a ``Feature`` object like ``exon1`` is to ask the sequence for it by biotype and / or name:

.. jupyter-execute::

    exons = list(s.get_features(biotype="exon"))
    print(exons)

If the sequence does not have a matching feature you get back an empty list.

.. jupyter-execute::

    dont_exist = list(s.get_features(biotype="dont_exist"))
    dont_exist

To construct a pseudo-feature covering (or excluding) multiple features, use ``Feature.union()``:

.. jupyter-execute::

    all_exons = exon1.union([exon2])
    print(all_exons)
    print(all_exons.shadow())

eg: all the exon sequence:

.. jupyter-execute::

    all_exons.get_slice()

or with slice notation:

.. jupyter-execute::

    s[all_exons]

.. note:: You can use features, maps, slices or integers. Since features know their parents you can't use a feature from one sequence to slice another.

When dealing with sequences that can be reverse complemented (e.g. ``DnaSequence``) features are **not** reversed. Features are considered to have strand specific meaning (.e.g CDS, exons) and so stay on their original strands. We create a sequence with a CDS that spans multiple exons, and show that after getting the reverse complement we have exactly the same result from getting the CDS annotation.

.. jupyter-execute::

    plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
    plus_cds = plus.add_feature(
        biotype="CDS", name="gene", spans=[(2, 6), (10, 15), (25, 35)]
    )
    print(plus_cds.get_slice())
    minus = plus.rc()
    minus_cds = list(minus.get_features(biotype="CDS"))[0]
    print(minus_cds.get_slice())

Sequence features can be accessed via a containing ``Alignment``:

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False, moltype="dna"
    )
    print(aln)

But these will be returned as **alignment** features with locations in alignment coordinates.

.. jupyter-execute::

    exon = aln.get_seq("x").add_feature(biotype="exon", name="fred", spans=[(3, 8)])
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    exon

.. jupyter-execute::

    aln_exons[0]

.. jupyter-execute::

    aln_exons[0].get_slice()

Similarly alignment features can be projected onto the aligned sequences, where they may end up falling across gaps:

.. jupyter-execute::

    exons = list(aln.get_projected_features(seqid="y", biotype="exon"))
    print(exons)
    print(aln.get_seq("y")[exons[0].map.without_gaps()])

We copy the annotations from another sequence, so long as the name is the same.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--CCCC"]], array_align=False, moltype="dna"
    )
    s = DNA.make_seq("AAAAAAAAA", name="x")
    exon = s.add_feature(biotype="exon", name="fred", spans=[(3, 8)])
    aln.copy_annotations(s.annotation_db)
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    print(aln_exons)

To include a feature that lies partially outside the sequence, you get a lost span for the missing piece

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAA"], ["y", "TTTTT"]], array_align=False, moltype="dna"
    )
    exon = aln.add_feature(
        seqid="x", biotype="exon", name="A", spans=[(2, 8)], on_alignment=False
    )
    feature = list(aln.get_features(seqid="x", biotype="exon", allow_partial=True))[0]
    feature

but you can still get a slice.

.. jupyter-execute::

    feature.get_slice()

We consider cases where there are terminal gaps.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "------TTTT"]], array_align=False, moltype="dna"
    )
    exon = aln.add_feature(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    print(aln_exons)
    print(aln_exons[0].get_slice())
    aln = make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--T---"]], array_align=False, moltype="dna"
    )
    exon = aln.add_feature(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    print(aln_exons[0].get_slice())

In this case, only those residues included within the feature are covered - note the omission of the T in ``y`` opposite the gap in ``x``.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]], moltype="dna", array_align=False
    )
    print(aln)
    exon = aln.add_feature(seqid="x", biotype="exon", name="ex1", spans=[(0, 4)])
    print(exon)
    print(exon.get_slice())
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    print(aln_exons)
    print(aln_exons[0].get_slice())

``Feature.as_one_span()``, is applied to the exon that straddles the gap in ``x``. The result is we preserve that feature.

.. jupyter-execute::

    print(aln_exons[0].as_one_span().get_slice())

These properties also are consistently replicated with reverse complemented sequences.

Features can provide their coordinates, useful for custom analyses.

.. jupyter-execute::

    all_exons = aln_exons[0].union(aln_exons[1:])
    coords = all_exons.get_coordinates()
    assert coords == [(0, 1), (2, 5)]

Annotated regions can be masked (observed sequence characters replaced by another), either through the sequence on which they reside or by projection from the alignment. Note that ``mask_char`` must be a valid character for the sequence ``MolType``. Either the features (multiple can be named), or their shadow, can be masked.

We create an alignment with a sequence that has two different annotation types.

.. jupyter-execute::

    aln = make_aligned_seqs(
        [["x", "C-CCCAAAAAGGGAA"], ["y", "-T----TTTTG-GTT"]],
        array_align=False,
        moltype="dna",
    )
    print(aln)
    exon = aln.add_feature(seqid="x", biotype="exon", name="norwegian", spans=[(0, 4)])
    print(exon.get_slice())
    repeat = aln.add_feature(seqid="x", biotype="repeat", name="blue", spans=[(9, 12)])
    print(repeat.get_slice())
    repeat = aln.add_feature(seqid="y", biotype="repeat", name="frog", spans=[(5, 7)])
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

    as_series.get_seq("human").add_feature(
        biotype="cpgsite", name="cpg", spans=[(0, 2), (5, 7)]
    )
    as_series.get_seq("mouse").add_feature(
        biotype="cpgsite", name="cpg", spans=[(5, 7), (8, 10)]
    )

We add the annotations to the sequences one segment at a time.

.. jupyter-execute::

    as_items.get_seq("human").add_feature(biotype="cpgsite", name="cpg", spans=[(0, 2)])
    as_items.get_seq("human").add_feature(biotype="cpgsite", name="cpg", spans=[(5, 7)])
    as_items.get_seq("mouse").add_feature(biotype="cpgsite", name="cpg", spans=[(5, 7)])
    as_items.get_seq("mouse").add_feature(
        biotype="cpgsite", name="cpg", spans=[(8, 10)]
    )

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

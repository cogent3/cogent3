from unittest import TestCase, main

import pytest

from cogent3 import ASCII, DNA, make_aligned_seqs
from cogent3.core.annotation import Feature, Variable
from cogent3.core.annotation_db import GffAnnotationDb

# Complete version of manipulating sequence annotations
from cogent3.util.deserialise import deserialise_object


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class FeaturesTest(TestCase):
    """Tests of features in core"""

    def setUp(self):
        # A Sequence with a couple of exons on it.
        self.s = DNA.make_seq(
            "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig"
        )
        self.exon1 = self.s.add_feature(biotype="exon", name="fred", spans=[(10, 15)])
        self.exon2 = self.s.add_feature(biotype="exon", name="trev", spans=[(30, 40)])
        self.nested_feature = self.s.add_feature(
            biotype="repeat", name="bob", spans=[(12, 17)]
        )

    def test_exon_extraction(self):
        """exon feature used to slice or directly access sequence"""
        # The corresponding sequence can be extracted either with
        # slice notation or by asking the feature to do it,
        # since the feature knows what sequence it belongs to.

        self.assertEqual(str(self.s[self.exon1]), "CCCCC")
        self.assertEqual(str(self.exon1.get_slice()), "CCCCC")

    def test_get_features(self):
        """correctly identifies all features of a given type"""

        # Usually the only way to get a Feature object like exon1
        # is to ask the sequence for it. There is one method for querying
        # annotations by type and optionally by name:

        exons = list(self.s.get_features(biotype="exon"))
        self.assertEqual(
            str(exons),
            '["Orig" exon "fred" at [10:15]/48, "Orig" exon "trev" at [30:40]/48]',
        )

    def test_union(self):
        """combines multiple features into one"""

        # To construct a pseudo-feature covering (or excluding)
        # multiple features, use get_region_covering_all:

        exons = list(self.s.get_features(biotype="exon"))
        exon1 = exons.pop(0)
        combined = exon1.union(exons)
        self.assertEqual(str(combined.get_slice()), "CCCCCTTTTTAAAAA")

    def test_shadow(self):
        """combines multiple features into shadow"""

        # To construct a pseudo-feature covering (or excluding)
        # multiple features, use get_region_covering_all:

        exons = list(self.s.get_features(biotype="exon"))
        expect = str(
            self.s[: exons[0].map.start]
            + self.s[exons[0].map.end : exons[1].map.start]
            + self.s[exons[1].map.end :]
        )
        exon1 = exons.pop(0)
        shadow = exon1.union(exons).shadow()
        assert str(shadow.get_slice()) == expect

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_feature_reverse(self):
        """reverse complement of features"""

        # When dealing with sequences that can be reverse complemented
        # (e.g. DnaSequence) features are **not** reversed.
        # Features are considered to have strand specific meaning
        # (.e.g CDS, exons) and so stay on their original strands.
        # We create a sequence with a CDS that spans multiple exons,
        # and show that after getting the reverse complement we have
        # exactly the same result from getting the CDS annotation.

        plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        plus_cds = plus.add_annotation(
            Feature, "CDS", "gene", [(2, 6), (10, 15), (25, 35)]
        )
        self.assertEqual(str(plus_cds.get_slice()), "GGGGCCCCCTTTTTTTTTT")
        minus = plus.rc()
        minus_cds = minus.get_features(biotype="CDS")[0]
        self.assertEqual(str(minus_cds.get_slice()), "GGGGCCCCCTTTTTTTTTT")

    @pytest.mark.xfail(reason="todo gah implement lost spans")
    def test_lost_spans(self):
        """features no longer included in an alignment represented by lost spans"""

        # If the feature lies outside the sequence being copied to, you get a
        # lost span

        aln = make_aligned_seqs(
            data=[["x", "-AAAA"], ["y", "TTTTT"]], array_align=False
        )
        seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
        seq.add_feature("exon", "A", [(5, 8)])
        aln.get_seq("x").copy_annotations(seq)
        copied = list(aln.get_features(seqid="x", name="exon"))
        self.assertEqual(str(copied), '[exon "A" at [5:5, -4-]/5]')
        self.assertEqual(str(copied[0].get_slice()), ">x\n----\n>y\n----\n")

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_annotated_separately_equivalence(self):
        """allow defining features as a series or individually"""

        # It shouldn't matter whether annotated coordinates are entered
        # separately, or as a series.

        data = [["human", "CGAAACGTTT"], ["mouse", "CTAAACGTCG"]]
        as_series = make_aligned_seqs(data=data, array_align=False)
        as_items = make_aligned_seqs(data=data, array_align=False)

        # We add annotations to the sequences as a series.

        self.assertEqual(
            str(
                as_series.get_seq("human").add_feature(
                    biotype="cpgsite", name="cpg", spans=[(0, 2), (5, 7)]
                )
            ),
            'cpgsite "cpg" at [0:2, 5:7]/10',
        )
        self.assertEqual(
            str(
                as_series.get_seq("mouse").add_feature(
                    biotype="cpgsite", name="cpg", spans=[(5, 7), (8, 10)]
                )
            ),
            'cpgsite "cpg" at [5:7, 8:10]/10',
        )

        # We add the annotations to the sequences one segment at a time.

        self.assertEqual(
            str(as_items.get_seq("human").add_feature("cpgsite", "cpg", [(0, 2)])),
            'cpgsite "cpg" at [0:2]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("human").add_feature("cpgsite", "cpg", [(5, 7)])),
            'cpgsite "cpg" at [5:7]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(5, 7)])),
            'cpgsite "cpg" at [5:7]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(8, 10)])),
            'cpgsite "cpg" at [8:10]/10',
        )

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_constructor_equivalence(self):
        # These different constructions should generate the same output.
        data = [["human", "CGAAACGTTT"], ["mouse", "CTAAACGTCG"]]
        as_series = make_aligned_seqs(data=data, array_align=False)
        as_items = make_aligned_seqs(data=data, array_align=False)

        serial = as_series.with_masked_annotations(["cpgsite"])
        itemwise = as_items.with_masked_annotations(["cpgsite"])
        self.assertEqual(str(serial), str(itemwise))

        # Annotations should be correctly masked,
        # whether the sequence has been reverse complemented or not.
        # We use the plus/minus strand CDS containing sequences created above.
        plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        _ = plus.add_annotation(Feature, "CDS", "gene", [(2, 6), (10, 15), (25, 35)])
        minus = plus.rc()
        self.assertEqual(
            str(plus.with_masked_annotations("CDS")),
            "AA????AAAA?????AAAAAAAAAA??????????AAA",
        )
        self.assertEqual(
            str(minus.with_masked_annotations("CDS")),
            "TTT??????????TTTTTTTTTT?????TTTT????TT",
        )

    @pytest.mark.xfail(
        reason="todo gah implement support for annotation_db serialisation"
    )
    def test_roundtrip_json(self):
        """features can roundtrip from json"""

        seq = DNA.make_seq("AAAAATATTATTGGGT")
        seq.add_annotation(Feature, "exon", "myname", [(0, 5)])
        got = seq.to_json()
        new = deserialise_object(got)
        feat = new.get_features(biotype="exon")[0]
        self.assertEqual(str(feat.get_slice()), "AAAAA")

        # now with a list span
        seq = seq[3:]
        feat = seq.get_features(biotype="exon")[0]
        got = seq.to_json()
        new = deserialise_object(got)
        feat = new.get_features(biotype="exon")[0]
        self.assertEqual(str(feat.get_slice(complete=False)), "AA")

    @pytest.mark.xfail(
        reason="todo gah update test to use latest API plus support serialisation"
    )
    def test_roundtripped_alignment(self):
        """Alignment with annotations roundtrips correctly"""
        # annotations just on member sequences
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
        )
        _ = aln.get_seq("x").add_annotation(Feature, "exon", "fred", [(3, 8)])
        seq_exon = list(aln.get_features(seqid="x", name="exon"))[0]
        expect = seq_exon.get_slice()

        json = aln.to_json()
        new = deserialise_object(json)
        got_exons = list(new.get_features(seqid="x", name="exon"))[0]
        self.assertEqual(got_exons.get_slice().to_dict(), expect.to_dict())

        # annotations just on alignment
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAGGGG"], ["y", "TTTT--CCCC"]], array_align=False
        )
        f = aln.add_annotation(Feature, "generic", "no name", [(1, 4), (6, 10)])
        expect = f.get_slice().to_dict()
        json = aln.to_json()
        new = deserialise_object(json)
        got = list(new.get_features(biotype="generic"))[0]
        self.assertEqual(got.get_slice().to_dict(), expect)

        # annotations on both alignment and sequence
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAGGGG"], ["y", "TTTT--CCCC"]], array_align=False
        )
        f = aln.add_annotation(Feature, "generic", "no name", [(1, 4), (6, 10)])
        _ = aln.get_seq("x").add_annotation(Feature, "exon", "1", [(3, 8)])
        json = aln.to_json()
        new = deserialise_object(json)
        ## get back the exon
        seq_exon = list(aln.get_features(seqid="x", name="exon"))[0]
        expect = seq_exon.get_slice().to_dict()
        got_exons = list(new.get_features(seqid="x", name="exon"))[0]
        self.assertEqual(got_exons.get_slice().to_dict(), expect)
        ## get back the generic
        expect = f.get_slice().to_dict()
        got = list(new.get_features(biotype="generic"))[0]
        self.assertEqual(got.get_slice().to_dict(), expect)

        # check masking of seq features still works
        new = new.with_masked_annotations("exon", mask_char="?")
        self.assertEqual(new[4:9].to_dict(), dict(x="?????", y="--CCC"))

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_roundtripped_alignment2(self):
        """Sliced Alignment with annotations roundtrips correctly"""
        # annotations just on member sequences
        aln = make_aligned_seqs(
            data=[["x", "-AAAGGGGGAACCCT"], ["y", "TTTT--TTTTAGGGA"]], array_align=False
        )
        aln.get_seq("x").add_annotation(Feature, "exon", "E1", [(3, 8)])
        aln.get_seq("x").add_annotation(Feature, "exon", "E2", [(10, 13)])
        # at the alignment level
        sub_aln = aln[:-3]
        s = sub_aln.named_seqs["x"]
        d = s.data[:11]
        json = s.to_json()
        new = deserialise_object(json)
        gf1, gf2 = list(new.data.get_features(biotype="exon"))
        self.assertEqual(str(gf1.get_slice()), "GGGGG")
        self.assertEqual(str(gf2.get_slice()), "C")
        # the sliced alignment
        json = sub_aln.to_json()
        got = deserialise_object(json)
        x = got.named_seqs["x"]
        self.assertEqual(str(x.data.annotations[0].get_slice()), "GGGGG")
        self.assertEqual(str(x.data.annotations[1].get_slice()), "C")

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_roundtrip_rc_annotated_align(self):
        """should work for an alignment that has been reverse complemented"""
        # the key that exposed the bug was a gap in the middle of the sequence
        aln = make_aligned_seqs(
            data=[["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
            array_align=False,
            moltype="dna",
        )
        aln.get_seq("x").add_annotation(Feature, "exon", "E1", [(3, 8)])
        aln.get_seq("x").add_annotation(Feature, "exon", "E2", [(10, 13)])

        raln = aln.rc()
        json = raln.to_json()
        got = deserialise_object(json)
        self.assertEqual(got.to_dict(), raln.to_dict())
        orig_annots = {
            a.name: a.get_slice() for a in raln.get_annotations_from_any_seq()
        }
        got_annots = {a.name: a.get_slice() for a in got.get_annotations_from_any_seq()}
        self.assertEqual(got_annots, orig_annots)

    @pytest.mark.xfail(reason="todo gah delete test not supporting Variable class")
    def test_roundtrip_variable(self):
        """should recover the Variable feature type"""
        seq = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        xx_y = [[[2, 6], 2.4], [[10, 15], 5.1], [[25, 35], 1.3]]
        y_valued = seq.add_annotation(Variable, "SNP", "freq", xx_y)
        json = seq.to_json()
        new = deserialise_object(json)
        got = list(new.get_features(biotype="SNP"))[0]
        # annoyingly, comes back as list of lists
        self.assertEqual(got.xxy_list, [[list(xx), y] for xx, y in y_valued.xxy_list])

    @pytest.mark.xfail(reason="todo gah change test to use latest API")
    def test_nested_to_rich_dict(self):
        """check the to_rich_dict method works with nested annotations"""
        self.assertEqual(
            self.exon1.to_rich_dict()["annotations"][0],
            self.nested_feature.to_rich_dict(),
        )

    @pytest.mark.xfail(reason="todo gah update test to use latest API")
    def test_nested_deserialise_annotation(self):
        """nested annotations can be deserialised"""
        got = self.s.to_json()
        new = deserialise_object(got)
        new_exon1 = new.annotations[0]
        new_nested_feature = new_exon1.annotations[0]
        self.assertEqual(
            new_nested_feature.to_rich_dict(), self.nested_feature.to_rich_dict()
        )

    def test_annotate_matches_to(self):
        """annotate_matches_to attaches annotations correctly to a Sequence"""
        seq = DNA.make_seq("TTCCACTTCCGCTT", name="x")
        pattern = "CCRC"
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=True
        )
        self.assertEqual([a.get_slice() for a in annot], ["CCAC", "CCGC"])
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=False
        )
        self.assertEqual(len(annot), 1)
        fred = annot[0].get_slice()
        self.assertEqual(str(fred), "CCAC")
        # For Sequence objects of a non-IUPAC MolType, annotate_matches_to
        # should return an empty annotation.
        seq = ASCII.make_seq(seq="TTCCACTTCCGCTT")
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=False
        )
        self.assertEqual(annot, [])


def test_copy_annotations():
    """copying features from a db"""
    aln = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--CCCT"]], array_align=False
    )
    db = GffAnnotationDb()
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    aln.copy_annotations(db)
    feat = list(aln.get_features(seqid="y", biotype="exon"))[0]
    assert feat.get_slice().to_dict() == dict(x="AAA", y="CCT")


def test_feature_residue():
    """seq features on alignment operate in sequence coordinates"""
    # In this case, only those residues included within the feature are
    # covered - note the omission of the T in y opposite the gap in x.

    aln = make_aligned_seqs(
        data=[["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]],
        moltype=DNA,
        array_align=False,
    )
    assert str(aln), ">x\nC-CCCAAAAA\n>y\n-T----TTTT\n"
    exon = aln.get_seq("x").add_feature(biotype="exon", name="ex1", spans=[(0, 4)])
    assert 'exon "ex1" at [0:4]/9' in str(exon)
    assert str(exon.get_slice()), "CCCC"
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert 'exon "ex1" at [0:1, 2:5]/10' in str(aln_exons)
    assert aln_exons[0].get_slice().to_dict() == dict(x="CCCC", y="----")
    # Feature.as_one_span(), is applied to the exon that
    # straddles the gap in x. The result is we preserve that feature.
    exon_full_aln = aln_exons[0].as_one_span()
    assert exon_full_aln.get_slice().to_dict() == dict(x="C-CCC", y="-T---")

    # These properties also are consistently replicated with reverse
    # complemented sequences.

    aln_rc = aln.rc()
    rc_exons = list(aln_rc.get_features(biotype="exon"))[0]
    assert rc_exons.get_slice().to_dict() == dict(x="CCCC", y="----")
    assert rc_exons.as_one_span().get_slice().to_dict() == dict(x="C-CCC", y="-T---")


@pytest.fixture()
def ann_seq():
    # A Sequence with a couple of exons on it.
    s = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig")
    s.add_feature(biotype="gene", name="a-gene", spans=[(10, 40)])
    s.add_feature(
        biotype="exon", name="fred", spans=[(10, 15), (30, 40)], parent_id="a-gene"
    )
    return s


def test_get_features_no_matches(ann_seq):
    """get_features returns empty list if no matches"""

    # If the sequence does not have a matching feature
    # you get back an empty list, and slicing the sequence
    # with that returns a sequence of length 0.

    dont_exist = list(ann_seq.get_features(biotype="dont_exist"))
    assert dont_exist == []


def _add_features(obj, on_alignment):
    kwargs = dict(on_alignment=on_alignment) if on_alignment else {}
    obj.add_feature(biotype="CDS", name="GG", spans=[(0, 10)], strand="+", **kwargs)
    obj.add_feature(
        biotype="exon",
        name="child",
        spans=[(3, 6)],
        strand="+",
        parent_id="GG",
        **kwargs,
    )
    obj.add_feature(
        biotype="exon",
        name="not-child",
        spans=[(3, 6)],
        strand="+",
        parent_id="AA",
        **kwargs,
    )
    return obj


def test_feature_query_child_seq():
    s = DNA.make_seq("AAAGGGAAAA", name="s1")
    s = _add_features(s, on_alignment=False)
    gene = list(s.get_features(biotype="CDS"))[0]
    child = list(gene.get_children())
    assert len(child) == 1
    child = child[0]
    assert child.name == "child"
    assert str(child.get_slice()) == str(s[3:6])


def test_feature_query_parent_seq():
    s = DNA.make_seq("AAAGGGAAAA", name="s1")
    s = _add_features(s, on_alignment=False)
    exon = list(s.get_features(name="child"))[0]
    parent = list(exon.get_parent())
    assert len(parent) == 1
    parent = parent[0]
    assert parent.name == "GG"
    assert str(parent.get_slice()) == str(s[0:10])


def test_feature_query_child_aln():
    aln = make_aligned_seqs(
        data=[["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
        array_align=False,
        moltype="dna",
    )
    aln = _add_features(aln, on_alignment=True)
    gene = list(aln.get_features(biotype="CDS"))[0]
    child = list(gene.get_children())
    assert len(child) == 1
    child = child[0]
    assert child.name == "child"
    assert child.get_slice().to_dict() == aln[3:6].to_dict()


def test_feature_query_parent_aln():
    aln = make_aligned_seqs(
        data=[["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
        array_align=False,
        moltype="dna",
    )
    aln = _add_features(aln, on_alignment=True)
    child = list(aln.get_features(name="child"))[0]
    parent = list(child.get_parent())
    assert len(parent) == 1
    parent = parent[0]
    assert parent.name == "GG"
    assert parent.get_slice().to_dict() == aln[0:10].to_dict()


def test_aln_feature_lost_spans():
    """features outside the sequence should not be returned"""
    db = GffAnnotationDb(data=[])
    db.add_feature(seqid="y", biotype="repeat", name="A", spans=[(12, 14)])
    # If the sequence is shorter, again you get a lost span.
    aln = make_aligned_seqs(
        data={"x": "-AAAAAAAAA", "y": "TTTT--TTTT"}, array_align=False
    )
    aln.annotation_db = db
    copied = list(aln.get_features(seqid="y", biotype="repeat"))
    assert not copied


def test_terminal_gaps():
    """features in cases of terminal gaps"""

    # We consider cases where there are terminal gaps.
    db = GffAnnotationDb()
    feat = dict(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    db.add_feature(**feat)
    aln = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "------TTTT"]], array_align=False
    )
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert 'exon "fred" at [4:9]/10' in str(aln_exons)
    assert aln_exons[0].get_slice().to_dict() == dict(x="AAAAA", y="--TTT")
    aln = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--T---"]], array_align=False
    )
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert aln_exons[0].get_slice().to_dict() == dict(x="AAAAA", y="--T--")


def test_annotated_region_masks():
    """masking a sequence with specific features"""

    # Annotated regions can be masked (observed sequence characters
    # replaced by another), either through the sequence on which they
    # reside or by projection from the alignment. Note that mask_char must
    # be a valid character for the sequence MolType. Either the features
    # (multiple can be named), or their shadow, can be masked.

    # We create an alignment with a sequence that has two different annotation types.
    orig_data = {"x": "C-CCCAAAAAGGGAA", "y": "-T----TTTTG-GTT"}
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="exon", name="norwegian", spans=[(0, 4)])
    db.add_feature(
        biotype="repeat",
        name="blue",
        spans=[(9, 12)],
        seqid="x",
    )
    db.add_feature(seqid="y", biotype="repeat", name="frog", spans=[(5, 7)])
    aln = make_aligned_seqs(data=orig_data, array_align=False, moltype="dna")
    aln.annotation_db = db

    assert aln.to_dict() == {"x": "C-CCCAAAAAGGGAA", "y": "-T----TTTTG-GTT"}
    x = aln.get_seq("x")
    y = aln.get_seq("y")
    exon = list(x.get_features(biotype="exon"))[0]
    assert str(exon.get_slice()) == "CCCC"
    repeat_x = list(x.get_features(biotype="repeat"))[0]
    assert str(repeat_x.get_slice()) == "GGG"
    repeat_y = list(y.get_features(biotype="repeat"))[0]
    assert str(repeat_y.get_slice()) == "GG"

    # Each sequence should correctly mask either the single feature,
    # it's shadow, or the multiple features, or shadow.

    assert (
        str(aln.get_seq("x").with_masked_annotations("exon", mask_char="?"))
        == "????AAAAAGGGAA"
    )
    assert (
        str(
            aln.get_seq("x").with_masked_annotations("exon", mask_char="?", shadow=True)
        )
        == "CCCC??????????"
    )
    assert (
        str(aln.get_seq("x").with_masked_annotations(["exon", "repeat"], mask_char="?"))
        == "????AAAAA???AA"
    )
    assert (
        str(
            aln.get_seq("x").with_masked_annotations(
                ["exon", "repeat"], mask_char="?", shadow=True
            )
        )
        == "CCCC?????GGG??"
    )
    assert (
        str(aln.get_seq("y").with_masked_annotations("exon", mask_char="?"))
        == "TTTTTGGTT"
    )
    assert (
        str(aln.get_seq("y").with_masked_annotations("repeat", mask_char="?"))
        == "TTTTT??TT"
    )
    assert (
        str(
            aln.get_seq("y").with_masked_annotations(
                "repeat", mask_char="?", shadow=True
            )
        )
        == "?????GG??"
    )

    # The same methods can be applied to annotated Alignment's.

    assert aln.with_masked_annotations("exon", mask_char="?").to_dict() == {
        "x": "?-???AAAAAGGGAA",
        "y": "-T----TTTTG-GTT",
    }
    assert aln.with_masked_annotations(
        "exon", mask_char="?", shadow=True
    ).to_dict() == {"x": "C-CCC??????????", "y": "-?----?????-???"}
    assert aln.with_masked_annotations("repeat", mask_char="?").to_dict() == {
        "x": "C-CCCAAAAA???AA",
        "y": "-T----TTTT?-?TT",
    }
    assert aln.with_masked_annotations(
        "repeat", mask_char="?", shadow=True
    ).to_dict() == {"x": "?-????????GGG??", "y": "-?----????G-G??"}
    assert aln.with_masked_annotations(["repeat", "exon"], mask_char="?").to_dict() == {
        "x": "?-???AAAAA???AA",
        "y": "-T----TTTT?-?TT",
    }
    assert aln.with_masked_annotations(["repeat", "exon"], shadow=True).to_dict() == {
        "x": "C-CCC?????GGG??",
        "y": "-?----????G-G??",
    }


def test_nested_annotated_region_masks():
    """masking a sequence with specific features when nested annotations"""
    # gene = aln.get_seq("x").add_feature("gene", "norwegian", [(0, 4)])
    # gene.add_feature("repeat", "blue", [(1, 3)])
    # exon = aln.get_seq("y").add_feature("repeat", "frog", [(1, 4)])
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="gene", name="norwegian", spans=[(0, 4)])
    db.add_feature(seqid="x", biotype="repeat", name="blue", spans=[(1, 3)])
    db.add_feature(seqid="y", biotype="repeat", name="frog", spans=[(1, 4)])
    aln = make_aligned_seqs(
        data=[["x", "C-GGCAAAAATTTAA"], ["y", "-T----TTTTG-GTT"]], array_align=False
    )
    aln.annotation_db = db
    gene = list(aln.get_seq("x").get_features(biotype="gene"))[0]
    assert str(gene.get_slice()) == "CGGC"

    # evaluate the sequence directly
    masked = str(aln.get_seq("x").with_masked_annotations("repeat", mask_char="?"))
    assert masked == "C??CAAAAATTTAA"

    exon = list(aln.get_seq("y").get_features(biotype="repeat", name="frog"))[0]
    assert str(exon.get_slice()) == "TTT"
    # evaluate the sequence directly
    masked = str(aln.get_seq("y").with_masked_annotations("repeat", mask_char="?"))
    assert masked == "T???TGGTT"
    masked = aln.with_masked_annotations("gene", mask_char="?")
    got = masked.to_dict()
    assert got["x"] == "?-???AAAAATTTAA"
    assert got["y"] == "-T----TTTTG-GTT"


def test_feature_from_alignment():
    """seq features obtained from the alignment"""

    # we no longer support copying annotations individually
    # nor do we provide a mechanism for copying annotations from one
    # sequence to another

    # Sequence features can be accessed via a containing Alignment:
    db = GffAnnotationDb()
    aln = make_aligned_seqs(
        data={"x": "-AAAAAAAAA", "y": "TTTT--TTTT"}, array_align=False
    )
    db.add_feature(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert len(aln_exons) == 1
    aln_exons = aln_exons[0]

    # But these will be returned as **alignment**
    # features with locations in alignment coordinates.
    assert aln_exons.get_slice().to_dict() == {"x": "AAAAA", "y": "--TTT"}

    # Similarly alignment features can be projected onto the aligned sequences,
    # where they may end up falling across gaps:

    exons = aln.get_projected_annotations(seqid="y", biotype="exon")
    assert len(exons) == 1
    assert str(aln.get_seq("y")[exons[0].map.without_gaps()]), "TTT"
    assert 'exon "fred" at [-2-, 4:7]/8' in str(exons[0])


def test_nested_get_slice():
    """check the get_slice method works on nested annotations"""
    s = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig")
    ex = s.add_feature(biotype="exon", name="fred", spans=[(10, 20)])
    s.add_feature(biotype="exon", name="trev", spans=[(30, 40)])
    s.add_feature(biotype="repeat", name="bob", spans=[(12, 17)], parent_id="fred")
    f = list(ex.get_children())[0]
    assert str(f.get_slice()) == str(s[12:17])

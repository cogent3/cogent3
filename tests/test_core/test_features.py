from unittest import TestCase

import numpy
import pytest

import cogent3
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import BasicAnnotationDb, GffAnnotationDb

# Complete version of manipulating sequence annotations
from cogent3.util.deserialise import deserialise_object

ASCII = cogent3.get_moltype("text")
DNA = cogent3.get_moltype("dna")


class FeaturesTest(TestCase):
    """Tests of features in core"""

    def setUp(self):
        # A Sequence with a couple of exons on it.
        self.s = DNA.make_seq(
            seq="AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA",
            name="Orig",
        )
        self.exon1 = self.s.add_feature(biotype="exon", name="fred", spans=[(10, 15)])
        self.exon2 = self.s.add_feature(biotype="exon", name="trev", spans=[(30, 40)])
        self.nested_feature = self.s.add_feature(
            biotype="repeat",
            name="bob",
            spans=[(12, 17)],
        )

    def test_exon_extraction(self):
        """exon feature used to slice or directly access sequence"""
        # The corresponding sequence can be extracted either with
        # slice notation or by asking the feature to do it,
        # since the feature knows what sequence it belongs to.

        assert str(self.s[self.exon1]) == "CCCCC"
        assert str(self.exon1.get_slice()) == "CCCCC"

    def test_get_features(self):
        """correctly identifies all features of a given type"""

        # Usually the only way to get a Feature object like exon1
        # is to ask the sequence for it. There is one method for querying
        # annotations by type and optionally by name:

        exons = list(self.s.get_features(biotype="exon"))
        assert str(exons).startswith(
            "[Feature(seqid='Orig', biotype='exon', name='fred', map=[10:15]/48, parent=DnaSequence",
        )

    def test_union(self):
        """combines multiple features into one"""

        # To construct a pseudo-feature covering (or excluding)
        # multiple features, use get_region_covering_all:

        exons = list(self.s.get_features(biotype="exon"))
        exon1 = exons.pop(0)
        combined = exon1.union(exons)
        assert str(combined.get_slice()) == "CCCCCTTTTTAAAAA"

    def test_shadow(self):
        """combines multiple features into shadow"""

        # To construct a pseudo-feature covering (or excluding)
        # multiple features, use get_region_covering_all:

        exons = list(self.s.get_features(biotype="exon"))
        expect = str(
            self.s[: exons[0].map.start]
            + self.s[exons[0].map.end : exons[1].map.start]
            + self.s[exons[1].map.end :],
        )
        exon1 = exons.pop(0)
        shadow = exon1.union(exons).shadow()
        assert str(shadow.get_slice()) == expect

    def test_annotate_matches_to(self):
        """annotate_matches_to attaches annotations correctly to a Sequence"""
        seq = DNA.make_seq(seq="TTCCACTTCCGCTT", name="x")
        pattern = "CCRC"
        annot = seq.annotate_matches_to(
            pattern=pattern,
            biotype="domain",
            name="fred",
            allow_multiple=True,
        )
        assert [a.get_slice() for a in annot] == ["CCAC", "CCGC"]
        annot = seq.annotate_matches_to(
            pattern=pattern,
            biotype="domain",
            name="fred",
            allow_multiple=False,
        )
        assert len(annot) == 1
        fred = annot[0].get_slice()
        assert str(fred) == "CCAC"
        # For Sequence objects of a non-IUPAC MolType, annotate_matches_to
        # should return an empty annotation.
        seq = ASCII.make_seq(seq="TTCCACTTCCGCTT")
        annot = seq.annotate_matches_to(
            pattern=pattern,
            biotype="domain",
            name="fred",
            allow_multiple=False,
        )
        assert annot == []


def test_copy_annotations():
    """copying features from a db"""
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--CCCT"]],
        moltype="dna",
    )
    db = GffAnnotationDb()
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    aln.copy_annotations(db)
    feat = next(iter(aln.get_features(seqid="y", biotype="exon")))
    assert feat.get_slice().to_dict() == {"x": "AAA", "y": "CCT"}


def test_copy_annotations_onto_seq():
    """copying features onto a sequence"""
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--CCCT"]],
        moltype="dna",
    )
    db = BasicAnnotationDb()
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    y = aln.get_seq("y")
    y.copy_annotations(db)
    feat = next(iter(aln.get_features(seqid="y", biotype="exon")))
    assert feat.get_slice().to_dict() == {"x": "AAA", "y": "CCT"}


def test_feature_residue():
    """seq features on alignment operate in sequence coordinates"""
    # In this case, only those residues included within the feature are
    # covered - note the omission of the T in y opposite the gap in x.

    aln = cogent3.make_aligned_seqs(
        [["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]],
        moltype=DNA,
    )
    db = aln.annotation_db
    assert str(aln), ">x\nC-CCCAAAAA\n>y\n-T----TTTT\n"
    db.add_feature(seqid="x", biotype="exon", name="ex1", spans=[(0, 4)])
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert "biotype='exon', name='ex1', map=[0:1, 2:5]/10" in str(aln_exons)
    exon = aln_exons[0]
    exon_seq = exon.get_slice()
    assert exon_seq.to_dict() == {"x": "CCCC", "y": "----"}
    # Feature.as_one_span(), is applied to the exon that
    # straddles the gap in x. The result is we preserve that feature.
    exon_full_aln = aln_exons[0].as_one_span()
    assert exon_full_aln.get_slice().to_dict() == {"x": "C-CCC", "y": "-T---"}

    # These properties also are consistently replicated with reverse
    # complemented sequences.

    aln_rc = aln.rc()
    rc_exons = next(iter(aln_rc.get_features(biotype="exon")))
    assert rc_exons.get_slice().to_dict() == {"x": "CCCC", "y": "----"}
    assert rc_exons.as_one_span().get_slice().to_dict() == {"x": "C-CCC", "y": "-T---"}


@pytest.fixture
def ann_seq():
    # A Sequence with a couple of exons on it.
    s = DNA.make_seq(
        seq="AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA",
        name="Orig",
    )
    s.add_feature(biotype="gene", name="a-gene", spans=[(10, 40)])
    s.add_feature(
        biotype="exon",
        name="fred",
        spans=[(10, 15), (30, 40)],
        parent_id="a-gene",
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
    kwargs = {"on_alignment": on_alignment} if on_alignment else {}
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
    s = DNA.make_seq(seq="AAAGGGAAAA", name="s1")
    s = _add_features(s, on_alignment=False)
    gene = next(iter(s.get_features(biotype="CDS")))
    child = list(gene.get_children())
    assert len(child) == 1
    child = child[0]
    assert child.name == "child"
    assert str(child.get_slice()) == str(s[3:6])


def test_feature_query_parent_seq():
    s = DNA.make_seq(seq="AAAGGGAAAA", name="s1")
    s = _add_features(s, on_alignment=False)
    exon = next(iter(s.get_features(name="child")))
    parent = list(exon.get_parent())
    assert len(parent) == 1
    parent = parent[0]
    assert parent.name == "GG"
    assert str(parent.get_slice()) == str(s[0:10])


def test_feature_query_child_aln():
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
        moltype="dna",
    )
    aln = _add_features(aln, on_alignment=True)
    gene = next(iter(aln.get_features(biotype="CDS")))
    child = list(gene.get_children())
    assert len(child) == 1
    child = child[0]
    assert child.name == "child"
    assert child.get_slice().to_dict() == aln[3:6].to_dict()


def test_feature_query_parent_aln():
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
        moltype="dna",
    )
    aln = _add_features(aln, on_alignment=True)
    child = next(iter(aln.get_features(name="child")))
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
    aln = cogent3.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--TTTT"},
        moltype="dna",
    )
    aln.annotation_db = db
    copied = list(aln.get_features(seqid="y", biotype="repeat"))
    assert not copied


def test_terminal_gaps():
    """features in cases of terminal gaps"""

    # We consider cases where there are terminal gaps.
    db = GffAnnotationDb()
    feat = {"seqid": "x", "biotype": "exon", "name": "fred", "spans": [(3, 8)]}
    db.add_feature(**feat)
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "------TTTT"]],
        moltype="dna",
    )
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert "biotype='exon', name='fred', map=[4:9]/10" in str(aln_exons)
    assert aln_exons[0].get_slice().to_dict() == {"x": "AAAAA", "y": "--TTT"}
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAAAAAAA"], ["y", "TTTT--T---"]],
        moltype="dna",
    )
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert aln_exons[0].get_slice().to_dict() == {"x": "AAAAA", "y": "--T--"}


def test_feature_from_alignment():
    """seq features obtained from the alignment"""

    # we no longer support copying annotations individually
    # nor do we provide a mechanism for copying annotations from one
    # sequence to another

    # Sequence features can be accessed via a containing Alignment:
    db = GffAnnotationDb()
    aln = cogent3.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--TTTT"},
        moltype="dna",
    )
    db.add_feature(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    aln.annotation_db = db
    aln_exons = list(aln.get_features(seqid="x", biotype="exon"))
    assert len(aln_exons) == 1
    aln_exons = aln_exons[0]

    # But these will be returned as **alignment**
    # features with locations in alignment coordinates.
    assert aln[aln_exons].to_dict() == {"x": "AAAAA", "y": "--TTT"}

    # Similarly alignment features can be projected onto the aligned sequences,
    # where they may end up falling across gaps:

    exons = aln.get_projected_features(seqid="y", biotype="exon")
    assert len(exons) == 1
    assert str(aln.get_seq("y")[exons[0].map.without_gaps()]), "TTT"
    assert "biotype='exon', name='fred', map=[-2-, 4:7]/8" in str(exons[0])


def test_nested_get_slice():
    """check the get_slice method works on nested annotations"""
    s = DNA.make_seq(
        seq="AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA",
        name="Orig",
    )
    ex = s.add_feature(biotype="exon", name="fred", spans=[(10, 20)])
    s.add_feature(biotype="exon", name="trev", spans=[(30, 40)])
    s.add_feature(biotype="repeat", name="bob", spans=[(12, 17)], parent_id="fred")
    f = next(iter(ex.get_children()))
    assert str(s[f]) == str(s[12:17])


def test_masking_strand_agnostic_seq():
    db = GffAnnotationDb()
    db.add_feature(
        seqid="plus",
        biotype="CDS",
        name="gene",
        spans=[(2, 6), (10, 15), (25, 35)],
    )

    # Annotations should be correctly masked,
    # whether the sequence has been reverse complemented or not.
    # We use the plus/minus strand CDS containing sequences created above.
    plus = DNA.make_seq(seq="AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
    plus.annotation_db = db
    masked = plus.with_masked_annotations("CDS")
    assert len(masked) == len(plus)
    assert str(masked) == "AA????AAAA?????AAAAAAAAAA??????????AAA"
    minus = plus.rc()
    masked = minus.with_masked_annotations("CDS")
    assert len(masked) == len(minus)
    assert str(masked) == "TTT??????????TTTTTTTTTT?????TTTT????TT"


def test_masking_strand_agnostic_aln():
    db = GffAnnotationDb()
    db.add_feature(
        seqid="x",
        biotype="CDS",
        name="gene",
        spans=[(2, 6), (10, 15), (25, 35)],
    )
    aln = cogent3.make_aligned_seqs(
        {
            "x": "AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA",
            "y": "AAGGGGAAAACCCCCGGGGGGGGGGTTTTTTTTTTAAA",
        },
        moltype="dna",
    )
    aln.annotation_db = db
    masked = aln.with_masked_annotations("CDS")
    assert masked.to_dict() == {
        "x": "AA????AAAA?????AAAAAAAAAA??????????AAA",
        "y": str(aln.get_seq("y")),
    }
    rc = aln.rc()
    masked = rc.with_masked_annotations("CDS")
    assert masked.to_dict() == {
        "x": "TTT??????????TTTTTTTTTT?????TTTT????TT",
        "y": str(rc.get_seq("y")),
    }


def test_feature_out_range():
    """features no longer included in an alignment will not be returned"""
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAA"], ["y", "TTTTT"]],
        moltype="dna",
    )
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="exon", name="A", spans=[(5, 8)])
    f = list(aln.get_features(seqid="x", biotype="exon"))
    assert not f


@pytest.mark.parametrize("cast", [list, numpy.array])
def test_search_with_ints(cast):
    """searching for features with numpy ints should work"""
    start, stop = cast([2, 5])
    seq = DNA.make_seq(seq="AAAGGGGGAACCCT", name="x")
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="exon", name="E1", spans=[(3, 8)])
    db.add_feature(seqid="x", biotype="exon", name="E2", spans=[(10, 13)])
    seq.annotation_db = db
    feats = list(
        seq.get_features(biotype="exon", allow_partial=True, start=start, stop=stop),
    )
    assert len(feats) == 1


def test_roundtripped_alignment_with_slices():
    """Sliced Alignment with annotations roundtrips correctly"""
    # annotations just on member sequences
    aln = cogent3.make_aligned_seqs(
        [["x", "-AAAGGGGGAACCCT"], ["y", "TTTT--TTTTAGGGA"]],
        moltype="dna",
    )
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="exon", name="E1", spans=[(3, 8)])
    db.add_feature(seqid="x", biotype="exon", name="E2", spans=[(10, 13)])
    aln.annotation_db = db
    # at the alignment level
    sl = aln.seqs[0][:-3]
    assert str(sl) == "-AAAGGGGGAACCCT"[:-3]
    sub_aln = aln[:-3]
    feats = list(sub_aln.get_features(biotype="exon", allow_partial=True))
    assert len(feats) == 2
    # new type alignments DO NOT support serialising annotation_db's
    new = deserialise_object(sub_aln.to_json())
    feats = list(new.get_features(biotype="exon", allow_partial=True))
    assert not feats


def test_feature_reverse():
    """reverse complement of features"""

    # When dealing with sequences that can be reverse complemented
    # (e.g. DnaSequence) features are **not** reversed.
    # Features are considered to have strand specific meaning
    # (.e.g CDS, exons) and so stay on their original strands.
    # We create a sequence with a CDS that spans multiple exons,
    # and show that after getting the reverse complement we have
    # exactly the same result from getting the CDS annotation.

    plus = DNA.make_seq(
        seq="AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA",
        name="plus",
    )
    plus_cds = plus.add_feature(
        biotype="CDS",
        name="gene",
        spans=[(2, 6), (10, 15), (25, 35)],
    )
    assert str(plus_cds.get_slice()) == "GGGGCCCCCTTTTTTTTTT"
    minus = plus.rc()
    minus_cds = next(iter(minus.get_features(biotype="CDS")))
    assert str(minus_cds.get_slice()) == "GGGGCCCCCTTTTTTTTTT"


@pytest.mark.parametrize("moltype", ["protein", "bytes", "text"])
def test_rc_feature_on_wrong_moltype(moltype):
    moltype = cogent3.get_moltype(moltype)
    seq = moltype.make_seq(seq="AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="s1")
    cds = seq.add_feature(
        biotype="CDS",
        name="gene",
        spans=[(2, 6), (10, 15), (25, 35)],
        strand="-",
    )
    with pytest.raises((TypeError, AttributeError)):
        cds.get_slice()


def test_feature_equal(ann_seq):
    (f1,) = list(ann_seq.get_features(biotype="gene"))
    (f2,) = list(ann_seq.get_features(biotype="gene"))
    assert f1 is not f2
    assert f1 == f2


def test_feature_not_equal(ann_seq):
    (f1,) = list(ann_seq.get_features(biotype="gene"))
    (f2,) = list(ann_seq.get_features(biotype="exon"))
    assert f1 is not f2
    assert f1 != f2
    # same attributes except parent different parent seq
    nseq = ann_seq.copy()
    (nf1,) = list(nseq.get_features(biotype="gene"))
    assert nf1 != f1


@pytest.mark.parametrize("attr", ["seqid", "biotype", "name", "map"])
def test_feature_not_equal_attr(ann_seq, attr):
    (f1,) = list(ann_seq.get_features(biotype="gene"))
    attrs = {
        "parent": f1.parent,
        "seqid": f1.seqid,
        "biotype": f1.biotype,
        "map": f1.map,
        "name": f1.name,
        "strand": "-" if f1.reversed else "+",
    }
    value = attrs["map"][:4] if attr == "map" else "different"
    attrs[attr] = value
    f2 = Feature(**attrs)
    assert f1 != f2


def test_hash_feature(ann_seq):
    (f1,) = list(ann_seq.get_features(biotype="gene"))
    (f2,) = list(ann_seq.get_features(biotype="gene"))
    got = {f1, f2}
    assert len(got) == 1


def test_seq_degap_preserves_annotations():
    s1 = DNA.make_seq(
        seq="AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT",
        name="seq1",
    )
    s1.add_feature(biotype="exon", name="A", spans=[(10, 15)])
    dg = s1.degap()
    assert len(s1.annotation_db) == len(dg.annotation_db)


@pytest.mark.parametrize("aligned", [True, False])
def test_align_degap_preserves_annotations(aligned):
    data = {"seq1": "GATN--", "seq2": "?GATCT"}
    kwargs = {
        "moltype": DNA,
    }
    coll = (
        cogent3.make_aligned_seqs(data, **kwargs)
        if aligned
        else cogent3.make_unaligned_seqs(data, **kwargs)
    )
    db = BasicAnnotationDb()
    db.add_feature(biotype="exon", name="exon1", spans=[(1, 2)], seqid="seq1")
    coll.annotation_db = db
    got = coll.degap()
    assert got.annotation_db is coll.annotation_db
    assert len(got.annotation_db) == 1


@pytest.fixture
def feature_data():
    from cogent3.core.location import FeatureMap

    seq = DNA.make_seq(seq="ACGGTG", name="demo")
    fmap = FeatureMap.from_locations(locations=[(0, 2), (4, 6)], parent_length=6)
    return {
        "seqid": "1",
        "strand": 0,
        "map": fmap,
        "biotype": "repeat",
        "name": "trf",
        "parent": seq,
    }


@pytest.fixture
def feature_data_xattr(feature_data):
    feature_data["xattr"] = {
        "repeat_type": "Tandem repeats",
        "repeat_class": "trf",
        "repeat_name": "trf",
    }
    return feature_data


def test_feature_xattr(feature_data_xattr):
    feature = Feature(**feature_data_xattr)
    assert feature.xattr == feature_data_xattr["xattr"]


def test_feature_xttar_none(feature_data):
    feature = Feature(**feature_data)
    assert feature.xattr is None


def test_feature_xattr_set(feature_data):
    feature = Feature(**feature_data)
    with pytest.raises(TypeError):
        feature.xattr = {}

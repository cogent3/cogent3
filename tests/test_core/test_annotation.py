import unittest

import pytest

from cogent3 import (
    get_moltype,
    load_seq,
    make_aligned_seqs,
    make_unaligned_seqs,
)
from cogent3.core.location import FeatureMap, Span

DNA = get_moltype("dna")


def makeSampleSequence(name, with_gaps=False):  # ported
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    if with_gaps:
        raw_seq = f"{raw_seq[:5]}-----{raw_seq[10:-2]}--"
    # name is required for creating annotations
    seq = DNA.make_seq(seq=raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


def makeSampleAlignment():  # ported
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {seq1.name: seq1, seq2.name: seq2}
    aln = make_aligned_seqs(data=seqs, array_align=False)
    aln.add_feature(
        biotype="misc_feature",
        name="misc",
        spans=[(12, 25)],
        on_alignment=True,
    )
    aln.add_feature(biotype="CDS", name="blue", spans=[(15, 25)], on_alignment=True)
    aln.add_feature(biotype="5'UTR", name="red", spans=[(2, 4)], on_alignment=True)
    aln.add_feature(biotype="LTR", name="fake", spans=[(2, 15)], on_alignment=True)
    return aln


class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.seq = makeSampleSequence("seq1")
        self.aln = makeSampleAlignment()

    def test_add_annotated_seqs_drops_annotations(self):
        # retain link to annotation db as long as a simple slice
        a = self.seq[:5]
        b = self.seq[10:]
        assert a.annotation_db is not None
        assert b.annotation_db is not None
        # but adding two seqs drops the annotation db since too
        # difficult to track coords
        newseq = a + b
        assert newseq.annotation_db is None

    def test_aln_annotations(self):  # ported
        """test that annotations to alignment and its' sequences"""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
        }
        seq_expecteds = {
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "TTT", "FAKE02": "TTT"},
        }
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = (
                list(self.aln.get_features(biotype=annot_type, on_alignment=True))[0]
                .get_slice()
                .to_dict()
            )
            expected = aln_expecteds[annot_type]
            assert observed == expected, (annot_type, expected, observed)
            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it

            observed = list(
                self.aln.get_features(seqid=self.aln.names, biotype=annot_type),
            )[0]
            observed = observed.get_slice().to_dict()
            expected = seq_expecteds[annot_type]
            assert observed == expected


class TestMapSpans(unittest.TestCase):
    """Test attributes of Map & Spans classes critical to annotation
    manipulation."""

    def test_span(self):
        forward = Span(20, 30)
        reverse = Span(70, 80, reverse=True)
        assert forward.reversed_relative_to(100) == reverse
        assert reverse.reversed_relative_to(100) == forward


@pytest.mark.parametrize("alignment", (False, True))
def test_constructing_collections(alignment):  # ported
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {"FAKE01": seq1, "FAKE02": seq2}
    expect = sum(s.annotation_db.num_matches() for s in seqs.values())
    if alignment:
        coll = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    else:
        coll = make_unaligned_seqs(data=seqs, moltype="dna")

    assert coll.annotation_db.num_matches() == expect


@pytest.mark.parametrize("reversed", (False, True))
@pytest.mark.parametrize("annot_type", ("LTR", "misc_feature", "CDS", "5'UTR"))
def test_region_union_on_alignment(annot_type, reversed):  # ported
    # >FAKE01
    # AACCCAAAATTTTTTGGGGGGGGGGCCCC
    # >FAKE02
    # AACCC-----TTTTTGGGGGGGGGGCC--
    #   **  5'UTR
    #   ************* LTR
    #                **********  CDS
    #             *************  misc_feature

    aln = makeSampleAlignment()
    aln_expecteds = {
        "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
        "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
        "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
        "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
    }
    newaln = aln.rc() if reversed else aln
    feature_list = list(newaln.get_features(biotype=annot_type, on_alignment=True))
    new = feature_list[0].union(feature_list[1:])
    new = new.get_slice().to_dict()
    expected = aln_expecteds[annot_type]
    assert expected == new, (annot_type, expected, new)


@pytest.fixture
def ann_aln():
    # synthetic annotated alignment
    return makeSampleAlignment()


def test_feature_projection_ungapped(ann_aln):  # ported
    # projection onto ungapped sequence
    expecteds = {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"}
    aln_ltr = list(ann_aln.get_features(biotype="LTR"))[0]
    seq_name = "FAKE01"
    expected = expecteds[seq_name]
    num = ann_aln.annotation_db.num_matches()
    seq_ltr = ann_aln.get_projected_feature(seqid=seq_name, feature=aln_ltr)
    assert ann_aln.annotation_db.num_matches() == num + 1
    assert str(seq_ltr.get_slice()) == expected
    assert seq_ltr.seqid == seq_name
    assert seq_ltr.parent == ann_aln.get_seq(seq_name)


def test_feature_projection_gapped(ann_aln):  # ported
    # projection onto gapped sequence
    expecteds = {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"}
    aln_ltr = list(ann_aln.get_features(biotype="LTR"))[0]
    seq_name = "FAKE02"
    expected = expecteds[seq_name]
    seq_ltr = ann_aln.get_projected_feature(seqid=seq_name, feature=aln_ltr)

    with pytest.raises(ValueError):
        seq_ltr.get_slice(complete=True)

    # to get the annotation on the seq coord
    seq_ltr = seq_ltr.without_lost_spans()
    expected = expected.replace("-", "")
    assert str(seq_ltr.get_slice()) == expected
    assert seq_ltr.seqid == seq_name
    assert seq_ltr.parent == ann_aln.get_seq(seq_name)


@pytest.fixture
def ann_seq():
    return makeSampleSequence("seq1")


@pytest.mark.parametrize("annot_type", ("CDS", "5'UTR"))
def test_slice_seq_with_full_annotations(ann_seq, annot_type):  # ported
    # this slice contains both features intact
    newseq = ann_seq[10:]
    orig = list(ann_seq.get_features(biotype=annot_type))[0]
    new = list(newseq.get_features(biotype=annot_type))[0]
    assert orig.name == new.name
    assert len(orig) == len(new)
    assert str(newseq[new]) == str(ann_seq[orig]), annot_type


@pytest.mark.parametrize("annot_type,num", (("CDS", 0), ("5'UTR", 1)))
def test_slice_seq_with_partial_end(ann_seq, annot_type, num):  # ported
    # this slice contains both features intact
    newseq = ann_seq[:14]
    # only UTR is present
    new = list(newseq.get_features(biotype=annot_type, allow_partial=True))
    assert len(new) == num, annot_type
    if num:
        feat = new[0]
        # length of the feature is the same as the original
        assert len(feat) == len(list(ann_seq.get_features(biotype=annot_type))[0])
        gapless = feat.without_lost_spans()
        # the sliced feature without gaps is shorter
        assert len(gapless) < len(feat)


@pytest.mark.parametrize("annot_type,num", (("CDS", 1), ("5'UTR", 0)))  # ported
def test_slice_seq_with_partial_start(ann_seq, annot_type, num):
    # this slice contains both features intact
    newseq = ann_seq[18:]
    # only UTR is present
    new = list(newseq.get_features(biotype=annot_type, allow_partial=True))
    assert len(new) == num, annot_type
    if num:
        feat = new[0]
        # length of the feature is the same as the original
        assert len(feat) == len(list(ann_seq.get_features(biotype=annot_type))[0])
        gapless = feat.without_lost_spans()
        # the sliced feature without gaps is shorter
        assert len(gapless) < len(feat)


def test_seq_feature_to_dict():  # ported
    """create the attributes necessary to write into the user table"""
    seq = DNA.make_seq(seq="ATTGTACGCCCCTGA", name="test_seq")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": "+",
        "seqid": "test_seq",
        "on_alignment": False,
    }
    expect = {k: v for k, v in feature_data.items() if k != "on_alignment"}
    f = seq.make_feature(feature_data)
    got = f.to_dict()
    assert got == expect
    c = seq.add_feature(**got)
    assert c.to_dict() == expect
    assert str(f.get_slice()) == str(c.get_slice())


def test_aln_feature_to_dict():  # ported
    seqs = [
        makeSampleSequence("s1", with_gaps=False),
        makeSampleSequence("s2", with_gaps=True),
    ]
    aln = make_aligned_seqs(seqs, array_align=False)
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": "+",
        "seqid": None,
        "on_alignment": True,
    }
    # copy this now because modified by the method
    expect = {k: v for k, v in feature_data.items() if k != "on_alignment"}
    f = aln.make_feature(feature=feature_data)
    d = f.to_dict()
    assert d == expect


@pytest.mark.parametrize("fix", ("ann_seq", "ann_aln"))
@pytest.mark.parametrize("type_", (list, tuple))
def test_aln_slice_feat_invalid(type_, fix, request):  # ported
    # incorrect parent
    obj = request.getfixturevalue(fix)
    with pytest.raises(TypeError):
        _ = obj[type_(obj.get_features(biotype="exon"))]


def test_seq_slice_seqfeat_invalid(ann_aln):  # ported
    # incorrect parent
    seq1 = ann_aln.get_seq("FAKE01")
    seq2 = ann_aln.get_seq("FAKE02")
    with pytest.raises(ValueError):
        _ = seq2[list(seq1.get_features(biotype="CDS"))[0]]


def test_gbdb_get_children_get_parent(DATA_DIR):  # ported
    seq = load_seq(DATA_DIR / "annotated_seq.gb")
    seq = seq[2900:6000]
    (orig,) = list(seq.get_features(biotype="gene", name="CNA00110"))
    (child,) = list(orig.get_children("CDS"))
    parent, *_ = list(child.get_parent())
    assert parent == orig


@pytest.mark.parametrize("rev", (False, True))
def test_features_survives_seq_rename(rev):  # ported
    segments = ["A" * 10, "C" * 10, "T" * 5, "C" * 5, "A" * 5]

    seq = DNA.make_seq(seq="".join(segments), name="original")
    gene = seq.add_feature(biotype="gene", name="gene1", spans=[(10, 20), (25, 30)])
    gene_expect = str(seq[10:20]) + str(seq[25:30])
    assert str(gene.get_slice()) == gene_expect
    domain = seq.add_feature(
        biotype="domain",
        name="domain1",
        spans=[(20, 25)],
        strand="-",
    )
    domain_expect = str(seq[20:25].rc())
    domain_got = domain.get_slice()
    assert str(domain_got) == domain_expect
    sliced = seq[5:-3]
    sliced.name = "sliced"
    sliced = sliced.rc() if rev else sliced

    got = list(sliced.get_features(name="gene1"))[0]
    got = got.get_slice()
    assert str(got) == gene_expect

    got = list(sliced.get_features(name="domain1"))[0]
    got = got.get_slice()
    assert str(got) == domain_expect


def make_aligned(**kwargs):
    return make_aligned_seqs(array_align=False, **kwargs)


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("make_cls", (make_unaligned_seqs, make_aligned))
def test_features_survives_aligned_seq_rename(rev, make_cls):  # ported
    segments = ["A" * 10, "C" * 10, "T" * 5, "C" * 5, "A" * 5]

    seqs = make_cls(data={"original": "".join(segments)}, moltype="dna")
    seqs.annotation_db.add_feature(
        seqid="original",
        biotype="gene",
        name="gene1",
        spans=[(10, 20), (25, 30)],
    )
    seqs.annotation_db.add_feature(
        seqid="original",
        biotype="domain",
        name="domain1",
        spans=[(20, 25)],
        strand="-",
    )
    seqs = seqs.rename_seqs(lambda x: "newname")
    assert seqs.names == ["newname"]
    seqs = seqs.rc() if rev else seqs
    # quite different behaviour from Alignment and SequenceCollection
    # so we convert to string to make comparison simpler
    got = list(seqs.get_features(name="gene1"))[0]
    sliced = str(got.get_slice()).splitlines()[-1]
    assert sliced == "C" * 15


@pytest.mark.parametrize("make", (make_unaligned_seqs, make_aligned))
def test_features_invalid_seqid(make):  # ported
    segments = ["A" * 10, "C" * 10, "T" * 5, "C" * 5, "A" * 5]

    seqs = make(data={"original": "".join(segments)}, moltype="dna")
    seqs.annotation_db.add_feature(
        seqid="original",
        biotype="domain",
        name="domain1",
        spans=[(20, 25)],
        strand="-",
    )
    with pytest.raises(ValueError):
        # seqid does not exist
        list(seqs.get_features(name="gene1", seqid="blah"))


def test_map():
    """reversing a map with multiple spans should match hand-crafted"""
    forward = [Span(20, 30), Span(40, 50)]
    fmap = FeatureMap(spans=forward, parent_length=100)
    fmap_reversed = fmap.nucleic_reversed()
    reverse = [Span(50, 60), Span(70, 80)]
    rmap = FeatureMap(spans=reverse, parent_length=100)
    assert fmap_reversed.get_coordinates() == rmap.get_coordinates()

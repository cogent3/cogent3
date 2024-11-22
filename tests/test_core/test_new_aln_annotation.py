import os

import pytest

from cogent3 import load_seq
from cogent3.core import new_alignment, new_genetic_code, new_moltype
from cogent3.core.annotation_db import GffAnnotationDb, load_annotations

DNA = new_moltype.get_moltype("dna")


@pytest.fixture(scope="function")
def gff_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "simple.gff")


def makeSampleSequence(name, with_gaps=False):
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


def makeSampleAlignment():
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {seq1.name: seq1, seq2.name: seq2}
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    aln.add_feature(
        biotype="misc_feature", name="misc", spans=[(12, 25)], on_alignment=True
    )
    aln.add_feature(biotype="CDS", name="blue", spans=[(15, 25)], on_alignment=True)
    aln.add_feature(biotype="5'UTR", name="red", spans=[(2, 4)], on_alignment=True)
    aln.add_feature(biotype="LTR", name="fake", spans=[(2, 15)], on_alignment=True)
    return aln


@pytest.fixture()
def ann_aln():
    # synthetic annotated alignment
    return makeSampleAlignment()


@pytest.fixture()
def ann_seq():
    return makeSampleSequence("seq1")


def _make_seq(name):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    seq = new_moltype.DNA.make_seq(seq=raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_init_seqs_have_annotations(mk_cls):
    """annotations on input seqs correctly merged and propagated"""
    data = {"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2")}
    seq_coll = mk_cls(data, moltype="dna")
    coll_db = seq_coll.annotation_db
    assert len(coll_db) == 4
    for name in seq_coll.names:
        seq = seq_coll.get_seq(name)
        db = seq.annotation_db
        assert db is coll_db


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_constructing_collections(mk_cls):
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {"FAKE01": seq1, "FAKE02": seq2}
    expect = sum(s.annotation_db.num_matches() for s in seqs.values())
    coll = mk_cls(seqs, moltype="dna")
    got = coll.annotation_db.num_matches()

    assert got == expect


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_init_annotated_seqs(mk_cls):
    """correctly construct from list with annotated seq"""
    seq = new_moltype.DNA.make_seq(seq="GCCAGGGGGGAAAG-GGAGAA", name="seq1")
    seq.add_feature(biotype="exon", name="name", spans=[(4, 10)])
    coll = mk_cls({"seq1": seq}, moltype="dna")
    features = list(coll.get_features(biotype="exon"))
    assert len(features) == 1


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_sequence_collection_add_feature(mk_cls):
    seqs = mk_cls({"seq1": "AAAAAA", "seq2": "TTTTTT", "seq3": "ATTCCC"}, moltype="dna")
    _ = seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    assert seqs.annotation_db.num_matches() == 1

    # if seqid is not in seqs, raise error
    with pytest.raises(ValueError):
        _ = seqs.add_feature(seqid="bad_seq", biotype="xyz", name="abc", spans=[(1, 2)])


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_sequence_collection_get_annotations_from_any_seq(mk_cls):
    """get_annotations_from_any_seq returns correct annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="dna")
    db = GffAnnotationDb()
    db.add_feature(seqid="seq1", biotype="exon", name="annotation1", spans=[(3, 8)])
    db.add_feature(seqid="seq2", biotype="exon", name="annotation2", spans=[(1, 2)])
    db.add_feature(seqid="seq3", biotype="exon", name="annotation3", spans=[(3, 6)])
    seqs.annotation_db = db
    got = list(seqs.get_features())
    assert len(got) == 3
    assert "biotype='exon', name='annotation1', map=[3:8]/9" in str(got[0])
    assert "biotype='exon', name='annotation2', map=[1:2]/9" in str(got[1])
    assert "biotype='exon', name='annotation3', map=[3:6]/9" in str(got[2])

    got = list(seqs.get_features(name="annotation1"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation1', map=[3:8]/9" in str(got[0])

    got = list(seqs.get_features(biotype="exon", name="annotation2"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation2', map=[1:2]/9" in str(got[0])

    got = list(seqs.get_features(name="annotation3"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation3', map=[3:6]/9" in str(got[0])


@pytest.mark.parametrize("rc", (False, True))
def test_alignment_get_slice(rc):
    """get_slice should return the same as slicing the sequence directly"""
    seq = DNA.make_seq(seq="ATTGTACGCCCCTGA", name="test_seq")
    feature_data = {
        "seqid": "test_seq",
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": "+",
    }
    aln = new_alignment.make_aligned_seqs([seq], moltype="dna")
    aln = aln.rc() if rc else aln
    feature = aln.make_feature(feature=feature_data)
    got_aln = feature.get_slice()
    got_seq = got_aln.get_seq("test_seq")
    assert str(got_seq) == str(seq[5:10])


def test_align_get_features():
    #                    0123456789   the positions
    seq1 = DNA.make_seq(seq="ACG--ACCGT", name="seq1")
    seq2 = DNA.make_seq(seq="ACGGGCCCGT", name="seq2")
    #                      *****      the CDS feature
    seq2.add_feature(biotype="CDS", name="fake01", spans=[(2, 7)], strand="+")
    aln = new_alignment.make_aligned_seqs([seq1, seq2], moltype="dna")
    feat = list(aln.get_features(biotype="CDS"))[0]
    sl = aln[feat]
    # slice is correct length
    assert len(sl) == (7 - 2)
    # returns correct value
    assert sl.to_dict() == dict(seq1="G--AC", seq2="GGGCC")


@pytest.fixture(scope="session")
def seqcoll_db():
    fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    seq = load_seq(fasta_path, moltype="dna", new_type=True, annotation_path=gff3_path)
    return new_alignment.make_unaligned_seqs({seq.name: seq}, moltype="dna")


def test_seqcoll_query(seqcoll_db):
    """querying from a SequenceCollection produces features bound to their seqs"""
    matches = list(seqcoll_db.get_features(seqid="I"))
    # 11 features
    assert len(matches) == 11
    matches = list(seqcoll_db.get_features(seqid="I", biotype="gene"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="mRNA"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="exon"))
    assert len(matches) == 3


def test_feature_projection_ungapped(ann_aln):
    # projection onto ungapped sequence
    expecteds = {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"}
    aln_ltr = list(ann_aln.get_features(biotype="LTR"))[0]
    seq_name = "FAKE01"
    expected = expecteds[seq_name]
    num = ann_aln.annotation_db.num_matches()
    seq_ltr = ann_aln.get_projected_feature(seqid=seq_name, feature=aln_ltr)
    assert ann_aln.annotation_db.num_matches() == num + 1
    got = seq_ltr.get_slice()
    assert str(got) == expected
    assert seq_ltr.seqid == seq_name
    assert seq_ltr.parent == ann_aln.get_seq(seq_name)


def test_feature_projection_gapped(ann_aln):
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


def test_get_feature():
    aln = new_alignment.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--CCCT"}, moltype="dna"
    )

    db = aln.annotation_db
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    feat = list(aln.get_features(seqid="y", biotype="exon", on_alignment=False))[0]
    sliced = feat.get_slice()
    assert sliced.to_dict() == dict(x="AAA", y="CCT")


def test_aln_feature_to_dict():
    seqs = [
        makeSampleSequence("s1", with_gaps=False),
        makeSampleSequence("s2", with_gaps=True),
    ]
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
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


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize(
    "make_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_features_survives_aligned_seq_rename(rev, make_cls):
    segments = ["A" * 10, "C" * 10, "T" * 5, "C" * 5, "A" * 5]

    seqs = make_cls({"original": "".join(segments)}, moltype="dna")
    seqs.annotation_db.add_feature(
        seqid="original", biotype="gene", name="gene1", spans=[(10, 20), (25, 30)]
    )
    seqs.annotation_db.add_feature(
        seqid="original", biotype="domain", name="domain1", spans=[(20, 25)], strand="-"
    )
    seqs = seqs.rename_seqs(lambda x: "newname")
    assert seqs.names == ["newname"]
    seqs = seqs.rc() if rev else seqs
    # quite different behaviour from Alignment and SequenceCollection
    # so we convert to string to make comparison simpler
    got = list(seqs.get_features(name="gene1"))[0]
    sliced = str(got.get_slice()).splitlines()[-1]
    assert sliced == "C" * 15


def test_alignment_annotations():
    """test that annotations to alignment and its' sequences"""

    aln = makeSampleAlignment()
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
            list(aln.get_features(biotype=annot_type, on_alignment=True))[0]
            .get_slice()
            .to_dict()
        )
        expected = aln_expecteds[annot_type]
        assert observed == expected, (annot_type, expected, observed)
        if annot_type in ["misc_feature", "LTR"]:
            continue  # because seqs haven't been annotated with it

        observed = list(aln.get_features(seqid=aln.names, biotype=annot_type))[0]
        observed = observed.get_slice().to_dict()
        expected = seq_expecteds[annot_type]
        assert observed == expected


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_add_to_seq_updates_coll(mk_cls):
    """annotating a seq updates the db of the propagated"""
    data = {
        "x": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
        "y": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
    }
    seq_coll = mk_cls(
        data,
        moltype="dna",
    )
    x = seq_coll.get_seq("x")
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 0
    x.add_feature(biotype="exon", name="E1", spans=[(3, 8)])
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 1


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_annotation_db_assign_none(mk_cls):
    """assigning None to annotation_db breaks conection"""
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = mk_cls(data, moltype="rna")
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    assert seq_coll.annotation_db is not None
    seq_coll.annotation_db = None
    assert seq_coll.annotation_db is None


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_annotation_db_assign_same(gff_db, mk_cls):
    """assigning the same annotation_db"""
    data = {"test_seq": "ACGT", "test_seq2": "CGTT"}
    seq_coll = mk_cls(data, moltype="dna", annotation_db=gff_db)
    seq_coll.annotation_db = gff_db
    assert seq_coll.annotation_db is gff_db


@pytest.mark.parametrize("rved", (False, True))
@pytest.mark.parametrize("annot_type", ("LTR", "misc_feature", "CDS", "5'UTR"))
def test_region_union_on_alignment(annot_type, rved):
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
    newaln = aln.rc() if rved else aln
    feature_list = list(newaln.get_features(biotype=annot_type, on_alignment=True))
    new = feature_list[0]
    new = new.get_slice()
    got = new.to_dict()
    expected = aln_expecteds[annot_type]
    assert expected == got, (annot_type, expected, got)


def test_annotate_matches_to():
    """Aligned.annotate_matches_to correctly delegates to sequence"""

    aln = new_alignment.make_aligned_seqs(dict(x="TTCCACTTCCGCTT"), moltype="dna")
    aln.annotation_db = GffAnnotationDb()
    seq = aln.seqs["x"]
    pattern = "CCRC"
    annot = seq.annotate_matches_to(
        pattern=pattern, biotype="domain", name="fred", allow_multiple=True
    )
    got = [str(a.get_slice()) for a in annot]
    matches = ["CCAC", "CCGC"]
    assert got == matches
    annot = seq.annotate_matches_to(
        pattern=pattern, biotype="domain", name="fred", allow_multiple=False
    )
    got = [a.get_slice() for a in annot]
    assert got == matches[:1]

    # handles regex from aa
    aln = new_alignment.make_aligned_seqs(dict(x="TTCCACTTCCGCTT"), moltype="dna")
    aln.annotation_db = GffAnnotationDb()
    gc = new_genetic_code.get_code(1)
    aa_regex = gc.to_regex("FHF")
    aln.seqs["x"].annotate_matches_to(aa_regex, "domain", "test", allow_multiple=False)
    a = list(aln.get_features(seqid="x"))[0]
    assert str(aln[a].seqs["x"]) == "TTCCACTTC"


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_features_invalid_seqid(mk_cls):
    segments = ["A" * 10, "C" * 10, "T" * 5, "C" * 5, "A" * 5]

    seqs = mk_cls({"original": "".join(segments)}, moltype="dna")
    seqs.annotation_db.add_feature(
        seqid="original", biotype="domain", name="domain1", spans=[(20, 25)], strand="-"
    )
    with pytest.raises(ValueError):
        # seqid does not exist
        list(seqs.get_features(name="gene1", seqid="blah"))


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_copy_annotations(gff_db, mk_cls):
    """copy_annotations copies records from annotation db"""
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = mk_cls(data, moltype="rna")
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches() + gff_db.num_matches()
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect

    # copy annotations with no current annotations
    seq_coll = mk_cls(data, moltype="rna")
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_copy_annotations_same_annotations(gff_db, mk_cls):
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = mk_cls(data, moltype="rna")
    seq_coll = mk_cls(data, moltype="rna")

    # copy annotations with the same annotation_db
    seq_coll.annotation_db = gff_db
    seq_coll.copy_annotations(gff_db)

    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_copy_annotations_none_matching(gff_db, mk_cls):
    """copy annotations should old copy annotations for matching seqids"""
    data = {"name1": "ACGU", "name2": "CGUA", "name_3": "CCGU"}
    seq_coll = mk_cls(data, moltype="rna")
    seq_coll.add_feature(seqid="name1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="name2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches()
    assert gff_db.num_matches() > 0
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_copy_annotations_no_db(gff_db, mk_cls):
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = mk_cls(data, moltype="rna")

    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


def test_project_features_onto_specified_seqid():
    db = GffAnnotationDb()
    aln = new_alignment.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--TTTT"}, moltype="dna"
    )
    db.add_feature(seqid="x", biotype="exon", name="fred", spans=[(3, 8)])
    aln.annotation_db = db

    exons = aln.get_projected_features(seqid="y", biotype="exon")
    assert len(exons) == 1
    assert str(aln.get_seq("y")[exons[0].map.without_gaps()]), "TTT"
    assert "biotype='exon', name='fred', map=[-2-, 4:7]/8" in str(exons[0])


def test_project_features_no_features_for_specified_seqid():
    db = GffAnnotationDb()
    aln = new_alignment.make_aligned_seqs(
        {"seq1": "ATGCGT", "seq2": "AT--GT"}, moltype="dna"
    )
    db.add_feature(seqid="seq1", biotype="gene", name="gene1", spans=[(0, 6)])
    aln.annotation_db = db

    projected_features = aln.get_projected_features(seqid="seq2", biotype="exon")

    assert len(projected_features) == 0


def test_annotated_region_masks():
    """masking a sequence with specific features"""
    # ported from test_features.py
    # refactor: break into multiple tests

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
    aln = new_alignment.make_aligned_seqs(orig_data, moltype="dna")
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


def test_masking_strand_agnostic_aln():
    # ported from test_features.py
    db = GffAnnotationDb()
    db.add_feature(
        seqid="x", biotype="CDS", name="gene", spans=[(2, 6), (10, 15), (25, 35)]
    )
    # Annotations should be correctly masked,
    # whether the sequence has been reverse complemented or not.
    # We use the plus/minus strand CDS containing sequences created above.

    aln = new_alignment.make_aligned_seqs(
        {
            "x": "AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA",
            "y": "AAGGGGAAAACCCCCGGGGGGGGGGTTTTTTTTTTAAA",
        },
        moltype="dna",
        annotation_db=db,
    )
    masked = aln.with_masked_annotations("CDS")
    assert masked.to_dict() == {
        "x": "AA????AAAA?????AAAAAAAAAA??????????AAA",
        "y": str(aln.seqs["y"]),
    }
    rc = aln.rc()
    masked = rc.with_masked_annotations("CDS")
    assert masked.to_dict() == {
        "x": "TTT??????????TTTTTTTTTT?????TTTT????TT",
        "y": str(rc.seqs["y"]),
    }


def test_nested_annotated_region_masks():
    """masking a sequence with specific features when nested annotations"""
    # ported from test_features.py
    db = GffAnnotationDb()
    db.add_feature(seqid="x", biotype="gene", name="norwegian", spans=[(0, 4)])
    db.add_feature(seqid="x", biotype="repeat", name="blue", spans=[(1, 3)])
    db.add_feature(seqid="y", biotype="repeat", name="frog", spans=[(1, 4)])
    aln = new_alignment.make_aligned_seqs(
        [["x", "C-GGCAAAAATTTAA"], ["y", "-T----TTTTG-GTT"]],
        annotation_db=db,
        moltype="text",
    )
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

import pytest

from cogent3 import load_seq
from cogent3.core import annotation_db as anndb_module
from cogent3.core import moltype as c3_moltype

DNA = c3_moltype.get_moltype("dna")
ASCII = c3_moltype.get_moltype("text")


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


@pytest.fixture
def ann_seq():
    return makeSampleSequence("seq1")


def test_seq_feature_to_dict():
    """create the attributes necessary to write into the user table"""
    seq = DNA.make_seq(seq="ATTGTACGCCCCTGA", name="test_seq")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": 1,
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


@pytest.mark.parametrize("rev", [False, True])
def test_features_survives_seq_rename(rev):
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

    got = next(iter(sliced.get_features(name="gene1")))
    got = got.get_slice()
    assert str(got) == gene_expect

    got = next(iter(sliced.get_features(name="domain1")))
    got = got.get_slice()
    assert str(got) == domain_expect


def test_annotate_matches_to():
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


@pytest.mark.parametrize("annot_type", ["CDS", "5'UTR"])
def test_slice_seq_with_full_annotations(ann_seq, annot_type):
    # this slice contains both features intact
    newseq = ann_seq[10:]
    orig = next(iter(ann_seq.get_features(biotype=annot_type)))
    new = next(iter(newseq.get_features(biotype=annot_type)))
    assert orig.name == new.name
    assert len(orig) == len(new)
    assert str(newseq[new]) == str(ann_seq[orig]), annot_type


@pytest.mark.parametrize(("annot_type", "num"), [("CDS", 0), ("5'UTR", 1)])
def test_slice_seq_with_partial_end(ann_seq, annot_type, num):
    # this slice contains both features intact
    newseq = ann_seq[:14]
    # only UTR is present
    new = list(newseq.get_features(biotype=annot_type, allow_partial=True))
    assert len(new) == num, annot_type
    if num:
        feat = new[0]
        # length of the feature is the same as the original
        assert len(feat) == len(next(iter(ann_seq.get_features(biotype=annot_type))))
        gapless = feat.without_lost_spans()
        # the sliced feature without gaps is shorter
        assert len(gapless) < len(feat)


@pytest.mark.parametrize(("annot_type", "num"), [("CDS", 1), ("5'UTR", 0)])
def test_slice_seq_with_partial_start(ann_seq, annot_type, num):
    # this slice contains both features intact
    newseq = ann_seq[18:]
    # only UTR is present
    new = list(newseq.get_features(biotype=annot_type, allow_partial=True))
    assert len(new) == num, annot_type
    if num:
        feat = new[0]
        # length of the feature is the same as the original
        assert len(feat) == len(next(iter(ann_seq.get_features(biotype=annot_type))))
        gapless = feat.without_lost_spans()
        # the sliced feature without gaps is shorter
        assert len(gapless) < len(feat)


def test_gbdb_get_children_get_parent(DATA_DIR):
    seq = load_seq(DATA_DIR / "annotated_seq.gb", moltype="dna")
    seq = seq[2900:6000]
    (orig,) = list(seq.get_features(biotype="gene", name="CNA00110"))
    (child,) = list(orig.get_children("CDS"))
    parent, *_ = list(child.get_parent())
    assert parent == orig


def test_feature_names_seq():
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    seq = DNA.make_seq(seq=raw_seq, name="s1")
    f = seq.add_feature(biotype="CDS", name="s1-cds", spans=[cds])
    assert f.name == "s1-cds"
    # apply the sequence name to the feature
    s = f.get_slice(apply_name=False)
    assert s.name == seq.name
    # or apply the feature name
    s = f.get_slice(apply_name=True)
    assert s.name == "s1-cds" != seq.name
    # default behaviour is to apply the feature name
    s = seq[f]
    assert s.name == "s1-cds"
    assert s.name == f.name


def test_seq_with_masked_annotations():
    seq = makeSampleSequence("seq1", with_gaps=False)
    raw_seq = str(seq)
    cds = next(iter(seq.annotation_db.get_records_matching(biotype="CDS")))
    start, stop = cds["start"], cds["stop"]
    expect = raw_seq[:start] + "?" * (stop - start) + raw_seq[stop:]
    masked = seq.with_masked_annotations(biotypes="CDS")
    assert str(masked) == expect
    shadow = seq.with_masked_annotations(biotypes="CDS", shadow=True)
    expect = "?" * start + raw_seq[start:stop] + "?" * (len(raw_seq) - stop)
    assert str(shadow) == expect


@pytest.mark.parametrize("shadow", [True, False])
def test_seq_with_masked_annotations_missing_feature(shadow):
    seq = makeSampleSequence("seq1", with_gaps=False)
    raw_seq = str(seq)
    masked = seq.with_masked_annotations(biotypes="not-present", shadow=shadow)
    assert str(masked) == raw_seq


def test_is_annotated():
    """is_annotated operates correctly"""
    s = c3_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG", name="s1")
    _ = s.add_feature(biotype="gene", name="blah", spans=[(0, 10)])
    assert s.is_annotated()


@pytest.mark.parametrize("biotype", ["gene", "exon", ("gene", "exon")])
def test_is_annotated_biotype(biotype):
    """is_annotated operates correctly"""
    s = c3_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG", name="s1")
    _ = s.add_feature(biotype="gene", name="blah", spans=[(0, 10)])
    _ = s.add_feature(biotype="exon", name="blah", spans=[(0, 10)])
    assert s.is_annotated(biotype=biotype)


def test_not_is_annotated():
    """is_annotated operates correctly"""
    s = c3_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG", name="s1")
    assert not s.is_annotated()
    # annotation on different seq
    s.annotation_db.add_feature(
        seqid="s2",
        biotype="gene",
        name="blah",
        spans=[(0, 10)],
    )
    assert not s.is_annotated()
    # annotation wrong biotype
    s.annotation_db.add_feature(
        seqid="s1",
        biotype="exon",
        name="blah",
        spans=[(0, 10)],
    )
    assert not s.is_annotated(biotype="gene")
    s.annotation_db = None
    assert not s.is_annotated()


def test_annotation_db_lazy_evaluation():
    s = c3_moltype.DNA.make_seq(seq="AC", name="s1")
    assert isinstance(s._annotation_db, list)
    # now if we invoke the property we get an actual db instance created
    assert isinstance(s.annotation_db, anndb_module.SupportsFeatures)


def test_init_with_annotationdb():
    anndb = anndb_module.GffAnnotationDb()
    s = c3_moltype.DNA.make_seq(seq="AC", name="s1", annotation_db=anndb)
    assert isinstance(s.annotation_db, anndb_module.GffAnnotationDb)
    assert s.annotation_db is anndb


def test_init_with_annotation_offset():
    s = c3_moltype.DNA.make_seq(seq="AC", name="s1", annotation_offset=2)
    assert s.annotation_offset == 2


def test_init_with_annotation_offset_sliced():
    s = c3_moltype.DNA.make_seq(seq="ACTTTGG", name="s1", annotation_offset=2)
    sl = s[2:5]
    assert sl.annotation_offset == 4
    _, start, stop, _ = sl.parent_coordinates()
    assert start == 4
    assert stop == 7


def test_init_with_annotation_offset_plus_strand():
    annotation_offset = 2000
    plus_raw = "ACTTTGGCC"
    s = c3_moltype.DNA.make_seq(
        seq=plus_raw, name="s1", annotation_offset=annotation_offset
    )
    rel_start = 2
    rel_stop = 5
    expect = plus_raw[rel_start:rel_stop]
    db = s.annotation_db
    db.add_feature(
        biotype="gene",
        name="blah",
        spans=[(annotation_offset + rel_start, annotation_offset + rel_stop)],
        seqid="s1",
    )
    ft = list(s.get_features(biotype="gene"))[0]
    assert str(ft.get_slice()) == expect


def test_init_with_annotation_offset_minus_strand():
    annotation_offset = 2000
    plus_raw = "ACTTTGGCC"
    rc = c3_moltype.DNA.make_seq(
        seq=plus_raw, name="s1", annotation_offset=annotation_offset
    ).rc()
    rel_start = 3
    rel_stop = 6
    expect = plus_raw[rel_start:rel_stop]
    db = rc.annotation_db
    db.add_feature(
        biotype="gene",
        name="blah",
        spans=[(annotation_offset + rel_start, annotation_offset + rel_stop)],
        strand=-1,
        seqid="s1",
    )
    ft = list(rc.get_features(biotype="gene"))[0]
    assert str(ft.get_slice()) == c3_moltype.DNA.rc(expect)

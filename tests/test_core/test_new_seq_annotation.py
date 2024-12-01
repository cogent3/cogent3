import pytest

from cogent3 import load_seq
from cogent3.core import new_moltype

DNA = new_moltype.get_moltype("dna")
ASCII = new_moltype.get_moltype("text")


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


@pytest.mark.parametrize("rev", (False, True))
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

    got = list(sliced.get_features(name="gene1"))[0]
    got = got.get_slice()
    assert str(got) == gene_expect

    got = list(sliced.get_features(name="domain1"))[0]
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


@pytest.mark.parametrize("annot_type", ("CDS", "5'UTR"))
def test_slice_seq_with_full_annotations(ann_seq, annot_type):
    # this slice contains both features intact
    newseq = ann_seq[10:]
    orig = list(ann_seq.get_features(biotype=annot_type))[0]
    new = list(newseq.get_features(biotype=annot_type))[0]
    assert orig.name == new.name
    assert len(orig) == len(new)
    assert str(newseq[new]) == str(ann_seq[orig]), annot_type


@pytest.mark.parametrize("annot_type,num", (("CDS", 0), ("5'UTR", 1)))
def test_slice_seq_with_partial_end(ann_seq, annot_type, num):
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


@pytest.mark.parametrize("annot_type,num", (("CDS", 1), ("5'UTR", 0)))
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


def test_gbdb_get_children_get_parent(DATA_DIR):
    seq = load_seq(DATA_DIR / "annotated_seq.gb", new_type=True)
    seq = seq[2900:6000]
    (orig,) = list(seq.get_features(biotype="gene", name="CNA00110"))
    (child,) = list(orig.get_children("CDS"))
    parent, *_ = list(child.get_parent())
    assert parent == orig

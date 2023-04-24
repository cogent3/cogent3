import pathlib

import numpy
import pytest

from cogent3 import DNA, RNA, SequenceCollection, _Table, load_seq
from cogent3.core.annotation import Annotation
from cogent3.core.annotation_db import (
    GffAnnotationDb,
    SupportsFeatures,
    _matching_conditions,
    load_annotations,
)
from cogent3.core.sequence import Sequence


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Katherine Caley"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

DATA_DIR = pathlib.Path(__file__).parent.parent / "data"


@pytest.fixture(scope="session")
def gff_db():
    paths = (
        "sample_data/Homo_sapiens.GRCh38.109.chromosome.22.gff3.gz",
        "data/c_elegans_WS199_shortened_gff.gff3",
    )
    return load_annotations(paths[1])


@pytest.fixture()
def seq_db():
    seq = load_seq(DATA_DIR / "c_elegans_WS199_dna_shortened.fasta", moltype="dna")
    db = load_annotations(DATA_DIR / "c_elegans_WS199_shortened_gff.gff3")

    seq.annotation_db = db

    return seq


@pytest.fixture()
def seq() -> Sequence:
    return Sequence("ATTGTACGCCTTTTTTATTATT", name="test_seq")


@pytest.fixture()
def anno_db() -> GffAnnotationDb:
    # an empty db that we can add to
    return GffAnnotationDb([])


@pytest.fixture()
def simple_seq_gff_db() -> Sequence:
    seq = Sequence("ATTGTACGCCTTTTTTATTATT", name="test_seq")
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    return seq


def test_gff_describe(gff_db):
    result = gff_db.describe
    assert isinstance(result, _Table)


def test_gff_features_matching(gff_db):
    result = list(gff_db.get_features_matching(biotype="CDS"))
    assert len(result)


def test_gff_get_children(gff_db):
    # an ID with 8 children
    children = list(gff_db.get_feature_children(name="Transcript:B0019.1"))
    assert len(children) == 8
    # this is parent of the transcript, and due to the funky structure of
    # this gff, also to one of the CDS entries
    children = list(gff_db.get_feature_children(name="Gene:WBGene00000138"))
    assert len(children) == 2


@pytest.mark.parametrize(
    "name,expected",
    (
        ("CDS:B0019.1", ("Transcript:B0019.1", "Gene:WBGene00000138")),
        ("Transcript:B0019.1", ("Gene:WBGene00000138",)),
    ),
)
def test_gff_get_parent(gff_db, name, expected):
    # the first id has two parents, which is weird
    got = list(gff_db.get_feature_parent(name=name))
    assert len(got) == len(expected)
    assert {g["name"] for g in got} == set(expected)


def test_gff_counts(gff_db):
    got = gff_db.biotype_counts()
    assert len(got) > 0


def test_gff_num_matches(gff_db):
    count = gff_db.num_matches()
    assert count == 11
    assert gff_db.num_matches(seqid="I") == 11
    assert gff_db.num_matches(seqid="IV") == 0


def test_gb_num_matches(gb_db):
    count = gb_db.num_matches()
    assert count == 10  # value from manual count from file
    assert gb_db.num_matches(seqid="AE017341") == 10
    assert gb_db.num_matches(seqid="IV") == 0


def test_gff_find_user_features(gff_db):
    record = dict(
        seqid="2", name="gene-01", biotype="gene", spans=[(23, 33)], strand="+"
    )
    gff_db.add_feature(**record)
    # by biotype
    found = any(
        gene["name"] == "gene-01"
        for gene in gff_db.get_features_matching(biotype="gene")
    )
    assert found
    # by name
    found = any(
        gene["name"] == "gene-01"
        for gene in gff_db.get_features_matching(name="gene-01")
    )
    assert found


def test_to_rich_dict(gff_db):
    # not a very good test! needs to be serialise / deserialize
    rd = gff_db.to_rich_dict()
    assert "type" in rd
    _ = GffAnnotationDb.from_dict(rd)


def test_empty_data():
    _ = GffAnnotationDb([])


# testing GenBank files
@pytest.fixture(scope="session")
def gb_db():
    return load_annotations(DATA_DIR / "annotated_seq.gb")


@pytest.mark.parametrize("parent_biotype, name", (("gene", "CNA00110"),))
def test_gb_get_children(gb_db, parent_biotype, name):
    parent = list(gb_db.get_features_matching(biotype=parent_biotype, name=name))[0]
    coords = numpy.array(parent["spans"])
    child = list(
        gb_db.get_feature_children(
            name=name,
            exclude_biotype=parent_biotype,
            start=coords.min(),
            end=coords.max(),
        )
    )[0]
    assert child["biotype"] != parent["biotype"]
    assert child["name"] == parent["name"]


def test_gb_get_parent(gb_db):
    cds_id = "CNA00110"
    cds = list(gb_db.get_features_matching(biotype="CDS", name=cds_id))[0]
    coords = numpy.array(cds["spans"])
    parent = list(
        gb_db.get_feature_parent(
            name=cds_id,
            exclude_biotype="CDS",
            start=coords.min(),
            end=coords.max(),
        )
    )[0]
    assert parent["biotype"] != cds["biotype"]
    assert parent["biotype"] == "gene"
    assert parent["name"] == cds["name"]


def test_protocol_adherence(gff_db, gb_db):
    for db in (gff_db, gb_db):
        assert isinstance(db, SupportsFeatures)


def test_get_features_matching_no_annotation_db(seq):
    """
    Test that `get_features_matching` returns an empty list when no annotation database is attached to the sequence.
    """
    assert not list(seq.get_features_matching())


def test_get_features_matching_no_matching_feature(seq, anno_db):
    """
    Test that `get_features_matching` returns an empty list when there is no matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )

    assert not list(seq.get_features_matching(feature_type="exon", name="non_matching"))
    assert not list(seq.get_features_matching(feature_type="CDS"))


def test_get_features_matching_matching_feature(seq, anno_db):
    """
    Test that `get_features_matching` returns a list with one matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    got = list(seq.get_features_matching(feature_type="exon"))

    assert got[0].biotype == "exon"
    assert got[0].name == "exon1"
    assert len(got) == 1


def test_gav():
    # gah todo also need to check a slice with stride
    from cogent3 import make_seq
    from cogent3.core.annotation import Feature

    #            ++   ++++
    #              ---    --
    raw_seq = "AACCTTTGGGGAATTT"

    plus_spans = [(2, 4), (7, 11)]
    plus_seq = "".join(raw_seq[s:e] for s, e in plus_spans)

    minus_spans = [(4, 7), (11, 13)]
    minus_seq = "".join(raw_seq[s:e] for s, e in minus_spans)
    minus_seq = "".join([{"T": "A", "A": "T"}[b] for b in minus_seq[::-1]])
    seq = make_seq(raw_seq, name="s1", moltype="dna")
    db = GffAnnotationDb(data=[])
    db.add_feature(
        seqid="s1",
        biotype="cds",
        name="plus",
        spans=plus_spans,
        strand="+",
        on_alignment=False,
    )
    db.add_feature(
        seqid="s1",
        biotype="cds",
        name="minus",
        spans=minus_spans,
        strand="-",
        on_alignment=False,
    )
    seq.annotation_db = db
    plus = list(seq.get_features_matching(name="plus"))[0]
    assert str(plus.get_slice()) == plus_seq
    minus = list(seq.get_features_matching(name="minus"))[0]
    assert str(minus.get_slice()) == minus_seq

    f = Feature(seq, "cds", "minus", reversed([reversed(s) for s in minus_spans]))
    s = f.get_slice()
    assert str(s) == minus_seq

    # now reverse complement the sequence
    rced = seq.rc()

    rced_spans = [(5, 9), (12, 14)]
    rf = Feature(rced, "cds", "minus", [reversed(s) for s in rced_spans])
    got = rf.get_slice()

    plus = list(rced.get_features_matching(name="plus"))[0]
    assert str(plus.get_slice()) == plus_seq
    minus = list(rced.get_features_matching(name="minus"))[0]
    print(rf, minus)
    assert str(minus.get_slice()) == minus_seq


def test_gav2():
    from cogent3 import make_seq
    from cogent3.core import location as loc

    seq = make_seq("AACCTTTGGGGAATTT", moltype="dna")
    mmap = loc.Map(locations=[(4, 7), (11, 13)], parent_length=16)
    expect = seq[mmap.reversed()]

    rcseq = seq.rc()
    rmap = mmap.nucleic_reversed().reversed()
    got = rcseq[rmap]
    assert str(got) == str(expect)


def test_get_features_matching_matching_features(anno_db: GffAnnotationDb, seq):
    """
    Test that `get_features_matching` returns a list with all matching features in the annotation database.
    """
    seq.annotation_db = anno_db

    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon2", spans=[(6, 10)], strand="+"
    )
    got = list(seq.get_features_matching(feature_type="exon"))

    assert len(got) == 2


def test_annotate_from_gff(seq):
    seq.annotate_from_gff(DATA_DIR / "simple.gff")

    got = list(seq.get_features_matching(feature_type="exon"))
    assert len(got) == 2

    got = list(seq.get_features_matching(feature_type="CDS"))
    assert len(got) == 2

    got = list(seq.get_features_matching(feature_type="CpG"))
    assert len(got) == 1


def test_get_features_matching_start_stop(seq):
    # todo: need to implement LostSpans()
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    got = list(seq.get_features_matching(start=2, stop=10))
    assert len(got) == 4


def test_matching_conditions():
    got, _ = _matching_conditions({"start": 1, "end": 5}, partial=True)
    expect = "((start <= 1 AND end > 1) OR (start < 5 AND end >= 5) OR (start <= 1 AND end >= 5))"
    assert got == expect


def test_get_features_matching_start_stop_seqview(seq):
    """testing that get_features_matching adjusts"""
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    seq_features = list(seq.get_features_matching(start=0, stop=3))
    assert len(seq_features) == 3

    # edge case, only 1 features that overlaps with index 12
    # is actually returning [exon2 at [11:20]/13, CpG1 at [2:12]/13]
    # possibly a bug in the SQL generating code
    subseq = seq[9:]
    seq_features_features = list(subseq.get_features_matching(start=3, stop=10))
    assert len(seq_features_features) == 1


def test_get_slice():
    """get_slice should return the same as slicing the sequence directly"""
    seq = Sequence("ATTGTACGCCCCTGA", name="test_seq")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "reversed": False,
    }
    feature = seq.make_feature(feature_data)
    got = feature.get_slice()
    assert str(got) == str(seq[5:10])


def test_feature_get_children(seq_db):
    feat = list(seq_db.get_features_matching(name="Transcript:B0019.1"))[0]
    new_feat_5pUTR = list(feat.get_children(biotype="five_prime_UTR"))
    assert len(new_feat_5pUTR) == 1

    new_feat_CDS = list(feat.get_children(biotype="CDS"))[0]
    assert new_feat_CDS.name == "CDS:B0019.1"


def test_db_persists_post_rc(seq_db):
    """assert that the db persists after the .rc() method call"""
    rc_seq = seq_db.rc()
    assert rc_seq.annotation_db is not None


def test_rc_get_slice_negative_feature(seq_db):
    """given a feature on the - strand, the feature.get_slice() should return
    the same sequence before and after the sequence is reverse complemented
    """

    feat = list(seq_db.get_features_matching(name="Transcript:B0019.1"))[0]
    rc_seq = seq_db.rc()
    r_feat = list(rc_seq.get_features_matching(name="Transcript:B0019.1"))[0]

    assert feat.get_slice() == r_feat.get_slice()


def test_rc_get_slice_positive_feature(anno_db):
    """given a feature on the + strand, the feature.get_slice() should return
    the same sequence before and after the sequence is reverse complemented
    """

    seq = DNA.make_seq("AAAAGGGG", name="seq1")

    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(2, 6)], strand="+"
    )

    feat = list(seq.get_features_matching(name="exon1"))[0]
    r_seq = seq.rc()
    r_feat = list(r_seq.get_features_matching(name="exon1"))[0]

    assert feat.get_slice() == r_feat.get_slice()


def test_add_feature(seq):
    """Sequence supports manual adding of features for a seq with no bound AnnotationDb"""
    record = dict(name="gene-01", biotype="gene", spans=[(12, 16)], strand="+")
    seq.add_feature(**record)
    feats = list(seq.get_features_matching(feature_type="gene"))

    assert seq.annotation_db is not None
    assert len(feats) == 1
    assert feats[0].biotype == "gene"


def test_add_feature_existing_db(simple_seq_gff_db):
    """Sequence supports manual adding of features for a seq with an existing AnnotationDb"""
    record = dict(name="gene-01", biotype="gene", spans=[(12, 16)], strand="+")
    simple_seq_gff_db.add_feature(**record)

    # total features should be 5+1=6
    all_feats = list(simple_seq_gff_db.get_features_matching())
    assert len(all_feats) == 6


def test__getitem__(simple_seq_gff_db):
    """Sequence.__getitem__ should keep the underlying seq in the SeqView
    and preserve any annotation_db"""

    seq_sliced = simple_seq_gff_db[4:6]
    assert seq_sliced == str(simple_seq_gff_db)[4:6]
    # check the underlying seq is still the original sequence data
    assert seq_sliced._seq.seq == str(simple_seq_gff_db)
    # check the annotation_db is still attached and the same instance
    assert (
        seq_sliced.annotation_db
        and seq_sliced.annotation_db is simple_seq_gff_db.annotation_db
    )


def test_to_moltype_dna():
    """to_moltype("dna") ensures conversion from T to U"""
    seq = DNA.make_seq("AAAAGGGGTTT", name="seq1")
    rna = seq.to_moltype("rna")

    assert "T" not in rna


def test_to_moltype_rna():
    """to_moltype("rna") ensures conversion from U to T"""
    seq = RNA.make_seq("AAAAGGGGUUU", name="seq1")
    rna = seq.to_moltype("dna")

    assert "U" not in rna


def test_annotate_from_gff_multiple_calls(seq):
    """5 records in each gff file, total features on seq should be 10"""
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    seq.annotate_from_gff(DATA_DIR / "simple2.gff")
    assert len(list(seq.get_features_matching())) == 10


def test_sequence_collection_annotate_from_gff():
    """providing a seqid to SequenceCollection.annotate_from_gff will
    annotate the SequenceCollection, and the Sequence. Both of these will point
    to the same AnnotationDb instance
    """
    seqs = {"test_seq": "ATCGATCGATCG", "test_seq2": "GATCGATCGATC"}
    seq_collection = SequenceCollection(seqs)
    seq_collection.annotate_from_gff(DATA_DIR / "simple.gff", seq_ids="test_seq")

    # the seq for which the seqid was provided is annotated
    seq = seq_collection.get_seq("test_seq")
    assert seq_collection.get_seq("test_seq").annotation_db is not None
    assert len(list(seq_collection.get_seq("test_seq").get_features_matching())) == 5
    # the seq for which the seqid was NOT provided is NOT annotated
    assert seq_collection.get_seq("test_seq2").annotation_db is None
    # the annotation_db on the seq and the seq collection are the same object
    assert (
        seq_collection.get_seq("test_seq").annotation_db is seq_collection.annotation_db
    )
    got = list(seq.get_features_matching(feature_type="CDS"))
    assert len(got) == 2

    got = list(seq.get_features_matching(feature_type="CpG"))
    assert len(got) == 1

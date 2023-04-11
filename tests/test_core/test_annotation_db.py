import pathlib

import numpy
import pytest

from cogent3 import _Table, load_seq
from cogent3.core.annotation import FeatureNew
from cogent3.core.annotation_db import (
    GffAnnotationDb,
    SupportsFeatures,
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
    got = GffAnnotationDb.from_dict(rd)


def test_empty_data():
    got = GffAnnotationDb([])


# testing GenBank files
@pytest.fixture(scope="session")
def gb_db():
    paths = ("data/annotated_seq.gb",)
    return load_annotations(paths[0])


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


@pytest.fixture()
def seq() -> Sequence:
    return Sequence("ATTGTACGCCTTTTTTATTATT", name="test_seq")


@pytest.fixture()
def anno_db() -> GffAnnotationDb:
    # an empty db that we can add to
    return GffAnnotationDb([])


def test_query_db_no_annotation_db(seq):
    """
    Test that `query_db` returns an empty list when no annotation database is attached to the sequence.
    """
    assert not list(seq.query_db(feature_type="exon", name="test"))


def test_query_db_no_matching_feature(seq, anno_db):
    """
    Test that `query_db` returns an empty list when there is no matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )

    assert not list(seq.query_db(feature_type="exon", name="non_matching"))
    assert not list(seq.query_db(feature_type="CDS"))


def test_query_db_matching_feature(seq, anno_db):
    """
    Test that `query_db` returns a list with one matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    got = list(seq.query_db(feature_type="exon"))

    assert got[0].biotype == "exon"
    assert got[0].name == "exon1"
    assert len(got) == 1
    assert got[0].reversed == False


def test_query_db_matching_features(anno_db: GffAnnotationDb, seq):
    """
    Test that `query_db` returns a list with all matching features in the annotation database.
    """
    seq.annotation_db = anno_db

    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon2", spans=[(6, 10)], strand="+"
    )
    got = list(seq.query_db(feature_type="exon"))

    assert len(got) == 2


def test_annotate_from(seq):
    seq.annotate_from("data/simple.gff", pre_parsed=False)

    got = list(seq.query_db(feature_type="exon"))
    assert len(got) == 2

    feature1 = got[0]
    assert feature1.name == "exon1"
    assert feature1.biotype == "exon"
    assert (feature1.map.start, feature1.map.end) == (1, 10)


def test_query_db_start_stop(seq):
    # todo: cannot forget the lost spans...
    seq.annotate_from("data/simple.gff")
    got = list(seq.query_db(start=2, stop=10))
    assert len(got) == 4


def test_query_db_start_stop_seqview(seq):
    seq.annotate_from("data/simple.gff", pre_parsed=False)

    subseq = seq[9:]
    got = list(subseq.query_db(start=2, stop=10))
    # the adjusted query should be .query_db(start=9+2, stop=9+10)
    assert (got[0].map.start, got[0].map.end) == (11, 20)


def test_feature_get_slice():
    seq = Sequence("ATTGTACGCCCCTGA", name="test_seq")
    feature_dict = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "reversed": False,
    }

    feature = FeatureNew(seq, **feature_dict)

    got = feature.get_slice()
    assert str(got) == str(seq[5:10])


def test_feature_query_get_slice(seq_db):
    feat = list(seq_db.query_db(name="Transcript:B0019.1"))[0]
    assert feat.name == "Transcript:B0019.1"
    assert str(feat.get_slice()) == str(seq_db)[9:70]


def test_feature_get_children(seq_db):
    feat = list(seq_db.query_db(name="Transcript:B0019.1"))[0]
    new_feat_5pUTR = list(feat.get_children(biotype="five_prime_UTR"))
    assert len(new_feat_5pUTR) == 1
    assert str(new_feat_5pUTR[0].get_slice()) == str(seq_db)[:9]

    new_feat_CDS = list(feat.get_children(biotype="CDS"))[0]
    assert new_feat_CDS.name == "CDS:B0019.1"


def test_db_rc_persists(seq_db):
    """assert that the db persists after the .rc() method call"""
    rc_seq = seq_db.rc()
    assert rc_seq.annotation_db is not None


def test_same_feature_rc(seq_db):
    # Transcript:B0019.1 is a feature on the reverse strand

    feat = list(seq_db.query_db(name="Transcript:B0019.1"))[0]
    rc_seq = seq_db.rc()
    r_feat = list(rc_seq.query_db(name="Transcript:B0019.1"))[0]

    assert feat.get_slice() == r_feat.get_slice()


def test_rc_features(anno_db):
    # adding the feature to the positive strand
    from cogent3 import DNA

    seq = DNA.make_seq("AAAAGGGG", name="seq1")

    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(2, 6)], strand="+"
    )

    feat = list(seq.query_db(name="exon1"))[0]

    r_seq = seq.rc()
    r_feat = list(r_seq.query_db(name="exon1"))[0]

    assert feat.get_slice() == r_feat.get_slice()


def test_seq_getitem():
    from cogent3 import DNA

    seq = DNA.make_seq("AAAAGGGG", name="seq1")

    seq_sliced = seq[4:6]
    assert seq_sliced == str(seq)[4:6:]
    assert seq_sliced._seq.seq == str(seq)


def test_to_moltype():
    from cogent3 import DNA

    seq = DNA.make_seq("AAAAGGGGTTT", name="seq1")
    s = DNA.make_seq(seq)
    rna = s.to_moltype("rna")

    assert "T" not in rna

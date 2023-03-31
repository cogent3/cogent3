import numpy
import pytest

from cogent3 import _Table
from cogent3.core.annotation_db import (
    FeatureDataType,
    GffAnnotationDb,
    load_annotations,
    SupportsFeatures,
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


@pytest.fixture(scope="session")
def gff_db():
    paths = (
        "sample_data/Homo_sapiens.GRCh38.109.chromosome.22.gff3.gz",
        "data/c_elegans_WS199_shortened_gff.gff3",
    )
    return load_annotations(paths[1])


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


@pytest.mark.xfail(reason="todo: upload data to be able to test this method")
@pytest.mark.parametrize("parent_biotype, name", (("gene", "GS_RS00055"),))
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


@pytest.mark.xfail(reason="todo: upload data to be able to test this method")
def test_gb_get_parent(gb_db):
    cds_id = "GS_RS00055"
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
    return Sequence("ATTGTACGC", name="test_seq")


@pytest.fixture()
def anno_db() -> GffAnnotationDb:
    return GffAnnotationDb([])


def test_query_db_no_annotation_db(seq):
    """
    Test that `query_db` returns an empty list when no annotation database is attached to the sequence.
    """
    assert not list(seq.query_db(biotype="exon", seqid="test"))


def test_query_db_no_matching_feature(seq, anno_db):
    """
    Test that `query_db` returns an empty list when there is no matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )

    assert not list(seq.query_db(biotype="exon", name="non_matching"))


def test_query_db_matching_feature(seq, anno_db):
    """
    Test that `query_db` returns a list with one matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    got = list(seq.query_db(biotype="exon", seqid=seq.name))
    expected = FeatureDataType(
        biotype="exon", name="exon1", spans=[(1, 4)], reverse=False
    )

    assert got == [expected]


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
    got = list(seq.query_db(biotype="exon"))
    expected = [
        FeatureDataType(biotype="exon", name="exon1", spans=[(1, 4)], reverse=False),
        FeatureDataType(biotype="exon", name="exon2", spans=[(6, 10)], reverse=False),
    ]

    assert got == expected


def test_annotate_from(seq):
    seq.annotate_from("data/simple.gff", pre_parsed=False)

    got = list(seq.query_db(biotype="exon"))
    assert len(got) == 2

    feature1 = got[0]
    assert feature1["name"] == "exon1"
    assert feature1["biotype"] == "exon"
    assert feature1["spans"] == [(1, 10)]


def test_query_db_start_stop():
    seq = Sequence("ATTGTACGCCCCTGA", name="test_seq")
    seq.annotate_from("data/simple.gff", pre_parsed=False)
    got = list(seq.query_db(start=2, stop=10))


def test_query_db_start_stop_seqview():
    seq = Sequence("A" * 22, name="test_seq")
    seq.annotate_from("data/simple.gff", pre_parsed=False)

    subseq = seq[9:]
    got = list(subseq.query_db(start=2, stop=10))
    # the adjust query should be .query_db(start=9+2, stop=9+10)
    expect = [
        {"biotype": "exon", "name": "exon2", "spans": [(11, 20)], "reverse": False}
    ]
    assert got == expect

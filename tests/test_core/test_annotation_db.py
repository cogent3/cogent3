import numpy
import pytest

from cogent3 import DNA, SequenceCollection, _Table, load_seq
from cogent3.core.annotation_db import (
    BasicAnnotationDb,
    GenbankAnnotationDb,
    GffAnnotationDb,
    SupportsFeatures,
    _matching_conditions,
    load_annotations,
)
from cogent3.core.sequence import Sequence
from cogent3.parse.genbank import MinimalGenbankParser
from cogent3.util import deserialise


@pytest.fixture(scope="function")
def gff_db(DATA_DIR):
    path = DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"
    return load_annotations(path=path)


@pytest.fixture(scope="function")
def gff_small_db(DATA_DIR):
    path = DATA_DIR / "simple.gff"
    return load_annotations(path=path)


@pytest.fixture()
def seq_db(DATA_DIR):
    seq = load_seq(DATA_DIR / "c_elegans_WS199_dna_shortened.fasta", moltype="dna")
    db = load_annotations(path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3")

    seq.annotation_db = db

    return seq


@pytest.fixture()
def seq() -> Sequence:
    return Sequence("ATTGTACGCCTTTTTTATTATT", name="test_seq")


@pytest.fixture(scope="function")
def anno_db() -> BasicAnnotationDb:
    # an empty db that we can add to
    return BasicAnnotationDb()


@pytest.fixture()
def simple_seq_gff_db(DATA_DIR) -> Sequence:
    seq = Sequence("ATTGTACGCCTTTTTTATTATT", name="test_seq")
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    return seq


def test_assign_valid_db(seq, anno_db):
    # should not fail
    seq.annotation_db = anno_db
    assert seq.annotation_db is anno_db


def test_replace_annotation_db_check_invalid(seq):
    with pytest.raises(TypeError):
        seq.replace_annotation_db(2, check=True)


def test_replace_annotation_db_nocheck_invalid(seq):
    seq.replace_annotation_db(2, check=False)
    assert seq.annotation_db == 2


@pytest.mark.parametrize(
    "db_name,cls", (("gff_db", GenbankAnnotationDb), ("gb_db", GffAnnotationDb))
)
def test_constructor_db_fail(db_name, cls, request):
    db = request.getfixturevalue(db_name)
    with pytest.raises(TypeError):
        cls(db=db)


@pytest.mark.parametrize(
    "db_name,cls", (("gff_db", GenbankAnnotationDb), ("gb_db", GffAnnotationDb))
)
def test_constructor_wrong_db_schema(db_name, cls, request):
    db = request.getfixturevalue(db_name)
    with pytest.raises(TypeError):
        cls(db=db.db)


@pytest.mark.parametrize(
    "db_name,cls",
    (
        ("gff_db", GffAnnotationDb),
        ("anno_db", GffAnnotationDb),
        ("gb_db", GenbankAnnotationDb),
        ("anno_db", GenbankAnnotationDb),
    ),
)
def test_constructor_db_instance_works(db_name, cls, request):
    # only compatible db's used to init
    db = request.getfixturevalue(db_name)
    cls(db=db)


@pytest.mark.parametrize(
    "db_name,cls",
    (
        ("gff_db", GffAnnotationDb),
        ("anno_db", GffAnnotationDb),
        ("gb_db", GenbankAnnotationDb),
        ("anno_db", GenbankAnnotationDb),
    ),
)
def test_constructor_db_connection_works(db_name, cls, request):
    # only compatible db's used to init
    db = request.getfixturevalue(db_name)
    cls(db=db.db)


def test_gff_describe(gff_db):
    result = gff_db.describe
    assert isinstance(result, _Table)


def test_count_distinct(gff_db):
    # no arguments, returns None
    assert gff_db.count_distinct() is None

    # there are 8 biotypes in the c.elegans gff sample, 2 columns
    # all arguments returns, from our example, all the rows
    got = gff_db.count_distinct(biotype=True)
    assert got.shape == (
        8,
        2,
    )
    # all names unique, 11 records, 4 columns
    got = gff_db.count_distinct(biotype=True, seqid=True, name=True)
    assert got.shape == (11, 4)


def test_count_distinct_values(gb_db):
    # there are 8 biotypes in the c.elegans gff sample, 2 columns
    # all arguments returns, from our example, all the rows
    got = {tuple(r) for r in gb_db.count_distinct(name=True).to_list()}
    expect = {("CNA00110", 4), ("CNA00120", 3), ("cgg", 1), ("cat", 1), ("JEC21", 1)}
    assert got == expect


def test_count_distinct_gene_name(gb_db):
    expect = {("CNA00110", 1), ("CNA00120", 1)}
    assert {
        tuple(r) for r in gb_db.count_distinct(biotype="gene", name=True).to_list()
    } == expect

    assert {
        tuple(r)
        for r in gb_db.count_distinct(
            seqid="AE017341", biotype="gene", name=True
        ).to_list()
    } == expect


def test_count_distinct_no_match(gb_db):
    # return a table with 0 rows, 2 columns
    got = gb_db.count_distinct(biotype=True, name="blah")
    assert got.shape == (0, 2)
    got = gb_db.count_distinct(biotype="madeup", name=True)
    assert got.shape == (0, 2)


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


def test_gff_get_children_empty(DATA_DIR):
    """if feature has no children then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_children(name="childless"))
    assert got == []


def test_gff_get_parent_empty(DATA_DIR):
    """if feature has no parent then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_parent(name="parentless"))
    assert got == []


def test_gff_get_children_non_existent(DATA_DIR):
    """if feature does not exist then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_children(name="nonexistendID"))
    assert got == []


def test_gff_get_parent_non_existent(DATA_DIR):
    """if feature does not exist then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_parent(name="nonexistendID"))
    assert got == []


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


def test_empty_data():
    _ = GffAnnotationDb()


# testing GenBank files
@pytest.fixture(scope="session")
def gb_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "annotated_seq.gb")


def test_load_annotations_multi(DATA_DIR):
    one = load_annotations(path=DATA_DIR / "simple.gff")
    two = load_annotations(path=DATA_DIR / "simple2.gff")
    expect = len(one) + len(two)
    got = load_annotations(path=DATA_DIR / "simple*.gff")
    assert len(got) == expect


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
    assert not list(seq.get_features())


def test_get_features_matching_no_matching_feature(seq, anno_db):
    """
    Test that `get_features_matching` returns an empty list when there is no matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )

    assert not list(seq.get_features(biotype="exon", name="non_matching"))
    assert not list(seq.get_features(biotype="CDS"))


def test_get_features_matching_matching_feature(seq, anno_db):
    """
    Test that `get_features_matching` returns a list with one matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(1, 4)], strand="+"
    )
    got = list(seq.get_features(biotype="exon"))

    assert got[0].biotype == "exon"
    assert got[0].name == "exon1"
    assert len(got) == 1


def test_feature_strand():
    from cogent3 import make_seq

    #            ++   ++++
    #              ---    --
    raw_seq = "AACCTTTGGGGAATTT"

    plus_spans = [(2, 4), (7, 11)]
    plus_seq = "".join(raw_seq[s:e] for s, e in plus_spans)

    minus_spans = [(4, 7), (11, 13)]
    minus_seq = "".join(raw_seq[s:e] for s, e in minus_spans)
    minus_seq = "".join([{"T": "A", "A": "T"}[b] for b in minus_seq[::-1]])
    seq = make_seq(raw_seq, name="s1", moltype="dna")
    db = GffAnnotationDb()
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
    plus = list(seq.get_features(name="plus"))[0]
    assert str(plus.get_slice()) == plus_seq
    minus = list(seq.get_features(name="minus"))[0]
    assert str(minus.get_slice()) == minus_seq

    # now reverse complement the sequence
    rced = seq.rc()
    plus = list(rced.get_features(name="plus"))[0]
    assert str(plus.get_slice()) == plus_seq
    minus = list(rced.get_features(name="minus"))[0]
    assert str(minus.get_slice()) == minus_seq


def test_feature_nucleic():
    from cogent3 import make_seq
    from cogent3.core import location as loc

    seq = make_seq("AACCTTTGGGGAATTT", moltype="dna")
    mmap = loc.Map(locations=[(4, 7), (11, 13)], parent_length=16)
    expect = seq[mmap.reversed()]

    rcseq = seq.rc()
    rmap = mmap.nucleic_reversed().reversed()
    got = rcseq[rmap]
    assert str(got) == str(expect)


def test_add_feature_with_parent():
    db = GffAnnotationDb()
    db.add_feature(
        seqid="s1",
        biotype="cds",
        name="GG",
        spans=[(0, 100)],
        strand="+",
    )
    db.add_feature(
        seqid="s1",
        biotype="exon",
        name="child",
        spans=[(10, 30)],
        strand="+",
        parent_id="GG",
    )
    got = list(db.get_features_matching(name="GG"))[0]
    child = list(db.get_feature_children(got["name"]))[0]
    assert child["name"] == "child"


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
    got = list(seq.get_features(biotype="exon"))

    assert len(got) == 2


def test_annotate_from_gff(DATA_DIR, seq):
    seq.annotate_from_gff(DATA_DIR / "simple.gff")

    got = list(seq.get_features(biotype="exon"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CDS"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CpG"))
    assert len(got) == 1


def test_get_features_matching_start_stop(DATA_DIR, seq):
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    got = list(seq.get_features(start=2, stop=10, allow_partial=True))
    assert len(got) == 4


def test_matching_conditions():
    got, _ = _matching_conditions({"start": 1, "end": 5}, allow_partial=True)
    expect = "((start >= 1 AND end <= 5) OR (start <= 1 AND end > 1) OR (start < 5 AND end >= 5) OR (start <= 1 AND end >= 5))"
    assert got == expect


def test_matching_conditions_IN():
    got, cond = _matching_conditions(
        {"biotype": ("CDS", "mRNA", "exon")}, allow_partial=True
    )
    expect = "biotype IN (?,?,?)"
    assert got == expect
    assert cond == ("CDS", "mRNA", "exon")


@pytest.mark.parametrize(
    "biotype_value_1", ["CDS", "mRNA", "exon", "three_prime_UTR", "intron"]
)
@pytest.mark.parametrize(
    "biotype_value_2", ["CDS", "mRNA", "exon", "five_prime_UTR", "intron"]
)
def test_get_features_matching_multiple_biotype_tuple(
    seq_db, biotype_value_1, biotype_value_2
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype=(biotype_value_1, biotype_value_2)))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both) and len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


@pytest.mark.parametrize("biotype_value_1", ["CDS", "mRNA", "exon", "intron"])
@pytest.mark.parametrize("biotype_value_2", ["CDS", "mRNA", "exon", "intron"])
def test_get_features_matching_multiple_biotype_list(
    seq_db, biotype_value_1, biotype_value_2
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype=[biotype_value_1, biotype_value_2]))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both) and len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


@pytest.mark.parametrize("biotype_value_1", ["CDS", "mRNA", "exon", "intron"])
@pytest.mark.parametrize("biotype_value_2", ["CDS", "mRNA", "exon", "intron"])
def test_get_features_matching_multiple_biotype_set(
    seq_db, biotype_value_1, biotype_value_2
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype={biotype_value_1, biotype_value_2}))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both) and len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


def test_get_features_matching_start_stop_seqview(DATA_DIR, seq):
    """testing that get_features_matching adjusts"""
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    seq_features = list(seq.get_features(start=0, stop=3, allow_partial=True))
    assert len(seq_features) == 3

    # edge case, only 1 features that overlaps with index 12
    # is actually returning [exon2 at [11:20]/13, CpG1 at [2:12]/13]
    # possibly a bug in the SQL generating code
    subseq = seq[9:]
    seq_features_features = list(
        subseq.get_features(start=3, stop=10, allow_partial=True)
    )
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
    feat = list(seq_db.get_features(name="Transcript:B0019.1"))[0]
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

    feat = list(seq_db.get_features(name="Transcript:B0019.1"))[0]
    rc_seq = seq_db.rc()
    r_feat = list(rc_seq.get_features(name="Transcript:B0019.1"))[0]

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

    feat = list(seq.get_features(name="exon1"))[0]
    r_seq = seq.rc()
    r_feat = list(r_seq.get_features(name="exon1"))[0]

    assert feat.get_slice() == r_feat.get_slice()


def test_add_feature(seq):
    """Sequence supports manual adding of features for a seq with no bound AnnotationDb"""
    record = dict(name="gene-01", biotype="gene", spans=[(12, 16)], strand="+")
    seq.add_feature(**record)
    feats = list(seq.get_features(biotype="gene"))

    assert seq.annotation_db is not None
    assert len(feats) == 1
    assert feats[0].biotype == "gene"


def test_add_feature_existing_db(simple_seq_gff_db):
    """Sequence supports manual adding of features for a seq with an existing AnnotationDb"""
    record = dict(name="gene-01", biotype="gene", spans=[(12, 16)], strand="+")
    simple_seq_gff_db.add_feature(**record)

    # total features should be 5+1=6
    all_feats = list(simple_seq_gff_db.get_features())
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


def test_annotate_from_gff_multiple_calls(DATA_DIR, seq):
    """5 records in each gff file, total features on seq should be 10"""
    seq.annotate_from_gff(DATA_DIR / "simple.gff")
    seq.annotate_from_gff(DATA_DIR / "simple2.gff")
    assert len(list(seq.get_features())) == 10


def test_sequence_collection_annotate_from_gff(DATA_DIR):
    """providing a seqid to SequenceCollection.annotate_from_gff will
    annotate the SequenceCollection, and the Sequence. Both of these will point
    to the same AnnotationDb instance
    """
    seqs = {"test_seq": "ATCGATCGATCG", "test_seq2": "GATCGATCGATC"}
    seq_coll = SequenceCollection(seqs)
    seq_coll.annotate_from_gff(DATA_DIR / "simple.gff", seq_ids="test_seq")

    # the seq for which the seqid was provided is annotated
    seq = seq_coll.get_seq("test_seq")
    assert seq_coll.get_seq("test_seq").annotation_db is not None
    # the annotation_db on the seq and the seq collection are the same object
    assert seq_coll.get_seq("test_seq").annotation_db is seq_coll.annotation_db
    got = list(seq_coll.get_seq("test_seq").get_features(allow_partial=True))
    assert len(got) == 5
    got = list(seq.get_features(biotype="CDS"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CpG"))
    assert len(got) == 1

    # the seq for which the seqid was NOT provided also has a reference to the same db
    seq2 = seq_coll.get_seq("test_seq2")
    assert seq2.annotation_db is not None
    # querying on that sequence returns []
    got = list(seq2.get_features(biotype="CDS"))
    assert not got


def test_seq_coll_query(DATA_DIR):
    """obtain same results when querying from collection as from seq"""
    seqs = {"test_seq": "ATCGATCGATCG", "test_seq2": "GATCGATCGATC"}
    seq_coll = SequenceCollection(seqs)
    seq_coll.annotate_from_gff(DATA_DIR / "simple.gff", seq_ids="test_seq")

    seq = seq_coll.get_seq("test_seq")
    # the seq for which the seqid was provided is annotated
    assert seq.annotation_db is not None
    # todo gah this test fails when allow_partial=False because start / stop
    # not used by seqcoll method
    expect = set(
        (f.seqid, f.biotype, f.name, str(f.map))
        for f in seq.get_features(allow_partial=True)
    )
    got = set(
        (f.seqid, f.biotype, f.name, str(f.map))
        for f in seq_coll.get_features(seqid="test_seq", allow_partial=True)
    )
    assert got == expect
    # the annotation_db on the seq and the seq collection are the same object
    assert seq.annotation_db is seq_coll.annotation_db
    got = list(seq.get_features(biotype="CDS"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CpG"))
    assert len(got) == 1


def test_gff_update_existing(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db.update(gff_small_db)
    assert gff_db.num_matches() == expect


@pytest.mark.parametrize("seqids", (None, "23", ["23"]))
def test_gff_update_existing_specify_seqid(gff_db, gff_small_db, seqids):
    expect = gff_db.num_matches() + gff_small_db.num_matches(
        seqid=seqids[0] if isinstance(seqids, list) else seqids
    )
    gff_db.update(gff_small_db, seqids=seqids)
    assert gff_db.num_matches() == expect


def test_gff_update_db_from_other_db_existing(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db._update_db_from_other_db(other_db=gff_small_db)
    assert gff_db.num_matches() == expect


@pytest.mark.parametrize(
    "seqids", ("test_seq", ["test_seq"], "test_seq2", ["test_seq2"], None)
)
def test_gff_update_db_from_other_db_existing_specify_seqid(
    gff_db, gff_small_db, seqids
):
    expect = gff_db.num_matches() + gff_small_db.num_matches(seqid=seqids)
    gff_db._update_db_from_other_db(other_db=gff_small_db, seqids=seqids)
    assert gff_db.num_matches() == expect


def test_gff_update_db_from_other_db_existing_none_seqid(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db._update_db_from_other_db(other_db=gff_small_db, seqids=None)
    assert gff_db.num_matches() == expect


def test_relative_position_negative_feature(seq_db):
    orig_feat_span = list(seq_db.get_features(name="Transcript:B0019.1"))[0].map

    view = seq_db[5:]
    view_feat_span = list(view.get_features(name="Transcript:B0019.1"))[0].map

    assert orig_feat_span[0].start - 5 == view_feat_span[0].start
    assert orig_feat_span[0].end - 5 == view_feat_span[0].end


def test_relative_position_positive_feature(anno_db):
    seq = DNA.make_seq("AAAAGGGG", name="seq1")

    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name, biotype="exon", name="exon1", spans=[(2, 6)], strand="+"
    )

    orig_feat_span = list(seq.get_features(name="exon1"))[0].map
    view = seq[2:]
    view_feat_span = list(view.get_features(name="exon1"))[0].map

    assert orig_feat_span[0].start - 2 == view_feat_span[0].start
    assert orig_feat_span[0].end - 2 == view_feat_span[0].end


def test_deepcopy(gff_db):
    import copy

    new = copy.deepcopy(gff_db)
    new.add_feature(
        seqid="s1", biotype="exon", name="copied-exon", spans=[(2, 6)], strand="+"
    )
    assert new.num_matches() == gff_db.num_matches() + 1
    assert len(list(new.get_features_matching(name="copied-exon"))) == 1
    assert len(list(gff_db.get_features_matching(name="copied-exon"))) == 0


def test_pickling(gff_db):
    import pickle

    recon = pickle.loads(pickle.dumps(gff_db))
    assert isinstance(recon, type(gff_db))
    recon.add_feature(
        seqid="s1", biotype="exon", name="copied-exon", spans=[(2, 6)], strand="+"
    )
    assert recon.num_matches() == gff_db.num_matches() + 1


@pytest.mark.parametrize("db_name", ("gff_db", "gb_db"))
def test_to_rich_dict(db_name, request):
    db = request.getfixturevalue(db_name)
    data = db.to_rich_dict()
    assert data["init_args"]["source"] == ":memory:"
    if db_name == "gb_db":
        assert "seqid" in data["init_args"]


@pytest.mark.parametrize("db_name", ("gff_db", "gb_db"))
def test_deserialise(db_name, request):
    db = request.getfixturevalue(db_name)
    data = db.to_json()
    got = deserialise.deserialise_object(data)
    assert got is not db
    assert isinstance(got, type(db))
    assert got.num_matches() == db.num_matches()


def test_querying_attributes_gb(gb_db):
    r = list(gb_db.get_records_matching(attributes="lysine biosynthesis"))
    assert "lysine biosynthesis" in r[0]["attributes"]["note"][0]


def test_querying_attributes_gff(gff_db):
    r = list(gff_db.get_records_matching(attributes="amx-2"))
    assert "amx-2" in r[0]["attributes"]


def test_writing_attributes_gff(gff_db):
    gff_db.add_feature(
        biotype="gene",
        seqid="blah",
        name="cancer-gene",
        attributes="description=cancer",
        spans=[(0, 10)],
    )
    r = list(gff_db.get_records_matching(attributes="cancer"))[0]
    assert r["name"] == "cancer-gene"


def test_equal():
    db1 = BasicAnnotationDb()
    db2 = BasicAnnotationDb()
    db3 = BasicAnnotationDb()
    assert db1 != db2
    # we define equality by same class AND same db instance
    db3._db = db2._db
    assert db2 == db3


@pytest.mark.parametrize("other", (GenbankAnnotationDb, GffAnnotationDb))
def test_compatible_symmetric(other):
    basic = BasicAnnotationDb()
    other = other()
    assert basic.compatible(basic)
    assert basic.compatible(other)
    assert other.compatible(other)
    assert other.compatible(basic)


@pytest.mark.parametrize("other", (GenbankAnnotationDb, GffAnnotationDb))
def test_compatible_not_symmetric(other):
    basic = BasicAnnotationDb()
    other = other()
    assert basic.compatible(basic, symmetric=False)
    assert not basic.compatible(other, symmetric=False)
    assert other.compatible(other, symmetric=False)
    assert other.compatible(basic, symmetric=False)


def test_incompatible():
    gff = GffAnnotationDb()
    gb = GenbankAnnotationDb()
    assert not gff.compatible(gb)
    assert not gb.compatible(gff)


@pytest.mark.parametrize("wrong_type", ({}, BasicAnnotationDb().db))
def test_incompatible_invalid_type(wrong_type):
    db = BasicAnnotationDb()
    with pytest.raises(TypeError):
        db.compatible(wrong_type)


def _custom_namer(data):
    for key in ("gene", "locus_tag", "strain"):
        if key in data:
            return data[key]
    return ["default name"]


def test_gb_namer(DATA_DIR):
    path = DATA_DIR / "annotated_seq.gb"
    got = list(MinimalGenbankParser(path.read_text().splitlines()))
    data = got[0]["features"]
    db = GenbankAnnotationDb(data=data, namer=_custom_namer, seqid=got[0]["locus"])
    # there are 2 repeat regions, which we don't catch with our namer
    assert db.num_matches(name="default name") == 2


def test_write(gb_db, tmp_path):
    outpath = tmp_path / "ondisk.gbkdb"
    gb_db.write(outpath)
    got = GenbankAnnotationDb(source=outpath)
    assert got.to_rich_dict()["tables"] == gb_db.to_rich_dict()["tables"]
    assert isinstance(got, GenbankAnnotationDb)


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("annotations")


def test_load_anns_with_write(DATA_DIR, tmp_dir):
    inpath = DATA_DIR / "simple.gff"
    outpath = tmp_dir / "simple.gffdb"
    orig = load_annotations(path=inpath, write_path=outpath)
    orig.db.close()
    expect = load_annotations(path=inpath)
    got = GffAnnotationDb(source=outpath)
    assert len(got) == len(expect)
    got_data = got.to_rich_dict()
    expect_data = expect.to_rich_dict()
    assert got_data["tables"] == expect_data["tables"]


def test_gbdb_get_children_fails_no_coords(gb_db):
    with pytest.raises(ValueError):
        _ = list(gb_db.get_feature_children(name="CNA00110"))


def test_gbdb_get_parent_fails_no_coords(gb_db):
    with pytest.raises(ValueError):
        _ = list(gb_db.get_feature_parent(name="CNA00110"))

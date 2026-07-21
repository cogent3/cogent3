import contextlib
import io
import pathlib

import numpy
import pytest
from scinexus import deserialise

import cogent3
from cogent3.core.annotation_db import (
    SCHEMA_VERSION,
    AnnotationDbABC,
    BasicAnnotationDb,
    GenbankAnnotationDb,
    GffAnnotationDb,
    LookupTableCache,
    _make_db_connection,
    _matching_conditions,
    _rename_column_if_exists,
    load_annotations,
    upgrade_annotation_db,
)
from cogent3.parse import genbank

DNA = cogent3.get_moltype("dna")


def close_dbs(*objs):
    for obj in objs:
        if not hasattr(obj, "has_annotation_db"):
            db = obj
        elif obj.has_annotation_db():
            db = obj.annotation_db
        else:
            continue

        db.close()


@pytest.fixture
def gff_db(DATA_DIR):
    path = DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"
    db = load_annotations(path=path)
    yield db
    db.close()


@pytest.fixture
def gff_small_db(DATA_DIR):
    path = DATA_DIR / "simple.gff"
    db = load_annotations(path=path)
    yield db
    db.close()


@pytest.fixture
def seq_db(DATA_DIR):
    seq = cogent3.load_seq(
        DATA_DIR / "c_elegans_WS199_dna_shortened.fasta",
        moltype="dna",
    )
    db = load_annotations(path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3")

    seq.annotation_db = db

    yield seq
    db.close()


@pytest.fixture
def seq():
    seq = cogent3.make_seq("ATTGTACGCCTTTTTTATTATT", name="test_seq", moltype="dna")
    yield seq
    if seq.has_annotation_db():
        with contextlib.suppress(AttributeError):
            seq.annotation_db.close()


@pytest.fixture
def anno_db() -> BasicAnnotationDb:
    # an empty db that we can add to
    db = BasicAnnotationDb()
    yield db
    db.close()


@pytest.fixture
def simple_seq_gff_db(DATA_DIR):
    seq = cogent3.make_seq("ATTGTACGCCTTTTTTATTATT", name="test_seq", moltype="dna")
    seq.annotation_db = load_annotations(path=DATA_DIR / "simple.gff")
    yield seq
    seq.annotation_db.close()


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
    ("db_name", "cls"),
    [("gff_db", GenbankAnnotationDb), ("gb_db", GffAnnotationDb)],
)
def test_constructor_db_fail(db_name, cls, request):
    db = request.getfixturevalue(db_name)
    with pytest.raises(TypeError):
        cls(db=db)


@pytest.mark.parametrize(
    ("db_name", "cls"),
    [("gff_db", GenbankAnnotationDb), ("gb_db", GffAnnotationDb)],
)
def test_constructor_wrong_db_schema(db_name, cls, request):
    db = request.getfixturevalue(db_name)
    with pytest.raises(TypeError):
        cls(db=db.db)


@pytest.mark.parametrize(
    ("db_name", "cls"),
    [
        ("gff_db", GffAnnotationDb),
        ("anno_db", GffAnnotationDb),
        ("gb_db", GenbankAnnotationDb),
        ("anno_db", GenbankAnnotationDb),
    ],
)
def test_constructor_db_instance_works(db_name, cls, request):
    # only compatible db's used to init
    db = request.getfixturevalue(db_name)
    n = cls(db=db)
    n.close()


@pytest.mark.parametrize(
    ("db_name", "cls"),
    [
        ("gff_db", GffAnnotationDb),
        ("anno_db", GffAnnotationDb),
        ("gb_db", GenbankAnnotationDb),
        ("anno_db", GenbankAnnotationDb),
    ],
)
def test_constructor_db_connection_works(db_name, cls, request):
    # only compatible db's used to init
    db = request.getfixturevalue(db_name)
    cls(db=db.db)


def test_gff_describe(gff_db):
    from cogent3.core.table import Table

    result = gff_db.describe
    assert isinstance(result, Table)


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
            seqid="AE017341",
            biotype="gene",
            name=True,
        ).to_list()
    } == expect


def test_count_distinct_no_match(gb_db):
    # return a table with 0 rows, 2 columns
    got = gb_db.count_distinct(biotype=True, name="blah")
    assert got.shape == (0, 2)
    got = gb_db.count_distinct(biotype="madeup", name=True)
    assert got.shape == (0, 2)


@pytest.mark.parametrize(
    ("db_fixture", "db_class"),
    [
        ("gb_db", GenbankAnnotationDb),
        ("gff_db", GffAnnotationDb),
        ("anno_db", BasicAnnotationDb),
    ],
)
def test_count_distinct_after_db_reuse(db_fixture, db_class, request):
    """count_distinct works when db class reuses existing db connection.

    Regression test: _setup_db was not copying _schema_version and
    _lookup_cache when binding to an existing db of the same class,
    causing "no such column: biotype" errors.
    """
    first_db = request.getfixturevalue(db_fixture)
    second_db = db_class(db=first_db)
    result = second_db.count_distinct(biotype=True)
    assert result is not None
    second_db.close()


def test_gff_features_matching(gff_db):
    result = list(gff_db.get_features_matching(biotype="CDS"))
    assert result


def test_gff_get_children(gff_db):
    # an ID with 8 children
    children = list(gff_db.get_feature_children(name="Transcript:B0019.1"))
    assert len(children) == 8
    # this is parent of the transcript, and due to the funky structure of
    # this gff, also to one of the CDS entries
    children = list(gff_db.get_feature_children(name="Gene:WBGene00000138"))
    assert len(children) == 2


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("CDS:B0019.1", ("Transcript:B0019.1", "Gene:WBGene00000138")),
        ("Transcript:B0019.1", ("Gene:WBGene00000138",)),
    ],
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
    db.close()


def test_gff_get_parent_empty(DATA_DIR):
    """if feature has no parent then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_parent(name="parentless"))
    assert got == []
    db.close()


def test_gff_get_children_non_existent(DATA_DIR):
    """if feature does not exist then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_children(name="nonexistendID"))
    assert got == []
    db.close()


def test_gff_get_parent_non_existent(DATA_DIR):
    """if feature does not exist then should return []"""
    db = load_annotations(path=DATA_DIR / "simple2.gff")
    got = list(db.get_feature_parent(name="nonexistendID"))
    assert got == []
    db.close()


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
    record = {
        "seqid": "2",
        "name": "gene-01",
        "biotype": "gene",
        "spans": [(23, 33)],
        "strand": "+",
    }
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
    db = GffAnnotationDb()
    db.close()


# testing GenBank files
@pytest.fixture
def gb_db(DATA_DIR):
    db = load_annotations(path=DATA_DIR / "annotated_seq.gb")
    yield db
    db.close()


def test_load_annotations_multi(DATA_DIR):
    one = load_annotations(path=DATA_DIR / "simple.gff")
    two = load_annotations(path=DATA_DIR / "simple2.gff")
    expect = len(one) + len(two)
    got = load_annotations(path=DATA_DIR / "simple*.gff")
    assert len(got) == expect
    close_dbs(one, two, got)


@pytest.fixture(params=["annotated_seq.gb", "simple2.gff"])
def compressed_flat_file(DATA_DIR, tmp_path, request):
    src = DATA_DIR / request.param
    out_path = tmp_path / f"{request.param}.gz"
    with cogent3.open_(out_path, mode="w") as out:
        out.write(src.read_bytes())
    return out_path


def test_load_annotations_compressed(compressed_flat_file):
    got = load_annotations(path=compressed_flat_file)
    assert isinstance(got, AnnotationDbABC)
    got.close()


def test_load_annotations_chunked(gff_db, DATA_DIR):
    path = DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"

    name = "CDS:B0019.1"
    expect = next(iter(gff_db.get_features_matching(name=name)))
    assert len(expect["spans"]) == 3
    # two lines splits a 3 line record into 2 and 1 line, so the
    # update record code is invoked
    db = load_annotations(path=path, write_path=":memory:", lines_per_block=2)
    got = next(iter(db.get_features_matching(name=name)))
    assert got.pop("spans") == expect.pop("spans")
    assert got == expect
    db.close()


@pytest.mark.parametrize(("parent_biotype", "name"), [("gene", "CNA00110")])
def test_gb_get_children(gb_db, parent_biotype, name):
    parent = next(iter(gb_db.get_features_matching(biotype=parent_biotype, name=name)))
    coords = numpy.array(parent["spans"])
    child = next(
        iter(
            gb_db.get_feature_children(
                name=name,
                exclude_biotype=parent_biotype,
                start=coords.min(),
                stop=coords.max(),
            ),
        ),
    )
    assert child["biotype"] != parent["biotype"]
    assert child["name"] == parent["name"]


def test_gb_get_parent(gb_db):
    cds_id = "CNA00110"
    cds = next(iter(gb_db.get_features_matching(biotype="CDS", name=cds_id)))
    coords = numpy.array(cds["spans"])
    parent = next(
        iter(
            gb_db.get_feature_parent(
                name=cds_id,
                exclude_biotype="CDS",
                start=coords.min(),
                stop=coords.max(),
            ),
        ),
    )
    assert parent["biotype"] != cds["biotype"]
    assert parent["biotype"] == "gene"
    assert parent["name"] == cds["name"]


def test_protocol_adherence(gff_db, gb_db):
    for db in (gff_db, gb_db):
        assert isinstance(db, AnnotationDbABC)


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
        seqid=seq.name,
        biotype="exon",
        name="exon1",
        spans=[(1, 4)],
        strand="+",
    )

    assert not list(seq.get_features(biotype="exon", name="non_matching"))
    assert not list(seq.get_features(biotype="CDS"))


def test_get_features_matching_matching_feature(seq, anno_db):
    """
    Test that `get_features_matching` returns a list with one matching feature in the annotation database.
    """
    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name,
        biotype="exon",
        name="exon1",
        spans=[(1, 4)],
        strand="+",
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
    seq = make_seq(seq=raw_seq, name="s1", moltype="dna")
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
    plus = next(iter(seq.get_features(name="plus")))
    assert str(plus.get_slice()) == plus_seq
    minus = next(iter(seq.get_features(name="minus")))
    assert str(minus.get_slice()) == minus_seq

    # now reverse complement the sequence
    rced = seq.rc()
    plus = next(iter(rced.get_features(name="plus")))
    assert str(plus.get_slice()) == plus_seq
    minus = next(iter(rced.get_features(name="minus")))
    assert str(minus.get_slice()) == minus_seq
    db.close()


def test_feature_nucleic():
    from cogent3 import make_seq
    from cogent3.core import location as loc

    #                              111111
    #                    0123456789012345
    seq = make_seq(seq="AACCTTTGGGGAATTT", moltype="dna")
    mmap = loc.FeatureMap.from_locations(
        locations=[(4, 7), (10, 12)],
        parent_length=len(seq),
    )
    expect = seq[mmap]

    rcseq = seq.rc()
    rmap = mmap.nucleic_reversed()
    got = rcseq[rmap].rc()
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
    got = next(iter(db.get_features_matching(name="GG")))
    child = next(iter(db.get_feature_children(got["name"])))
    assert child["name"] == "child"
    db.close()


def test_get_features_matching_matching_features(anno_db: GffAnnotationDb, seq):
    """
    Test that `get_features_matching` returns a list with all matching features in the annotation database.
    """
    seq.annotation_db = anno_db

    anno_db.add_feature(
        seqid=seq.name,
        biotype="exon",
        name="exon1",
        spans=[(1, 4)],
        strand="+",
    )
    anno_db.add_feature(
        seqid=seq.name,
        biotype="exon",
        name="exon2",
        spans=[(6, 10)],
        strand="+",
    )
    got = list(seq.get_features(biotype="exon"))

    assert len(got) == 2


def test_annotate_from_gff(DATA_DIR, seq):
    seq.annotation_db = load_annotations(path=DATA_DIR / "simple.gff")

    got = list(seq.get_features(biotype="exon"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CDS"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CpG"))
    assert len(got) == 1


def test_get_features_matching_start_stop(DATA_DIR, seq):
    seq.annotation_db = load_annotations(path=DATA_DIR / "simple.gff")
    got = list(seq.get_features(start=2, stop=10, allow_partial=True))
    assert len(got) == 4


def test_matching_conditions():
    got, _ = _matching_conditions({"start": 1, "stop": 5}, allow_partial=True)
    expect = "((start >= 1 AND stop <= 5) OR (start <= 1 AND stop > 1) OR (start < 5 AND stop >= 5) OR (start <= 1 AND stop >= 5))"
    assert got == expect


def test_matching_conditions_IN():
    got, cond = _matching_conditions(
        {"biotype": ("CDS", "mRNA", "exon")},
        allow_partial=True,
    )
    expect = "biotype IN (?,?,?)"
    assert got == expect
    assert cond == ("CDS", "mRNA", "exon")


def test_matching_conditions_stop_only():
    got, _ = _matching_conditions({"stop": 20000}, allow_partial=False)
    assert got == "(stop <= 20000)"


def test_matching_conditions_stop_only_partial():
    got, _ = _matching_conditions({"stop": 20000}, allow_partial=True)
    assert got == "(start < 20000)"


def test_matching_conditions_start_only():
    got, _ = _matching_conditions({"start": 500}, allow_partial=False)
    assert got == "(start >= 500)"


def test_matching_conditions_start_only_partial():
    got, _ = _matching_conditions({"start": 500}, allow_partial=True)
    assert got == "(stop > 500)"


def test_get_features_stop_only():
    """querying with stop-only should return features in [0, stop), not
    just features that straddle the stop position"""
    from cogent3.core.annotation_db import GenbankAnnotationDb

    db = GenbankAnnotationDb()
    db.add_feature(seqid="seq1", biotype="CDS", name="a", spans=[(100, 500)])
    db.add_feature(seqid="seq1", biotype="CDS", name="b", spans=[(1000, 5000)])
    db.add_feature(seqid="seq1", biotype="CDS", name="c", spans=[(10000, 15000)])
    db.add_feature(seqid="seq1", biotype="CDS", name="d", spans=[(19000, 25000)])

    # stop=20000, allow_partial=False: features fully within [0, 20000)
    results = list(db.get_features_matching(seqid="seq1", biotype="CDS", stop=20000))
    names = {r["name"] for r in results}
    assert names == {"a", "b", "c"}

    # stop=20000, allow_partial=True: features overlapping [0, 20000)
    results = list(
        db.get_features_matching(
            seqid="seq1",
            biotype="CDS",
            stop=20000,
            allow_partial=True,
        ),
    )
    names = {r["name"] for r in results}
    assert names == {"a", "b", "c", "d"}


@pytest.mark.parametrize(
    "biotype_value_1",
    ["CDS", "mRNA", "exon", "three_prime_UTR", "intron"],
)
@pytest.mark.parametrize(
    "biotype_value_2",
    ["CDS", "mRNA", "exon", "five_prime_UTR", "intron"],
)
def test_get_features_matching_multiple_biotype_tuple(
    seq_db,
    biotype_value_1,
    biotype_value_2,
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype=(biotype_value_1, biotype_value_2)))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both)
        assert len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


@pytest.mark.parametrize("biotype_value_1", ["CDS", "mRNA", "exon", "intron"])
@pytest.mark.parametrize("biotype_value_2", ["CDS", "mRNA", "exon", "intron"])
def test_get_features_matching_multiple_biotype_list(
    seq_db,
    biotype_value_1,
    biotype_value_2,
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype=[biotype_value_1, biotype_value_2]))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both)
        assert len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


@pytest.mark.parametrize("biotype_value_1", ["CDS", "mRNA", "exon", "intron"])
@pytest.mark.parametrize("biotype_value_2", ["CDS", "mRNA", "exon", "intron"])
def test_get_features_matching_multiple_biotype_set(
    seq_db,
    biotype_value_1,
    biotype_value_2,
):
    """querying for features with multiple values should return the
    same result as the sum of querying for each value seperately"""
    where_1 = list(seq_db.get_features(biotype=biotype_value_1))
    where_2 = list(seq_db.get_features(biotype=biotype_value_2))
    in_both = list(seq_db.get_features(biotype={biotype_value_1, biotype_value_2}))

    if biotype_value_1 == biotype_value_2:
        assert len(where_1) == len(in_both)
        assert len(where_2) == len(in_both)
    else:
        assert len(where_1) + len(where_2) == len(in_both)


def test_get_features_matching_start_stop_seqview(DATA_DIR, seq):
    """testing that get_features_matching adjusts"""
    db = load_annotations(path=DATA_DIR / "simple.gff")
    seq.annotation_db = db
    seq_features = list(seq.get_features(start=0, stop=3, allow_partial=True))
    assert len(seq_features) == 3

    # edge case, only 1 features that overlaps with index 12
    # is actually returning [exon2 at [11:20]/13, CpG1 at [2:12]/13]
    # possibly a bug in the SQL generating code
    subseq = seq[9:]
    seq_features_features = list(
        subseq.get_features(start=3, stop=10, allow_partial=True),
    )
    assert len(seq_features_features) == 1
    close_dbs(db)


@pytest.fixture
def populated_basic_db() -> BasicAnnotationDb:
    db = BasicAnnotationDb()
    for i in range(5):
        db.add_feature(
            seqid="seq1",
            biotype="exon",
            name=f"exon{i}",
            spans=[(i * 10, i * 10 + 5)],
            strand="+",
        )
    yield db
    db.close()


@pytest.mark.parametrize("limit", [1, 2, 5, 10])
def test_get_features_matching_limit_basic(populated_basic_db, limit):
    got = list(populated_basic_db.get_features_matching(biotype="exon", limit=limit))
    assert len(got) == min(limit, 5)


@pytest.mark.parametrize("limit", [1, 2, 5, 10])
def test_get_records_matching_limit_basic(populated_basic_db, limit):
    got = list(populated_basic_db.get_records_matching(biotype="exon", limit=limit))
    assert len(got) == min(limit, 5)


def test_get_features_matching_limit_none_equivalent_to_unlimited(populated_basic_db):
    default = list(populated_basic_db.get_features_matching(biotype="exon"))
    with_none = list(
        populated_basic_db.get_features_matching(biotype="exon", limit=None)
    )
    assert len(default) == len(with_none) == 5


@pytest.mark.parametrize("limit", [0, -1, -5])
def test_get_features_matching_limit_invalid(populated_basic_db, limit):
    with pytest.raises(ValueError):
        # generator must be consumed for the exception to fire
        list(populated_basic_db.get_features_matching(limit=limit))


@pytest.mark.parametrize("limit", [0, -1])
def test_get_records_matching_limit_invalid(populated_basic_db, limit):
    with pytest.raises(ValueError):
        list(populated_basic_db.get_records_matching(limit=limit))


def test_get_features_matching_limit_with_filters(populated_basic_db):
    got = list(
        populated_basic_db.get_features_matching(
            biotype="exon",
            start=0,
            stop=25,
            allow_partial=True,
            limit=2,
        )
    )
    assert len(got) == 2
    # confirm filter still applies (the unlimited result would be 3)
    full = list(
        populated_basic_db.get_features_matching(
            biotype="exon",
            start=0,
            stop=25,
            allow_partial=True,
        )
    )
    assert len(full) == 3


@pytest.fixture
def gff_db_with_user_exons(DATA_DIR):
    # simple.gff loads 2 exon rows into the 'gff' table; add 2 more into
    # the 'user' table to exercise limit-spanning behaviour across tables
    db = load_annotations(path=DATA_DIR / "simple.gff")
    for i in range(2):
        db.add_feature(
            seqid="test_seq",
            biotype="exon",
            name=f"user_exon{i}",
            spans=[(100 + i * 10, 105 + i * 10)],
            strand="+",
        )
    yield db
    db.close()


def test_get_features_matching_spans_tables_total(gff_db_with_user_exons):
    got = list(gff_db_with_user_exons.get_features_matching(biotype="exon"))
    assert len(got) == 4


@pytest.mark.parametrize(("limit", "expected"), [(1, 1), (2, 2), (3, 3), (100, 4)])
def test_get_features_matching_limit_spans_tables(
    gff_db_with_user_exons,
    limit,
    expected,
):
    got = list(
        gff_db_with_user_exons.get_features_matching(biotype="exon", limit=limit)
    )
    assert len(got) == expected


def test_get_records_matching_spans_tables_total(gff_db_with_user_exons):
    total = sum(1 for _ in gff_db_with_user_exons.get_records_matching(biotype="exon"))
    assert total == 4


def test_get_records_matching_limit_spans_tables(gff_db_with_user_exons):
    got = list(gff_db_with_user_exons.get_records_matching(biotype="exon", limit=2))
    assert len(got) == 2


@pytest.fixture
def seq_with_exons(seq, anno_db):
    seq.annotation_db = anno_db
    for i in range(4):
        anno_db.add_feature(
            seqid=seq.name,
            biotype="exon",
            name=f"exon{i}",
            spans=[(i * 3, i * 3 + 2)],
            strand="+",
        )
    return seq


def test_sequence_get_features_unlimited(seq_with_exons):
    full = list(seq_with_exons.get_features(biotype="exon"))
    assert len(full) == 4


def test_sequence_get_features_limit(seq_with_exons):
    got = list(seq_with_exons.get_features(biotype="exon", limit=2))
    assert len(got) == 2


def test_get_slice():
    """get_slice should return the same as slicing the sequence directly"""
    seq = cogent3.make_seq("ATTGTACGCCCCTGA", name="test_seq", moltype="dna")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": "+",
    }
    feature = seq.make_feature(feature_data)
    got = feature.get_slice()
    assert str(got) == str(seq[5:10])


def test_get_slice_annotation_offset():
    """get_slice should return the same as slicing the sequence directly"""
    seq = cogent3.make_seq("ATTGTACGCCCCTGA", name="test_seq", moltype="dna")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 10),
        ],
        "strand": "+",
    }
    feature = seq.make_feature(feature_data)
    assert feature.map.num_spans == 1
    got = feature.get_slice()
    assert got.annotation_offset == 5


def test_get_slice_annotation_offset_not_set():
    # annotation offset is set to zero if a feature spans disjoint segments
    # plus the annotation db is set to None
    seq = cogent3.make_seq("ATTGTACGCCCCTGA", name="test_seq", moltype="dna")
    feature_data = {
        "biotype": "CDS",
        "name": "fake",
        "spans": [
            (5, 7),
            (9, 11),
        ],
        "strand": "+",
    }
    feature = seq.make_feature(feature_data)
    assert feature.map.num_spans == 2
    got = feature.get_slice()
    assert got.annotation_offset == 0
    assert not len(got.annotation_db)
    f = list(got.get_features(biotype="gene"))
    assert not f
    close_dbs(got)


def test_feature_get_children(seq_db):
    feat = next(iter(seq_db.get_features(name="Transcript:B0019.1")))
    new_feat_5pUTR = list(feat.get_children(biotype="five_prime_UTR"))
    assert len(new_feat_5pUTR) == 1

    new_feat_CDS = next(iter(feat.get_children(biotype="CDS")))
    assert new_feat_CDS.name == "CDS:B0019.1"


def test_db_persists_post_rc(seq_db):
    """assert that the db persists after the .rc() method call"""
    rc_seq = seq_db.rc()
    assert rc_seq.annotation_db is not None


def test_rc_get_slice_negative_feature(seq_db):
    """given a feature on the - strand, the feature.get_slice() should return
    the same sequence before and after the sequence is reverse complemented
    """

    feat = next(iter(seq_db.get_features(name="Transcript:B0019.1")))
    rc_seq = seq_db.rc()
    r_feat = next(iter(rc_seq.get_features(name="Transcript:B0019.1")))

    assert feat.get_slice() == r_feat.get_slice()


def test_rc_get_slice_positive_feature(anno_db):
    """given a feature on the + strand, the feature.get_slice() should return
    the same sequence before and after the sequence is reverse complemented
    """

    seq = DNA.make_seq(seq="AAAAGGGG", name="seq1")

    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name,
        biotype="exon",
        name="exon1",
        spans=[(2, 6)],
        strand="+",
    )

    feat = next(iter(seq.get_features(name="exon1")))
    r_seq = seq.rc()
    r_feat = next(iter(r_seq.get_features(name="exon1")))

    assert feat.get_slice() == r_feat.get_slice()


def test_add_feature(seq):
    """Sequence supports manual adding of features for a seq with no bound AnnotationDb"""
    record = {"name": "gene-01", "biotype": "gene", "spans": [(12, 16)], "strand": "+"}
    seq.add_feature(**record)
    feats = list(seq.get_features(biotype="gene"))

    assert seq.annotation_db is not None
    assert len(feats) == 1
    assert feats[0].biotype == "gene"


def test_add_feature_existing_db(simple_seq_gff_db):
    """Sequence supports manual adding of features for a seq with an existing AnnotationDb"""
    record = {"name": "gene-01", "biotype": "gene", "spans": [(12, 16)], "strand": "+"}
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
    assert seq_sliced._seq.parent == str(simple_seq_gff_db)
    # check the annotation_db is still attached and the same instance
    assert seq_sliced.annotation_db
    assert seq_sliced.annotation_db is simple_seq_gff_db.annotation_db


def test_sequence_collection_annotate_from_gff(DATA_DIR):
    """providing a seqid to SequenceCollection.annotate_from_gff will
    annotate the SequenceCollection, and the Sequence. Both of these will point
    to the same AnnotationDb instance
    """
    seqs = {"test_seq": "ATCGATCGATCG", "test_seq2": "GATCGATCGATC"}
    seq_coll = cogent3.make_unaligned_seqs(seqs, moltype="dna")
    db = cogent3.load_annotations(path=DATA_DIR / "simple.gff")
    seq_coll.annotation_db = db

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
    # querying on that sequence returns []
    got = list(seq2.get_features(biotype="CDS"))
    assert not got
    db.close()


def test_seq_coll_query(DATA_DIR):
    """obtain same results when querying from collection as from seq"""
    seqs = {"test_seq": "ATCGATCGATCG", "test_seq2": "GATCGATCGATC"}
    seq_coll = cogent3.make_unaligned_seqs(seqs, moltype="dna")
    db = cogent3.load_annotations(path=DATA_DIR / "simple.gff")
    seq_coll.annotation_db = db

    seq = seq_coll.get_seq("test_seq")
    # the seq for which the seqid was provided is annotated
    assert seq.has_annotation_db()
    # TODO gah this test fails when allow_partial=False because start / stop
    # not used by seqcoll method
    expect = {
        (f.seqid, f.biotype, f.name, str(f.map))
        for f in seq.get_features(allow_partial=True)
    }
    got = {
        (f.seqid, f.biotype, f.name, str(f.map))
        for f in seq_coll.get_features(seqid="test_seq", allow_partial=True)
    }
    assert got == expect
    # the annotation_db on the seq and the seq collection are the same object
    assert seq.annotation_db is seq_coll.annotation_db
    got = list(seq.get_features(biotype="CDS"))
    assert len(got) == 2

    got = list(seq.get_features(biotype="CpG"))
    assert len(got) == 1
    db.close()


def test_gff_update_existing(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db.update(gff_small_db)
    assert gff_db.num_matches() == expect


@pytest.mark.parametrize("seqids", [None, "23", ["23"]])
def test_gff_update_existing_specify_seqid(gff_db, gff_small_db, seqids):
    expect = gff_db.num_matches() + gff_small_db.num_matches(
        seqid=seqids[0] if isinstance(seqids, list) else seqids,
    )
    gff_db.update(gff_small_db, seqids=seqids)
    assert gff_db.num_matches() == expect


def test_gff_update_db_from_other_db_existing(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db._update_db_from_other_db(other_db=gff_small_db)
    assert gff_db.num_matches() == expect


@pytest.mark.parametrize(
    "seqids",
    ["test_seq", ["test_seq"], "test_seq2", ["test_seq2"], None],
)
def test_gff_update_db_from_other_db_existing_specify_seqid(
    gff_db,
    gff_small_db,
    seqids,
):
    expect = gff_db.num_matches() + gff_small_db.num_matches(seqid=seqids)
    gff_db._update_db_from_other_db(other_db=gff_small_db, seqids=seqids)
    assert gff_db.num_matches() == expect


def test_gff_update_db_from_other_db_existing_none_seqid(gff_db, gff_small_db):
    expect = gff_db.num_matches() + gff_small_db.num_matches()
    gff_db._update_db_from_other_db(other_db=gff_small_db, seqids=None)
    assert gff_db.num_matches() == expect


def test_relative_position_negative_feature(seq_db):
    orig_feat_span = next(iter(seq_db.get_features(name="Transcript:B0019.1"))).map

    view = seq_db[5:]
    view_feat_span = next(iter(view.get_features(name="Transcript:B0019.1"))).map

    assert orig_feat_span[0].start - 5 == view_feat_span[0].start
    assert orig_feat_span[0].end - 5 == view_feat_span[0].end


def test_relative_position_positive_feature(anno_db):
    seq = DNA.make_seq(seq="AAAAGGGG", name="seq1")

    seq.annotation_db = anno_db
    anno_db.add_feature(
        seqid=seq.name,
        biotype="exon",
        name="exon1",
        spans=[(2, 6)],
        strand="+",
    )

    orig_feat_span = next(iter(seq.get_features(name="exon1"))).map
    view = seq[2:]
    view_feat_span = next(iter(view.get_features(name="exon1"))).map

    assert orig_feat_span[0].start - 2 == view_feat_span[0].start
    assert orig_feat_span[0].end - 2 == view_feat_span[0].end


def test_deepcopy(gff_db):
    import copy

    new = copy.deepcopy(gff_db)
    new.add_feature(
        seqid="s1",
        biotype="exon",
        name="copied-exon",
        spans=[(2, 6)],
        strand="+",
    )
    assert new.num_matches() == gff_db.num_matches() + 1
    assert len(list(new.get_features_matching(name="copied-exon"))) == 1
    assert len(list(gff_db.get_features_matching(name="copied-exon"))) == 0
    close_dbs(new)


def test_pickling(gff_db):
    import pickle

    recon = pickle.loads(pickle.dumps(gff_db))
    assert isinstance(recon, type(gff_db))
    recon.add_feature(
        seqid="s1",
        biotype="exon",
        name="copied-exon",
        spans=[(2, 6)],
        strand="+",
    )
    assert recon.num_matches() == gff_db.num_matches() + 1
    close_dbs(recon)


@pytest.mark.parametrize("db_name", ["gff_db", "gb_db"])
def test_to_rich_dict(db_name, request):
    db = request.getfixturevalue(db_name)
    data = db.to_rich_dict()
    assert data["init_args"]["source"] == ":memory:"
    if db_name == "gb_db":
        assert "seqid" in data["init_args"]


@pytest.mark.parametrize("db_name", ["gff_db", "gb_db"])
def test_deserialise(db_name, request):
    db = request.getfixturevalue(db_name)
    data = db.to_json()
    got = deserialise.deserialise_object(data)
    assert got is not db
    assert isinstance(got, type(db))
    assert got.num_matches() == db.num_matches()
    close_dbs(got, db)


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
        name="agene",
        attributes="description=cancer",
        spans=[(0, 10)],
    )
    r = next(iter(gff_db.get_records_matching(attributes="cancer")))
    assert r["name"] == "agene"


def test_equal():
    db1 = BasicAnnotationDb()
    db2 = BasicAnnotationDb()
    db3 = BasicAnnotationDb()
    assert db1 != db2
    # we define equality by same class AND same db instance, so we close the db
    db3.close()
    # then override a private attribute for this test
    db3._db = db2._db
    assert db2 == db3
    close_dbs(db1, db2, db3)


@pytest.mark.parametrize("other", [GenbankAnnotationDb, GffAnnotationDb])
def test_compatible_symmetric(other):
    basic = BasicAnnotationDb()
    other = other()
    assert basic.compatible(basic)
    assert basic.compatible(other)
    assert other.compatible(other)
    assert other.compatible(basic)
    close_dbs(basic, other)


@pytest.mark.parametrize("other", [GenbankAnnotationDb, GffAnnotationDb])
def test_compatible_not_symmetric(other):
    basic = BasicAnnotationDb()
    other = other()
    assert basic.compatible(basic, symmetric=False)
    assert not basic.compatible(other, symmetric=False)
    assert other.compatible(other, symmetric=False)
    assert other.compatible(basic, symmetric=False)
    close_dbs(basic, other)


def test_incompatible():
    gff = GffAnnotationDb()
    gb = GenbankAnnotationDb()
    assert not gff.compatible(gb)
    assert not gb.compatible(gff)
    close_dbs(gff, gb)


@pytest.mark.parametrize("wrong_type", [{}, BasicAnnotationDb().db])
def test_incompatible_invalid_type(wrong_type):
    db = BasicAnnotationDb()
    with pytest.raises(TypeError):
        db.compatible(wrong_type)
    db.close()


def _custom_namer(data):
    return next(
        (data[key] for key in ("gene", "locus_tag", "strain") if key in data),
        ["default name"],
    )


def test_gb_namer(DATA_DIR):
    path = DATA_DIR / "annotated_seq.gb"
    got = list(genbank.minimal_parser(path))
    data = got[0]["features"]
    db = GenbankAnnotationDb(data=data, namer=_custom_namer, seqid=got[0]["locus"])
    # there are 2 repeat regions, which we don't catch with our namer
    assert db.num_matches(name="default name") == 2
    close_dbs(db)


def test_write(gb_db, tmp_path):
    outpath = tmp_path / "ondisk.c3gbdb"
    gb_db.write(outpath)
    got = GenbankAnnotationDb.from_file(outpath)
    assert got.to_rich_dict()["tables"] == gb_db.to_rich_dict()["tables"]
    assert isinstance(got, GenbankAnnotationDb)
    close_dbs(got, gb_db)


@pytest.mark.parametrize("transform", [str, pathlib.Path])
def test_write_home_dir(gb_db, HOME_TMP_DIR, transform):
    # write to a "~/" home directory style path
    gb_db.write(transform(f"~/{HOME_TMP_DIR.name}/ondisk.c3gbdb"))
    got = GenbankAnnotationDb.from_file(HOME_TMP_DIR / "ondisk.c3gbdb")
    assert got.to_rich_dict()["tables"] == gb_db.to_rich_dict()["tables"]
    assert isinstance(got, GenbankAnnotationDb)
    close_dbs(got, gb_db)


def convert_to_old_np_format(data: bytes) -> bytes:
    with io.BytesIO(data) as stream:
        return numpy.load(stream).tobytes()


def convert_spans_column(db, table_name):
    # Convert the database to the old numpy format.
    conn = db.db
    conn.create_function("convert_to_old_np_format", 1, convert_to_old_np_format)
    cursor = conn.cursor()
    cursor.execute(f"UPDATE {table_name} SET spans = convert_to_old_np_format(spans);")
    conn.commit()


@pytest.mark.parametrize(
    ("db_name", "cls"),
    [("gb_db", GenbankAnnotationDb), ("gff_db", GffAnnotationDb)],
)
def test_read_old_format(db_name, cls, tmp_path, request):
    db = request.getfixturevalue(db_name)

    path = tmp_path / f"old_format_{db_name}.{cls._suffix}"

    old_tables = db.to_rich_dict()["tables"]

    table_name = db_name.split("_")[0]
    convert_spans_column(db, table_name)

    db.write(path)

    new_db = cls.from_file(path)
    with pytest.warns(UserWarning):
        new_tables = new_db.to_rich_dict()["tables"]
    assert new_tables == old_tables
    close_dbs(new_db, db)


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("annotations")


def test_load_anns_with_write(DATA_DIR, tmp_dir):
    inpath = DATA_DIR / "simple.gff"
    outpath = tmp_dir / "simple.c3gffdb"
    orig = load_annotations(path=inpath, write_path=outpath)
    expect = load_annotations(path=inpath)
    got = GffAnnotationDb.from_file(outpath)
    assert len(got) == len(expect)
    got_data = got.to_rich_dict()
    expect_data = expect.to_rich_dict()
    assert got_data["tables"] == expect_data["tables"]
    close_dbs(got, orig, expect)


def test_load_anns_from_json(DATA_DIR, tmp_dir):
    inpath = DATA_DIR / "simple.gff"
    orig = load_annotations(path=inpath)

    outpath = tmp_dir / "simple.json"
    with cogent3.open_(outpath, "w") as json_file:
        json_file.write(orig.to_json())

    got = load_annotations(path=outpath)
    assert len(got) == len(orig)
    assert got.to_rich_dict() == orig.to_rich_dict()
    close_dbs(got, orig)


def test_gff_end_renamed_to_stop(gff_db, tmp_path):
    bad_col_path = tmp_path / "bad_col.c3gffdb"

    correct_rich_dict = gff_db.to_rich_dict()
    del correct_rich_dict["init_args"]

    for table_name in gff_db.table_names:
        # Convert in memory db stop column to the old end name
        _rename_column_if_exists(gff_db.db, table_name, "stop", "end")

    gff_db.write(bad_col_path)

    # Test that end column is renamed to stop on construction
    loaded_gff_db = GffAnnotationDb.from_file(bad_col_path)

    new_rich_dict = loaded_gff_db.to_rich_dict()
    del new_rich_dict["init_args"]

    assert new_rich_dict == correct_rich_dict
    close_dbs(loaded_gff_db, gff_db)


def test_gb_end_renamed_to_stop(gb_db, tmp_path):
    bad_col_path = tmp_path / "bad_col.c3gbdb"

    correct_rich_dict = gb_db.to_rich_dict()
    del correct_rich_dict["init_args"]

    for table_name in gb_db.table_names:
        # Convert in memory db stop column to the old end name
        _rename_column_if_exists(gb_db.db, table_name, "stop", "end")

    gb_db.write(bad_col_path)

    # Test that end column is renamed to stop on construction
    loaded_gb_db = GenbankAnnotationDb.from_file(bad_col_path)

    new_rich_dict = loaded_gb_db.to_rich_dict()
    del new_rich_dict["init_args"]

    assert new_rich_dict == correct_rich_dict
    close_dbs(loaded_gb_db, gb_db)


def test_gbdb_get_children_fails_no_coords(gb_db):
    with pytest.raises(ValueError):
        _ = list(gb_db.get_feature_children(name="CNA00110"))


def test_gbdb_get_parent_fails_no_coords(gb_db):
    with pytest.raises(ValueError):
        _ = list(gb_db.get_feature_parent(name="CNA00110"))


@pytest.mark.parametrize("suffix", ["gff3", "gb"])
def test_load_annotations_invalid_path(suffix):
    with pytest.raises(IOError):
        load_annotations(path=f"invalidfile.{suffix}")


@pytest.mark.parametrize("integer", [int, numpy.int64])
def test_subset_gff3_db(gff_db, integer):
    subset = gff_db.subset(
        seqid="I",
        start=integer(40),
        stop=integer(70),
        allow_partial=True,
    )
    # manual inspection of the original GFF3 file indicates 7 records
    # BUT the CDS records get merged into a single row
    assert len(subset) == 6
    subset.close()


def test_subset_empty_db(gff_db):
    subset = gff_db.subset(seqid="X", start=40, stop=70, allow_partial=True)
    # no records
    assert not len(subset)
    subset.close()


def test_subset_gff3_db_with_user(gff_db):
    record = {
        "seqid": "I",
        "name": "gene-01",
        "biotype": "gene",
        "spans": [(23, 43)],
        "strand": "+",
    }
    gff_db.add_feature(**record)
    subset = gff_db.subset(seqid="I", start=40, stop=70, allow_partial=True)
    # manual inspection of the original GFF3 file indicates 7 records
    # BUT the CDS records get merged into a single row
    assert len(subset) == 7
    subset.close()


def test_subset_gb_db(gb_db):
    subset = gb_db.subset(biotype="gene")
    # manual inspection of the original annotated gb file indicates 2 genes
    assert len(subset) == 2
    subset.close()


def test_subset_gff3_db_source(gff_db, tmp_dir):
    outpath = tmp_dir / "subset.c3gffdb"
    subset = gff_db.subset(
        seqid="I",
        start=40,
        stop=70,
        allow_partial=True,
        source=outpath,
    )
    subset.db.close()

    # reload
    subset = type(gff_db).from_file(outpath)
    # manual inspection of the original GFF3 file indicates 7 records
    # BUT the CDS records get merged into a single row
    assert len(subset) == 6
    subset.close()


@pytest.mark.parametrize(
    "new_span",
    [[[22, 24]], [[9, 20]], [[0, 1]], [[9, 20], [29, 45], [59, 70]]],
)
def test_gffdb_update_record(gff_db, new_span):
    name = "CDS:B0019.1"
    init_value = [[9, 20], [29, 45], [59, 70]]
    combined = {tuple(c) for c in new_span + init_value}
    expect = [list(p) for p in sorted(combined)]
    gff_db.update_record_spans(name=name, spans=new_span)
    got = next(iter(gff_db.get_records_matching(name=name)))
    assert got["spans"].tolist() == expect


def test_gffdb_update_record_empty(gff_db):
    name = "CDS:B0019.1"
    init_value = [[9, 20], [29, 45], [59, 70]]
    gff_db.update_record_spans(name=name, spans=[])
    got = next(iter(gff_db.get_records_matching(name=name)))
    assert got["spans"].tolist() == init_value


def test_gffdb_update_absent_record(gff_db):
    name = "qwerty"  # does not exist
    init_value = [[9, 20], [29, 45], [59, 70]]
    gff_db.update_record_spans(name=name, spans=init_value)
    got = list(gff_db.get_records_matching(name=name))
    assert not got


@pytest.mark.parametrize("db", ["gff_db", "gb_db"])
def test_db_close(request, db):
    import sqlite3

    db = request.getfixturevalue(db)
    db.close()
    with pytest.raises((sqlite3.ProgrammingError, sqlite3.DatabaseError)):
        _ = list(db.get_features_matching(biotype="gene"))


@pytest.mark.parametrize("db", ["gff_db", "gb_db"])
@pytest.mark.parametrize(
    "col",
    ["biotype_id", "seqid_id", "name", "start", "stop", "parent_id"],
)
def test_db_make_index(request, db, col):
    ann_db = request.getfixturevalue(db)
    expect = {("index", col, tn) for tn in ann_db.table_names}
    ann_db.make_indexes()
    sql_template = (
        "SELECT * FROM sqlite_master WHERE type = 'index' AND"  # nosec B608
        f" tbl_name = '%s' and name = {col!r}"  # nosec B608
    )

    result = ann_db._execute_sql(sql_template % ann_db.table_names[0]).fetchone()
    got = tuple(result)[:3]
    assert got in expect
    ann_db.close()


@pytest.mark.parametrize("db", ["gff_db", "gb_db"])
def test_db_repr(request, db):
    ann_db = request.getfixturevalue(db)
    got = repr(ann_db)
    assert isinstance(got, str)
    ann_db.close()


@pytest.fixture
def home_gff(HOME_TMP_DIR, DATA_DIR) -> str:
    gff = DATA_DIR / "simple.gff"
    hgff = HOME_TMP_DIR / "simple.gff"
    hgff.expanduser().write_text(gff.read_text())
    return f"~/{HOME_TMP_DIR.name}/simple.gff"


@pytest.mark.parametrize("transform", [str, pathlib.Path])
def test_load_annotations_home_dir(home_gff, transform):
    got = load_annotations(path=transform(home_gff))
    assert len(got) == 6
    assert isinstance(got, GffAnnotationDb)
    got.close()


def test_load_annotations_invalid_backend():
    with pytest.raises(ValueError):
        load_annotations(path="somefile.txt", storage_backend="invalid_backend")


def test_load_annotations_invalid_file_format():
    with pytest.raises(ValueError):
        load_annotations(path="somefile.txt", format_name="invalid_format")


# Tests for optimized schema and upgrade_annotation_db function


def test_schema_version_new_db():
    """New databases should have the current schema version."""
    db = BasicAnnotationDb()
    assert db.schema_version == SCHEMA_VERSION
    db.close()

    db = GffAnnotationDb()
    assert db.schema_version == SCHEMA_VERSION
    db.close()

    db = GenbankAnnotationDb()
    assert db.schema_version == SCHEMA_VERSION
    db.close()


def test_lookup_tables_exist(DATA_DIR):
    """Lookup tables should exist in the schema."""
    db = load_annotations(path=DATA_DIR / "simple.gff")

    # Check that lookup tables exist (they may or may not have entries)
    conn = db.db
    cursor = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name IN ('seqids', 'biotypes')"
    )
    tables = [row[0] for row in cursor.fetchall()]
    assert "seqids" in tables
    assert "biotypes" in tables

    db.close()


def test_hierarchy_table_populated(DATA_DIR):
    """Hierarchy table should be populated for features with parent_id."""
    db = load_annotations(path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3")

    # This GFF has parent-child relationships
    conn = db.db
    cursor = conn.execute("SELECT COUNT(*) FROM feature_hierarchy")
    count = cursor.fetchone()[0]
    # The shortened GFF should have some hierarchy entries
    assert count > 0  # May be 0 if no parent_id in data
    db.close()


def test_upgrade_annotation_db(DATA_DIR, tmp_path):
    """Test upgrading an old-format database to the optimized schema."""
    import sqlite3

    # Create a database and write it to file
    db_path = tmp_path / "test_upgrade.c3gffdb"
    db = load_annotations(path=DATA_DIR / "simple.gff")
    db.write(db_path)
    db.close()

    # Simulate an old-format database by removing schema_version table
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    # Also drop optimized tables to simulate old format
    conn.execute("DROP TABLE IF EXISTS seqids")
    conn.execute("DROP TABLE IF EXISTS biotypes")
    conn.execute("DROP TABLE IF EXISTS sources")
    conn.execute("DROP TABLE IF EXISTS feature_spans")
    conn.execute("DROP TABLE IF EXISTS feature_hierarchy")
    conn.execute("DROP TABLE IF EXISTS feature_rtree")
    conn.commit()
    conn.close()

    # Upgrade the database
    upgrade_annotation_db(db_path, backup=True)

    # Verify backup was created
    backup_path = tmp_path / "test_upgrade.c3gffdb.bak"
    assert backup_path.exists()

    # Verify schema version is updated
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT version FROM schema_version ORDER BY version DESC LIMIT 1"
    )
    version = cursor.fetchone()[0]
    assert version == SCHEMA_VERSION

    # Verify optimized tables exist
    for table in ["seqids", "biotypes", "sources", "feature_hierarchy"]:
        cursor = conn.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table}'"
        )
        assert cursor.fetchone() is not None, f"Table {table} should exist"

    conn.close()


def test_upgrade_annotation_db_idempotent(DATA_DIR, tmp_path):
    """Test that upgrading twice is safe (idempotent)."""
    import sqlite3

    # Create a database and write it to file
    db_path = tmp_path / "test_idempotent.c3gffdb"
    db = load_annotations(path=DATA_DIR / "simple.gff")
    db.write(db_path)
    db.close()

    # Simulate old format
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    conn.commit()
    conn.close()

    # Upgrade once
    upgrade_annotation_db(db_path, backup=True)

    # Upgrade again - should not fail and should not create another backup
    upgrade_annotation_db(db_path, backup=False)

    # Verify still at current version
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT version FROM schema_version ORDER BY version DESC LIMIT 1"
    )
    version = cursor.fetchone()[0]
    assert version == SCHEMA_VERSION
    conn.close()


def test_upgrade_annotation_db_file_not_exists(tmp_path):
    """Test that upgrade fails gracefully when file doesn't exist."""
    with pytest.raises(OSError):
        upgrade_annotation_db(tmp_path / "nonexistent.gffdb")


def test_upgrade_annotation_db_backup_exists(DATA_DIR, tmp_path):
    """Test that upgrade fails when backup already exists."""
    # Create database and backup file
    db_path = tmp_path / "test_backup.c3gffdb"
    backup_path = tmp_path / "test_backup.c3gffdb.bak"

    db = load_annotations(path=DATA_DIR / "simple.gff")
    db.write(db_path)
    db.close()

    # Create a fake backup file
    backup_path.touch()

    # Simulate old format
    import sqlite3

    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    conn.commit()
    conn.close()

    with pytest.raises(FileExistsError):
        upgrade_annotation_db(db_path, backup=True)


def test_upgrade_annotation_db_genbank(DATA_DIR, tmp_path):
    """Test upgrading a GenBank annotation database."""
    import sqlite3

    # Create a GenBank database
    db_path = tmp_path / "test_gb.c3gbdb"
    db = load_annotations(path=DATA_DIR / "annotated_seq.gb")
    db.write(db_path)
    db.close()

    # Simulate old format
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    conn.commit()
    conn.close()

    # Upgrade
    upgrade_annotation_db(db_path, backup=True)

    # Verify schema version
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT version FROM schema_version ORDER BY version DESC LIMIT 1"
    )
    version = cursor.fetchone()[0]
    assert version == SCHEMA_VERSION
    conn.close()


def test_rtree_range_query_matches_original(DATA_DIR):
    """R*Tree range queries should return same results as original queries."""
    db = load_annotations(path=DATA_DIR / "simple.gff")

    # Get features with R*Tree
    features_rtree = list(
        db.get_features_matching(start=0, stop=10, allow_partial=True)
    )

    # Disable R*Tree and get features again
    db._use_rtree = False
    features_no_rtree = list(
        db.get_features_matching(start=0, stop=10, allow_partial=True)
    )

    # Both should return the same results
    assert len(features_rtree) == len(features_no_rtree)

    # Compare feature names
    rtree_names = {f["name"] for f in features_rtree}
    no_rtree_names = {f["name"] for f in features_no_rtree}
    assert rtree_names == no_rtree_names

    db.close()


def test_hierarchy_children_query_matches_original(DATA_DIR):
    """Hierarchy-based children queries should return same results as LIKE-based."""
    db = load_annotations(path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3")

    # Get a gene name that has children
    gene_name = "Gene:WBGene00000001"

    # Get children with hierarchy table
    children_hierarchy = list(db.get_feature_children(gene_name))

    # Disable hierarchy and get children again
    original_version = db._schema_version
    db._schema_version = 1  # Force fallback to LIKE-based query
    children_like = list(db.get_feature_children(gene_name))
    db._schema_version = original_version

    # Both should return the same results (or both empty if no children)
    assert len(children_hierarchy) == len(children_like)

    db.close()


def test_data_integrity_after_upgrade(DATA_DIR, tmp_path):
    """Data should be unchanged after upgrade."""
    import sqlite3

    # Create a database and write it to file
    db_path = tmp_path / "test_integrity.c3gffdb"
    db = load_annotations(path=DATA_DIR / "simple.gff")

    # Get feature count before
    features_before = list(db.get_features_matching())
    count_before = len(features_before)

    db.write(db_path)
    db.close()

    # Simulate old format
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    conn.commit()
    conn.close()

    # Upgrade
    upgrade_annotation_db(db_path, backup=True)

    # Load upgraded database and count features
    db = GffAnnotationDb.from_file(path=db_path)
    features_after = list(db.get_features_matching())
    count_after = len(features_after)

    assert count_before == count_after

    # Check that specific features are preserved
    names_before = {f["name"] for f in features_before}
    names_after = {f["name"] for f in features_after}
    assert names_before == names_after

    db.close()


@pytest.mark.parametrize(
    ("db_name", "cls", "suffix"),
    [
        ("gff_db", GffAnnotationDb, "c3gffdb"),
        ("gb_db", GenbankAnnotationDb, "c3gbdb"),
        ("anno_db", BasicAnnotationDb, "c3andb"),
    ],
)
def test_constructor_rejects_existing_db_file(db_name, cls, suffix, tmp_path, request):
    """Constructing a db class with an existing file path should raise ValueError.

    Users should use <ClassName>.from_file() to load existing databases,
    not the constructor with source=path.
    """
    db = request.getfixturevalue(db_name)
    outpath = tmp_path / f"existing_db.{suffix}"
    db.write(outpath)

    # This should raise ValueError - users should use from_file() instead
    with pytest.raises(ValueError):
        cls(source=outpath)


# Tests for load_annotations plugin delegation with custom suffixes


@pytest.mark.parametrize(
    ("db_name", "cls", "suffix"),
    [
        ("gff_db", GffAnnotationDb, "c3gffdb"),
        ("gb_db", GenbankAnnotationDb, "c3gbdb"),
        ("anno_db", BasicAnnotationDb, "c3andb"),
    ],
)
def test_load_annotations_saved_db_file(db_name, cls, suffix, tmp_path, request):
    """load_annotations should correctly load saved database files by suffix."""
    db = request.getfixturevalue(db_name)
    outpath = tmp_path / f"saved_db.{suffix}"
    db.write(outpath)

    # load_annotations should delegate to the correct from_file method
    loaded = load_annotations(path=outpath)
    assert isinstance(loaded, cls)
    assert loaded.num_matches() == db.num_matches()
    close_dbs(loaded)


@pytest.mark.parametrize(
    ("db_name", "cls", "suffix"),
    [
        ("gff_db", GffAnnotationDb, "c3gffdb"),
        ("gb_db", GenbankAnnotationDb, "c3gbdb"),
        ("anno_db", BasicAnnotationDb, "c3andb"),
    ],
)
def test_load_annotations_with_storage_backend(db_name, cls, suffix, tmp_path, request):
    """load_annotations with storage_backend='c3anndb' should load saved db files."""
    db = request.getfixturevalue(db_name)
    outpath = tmp_path / f"saved_db.{suffix}"
    db.write(outpath)

    # Explicitly specify the storage backend
    loaded = load_annotations(path=outpath, storage_backend="c3anndb")
    assert isinstance(loaded, cls)
    assert loaded.num_matches() == db.num_matches()
    close_dbs(loaded)


def test_load_annotations_gff_with_c3anndb_backend(DATA_DIR):
    """load_annotations with storage_backend='c3anndb' should load GFF files."""
    path = DATA_DIR / "simple.gff"
    db = load_annotations(path=path, storage_backend="c3anndb")
    assert isinstance(db, GffAnnotationDb)
    assert db.num_matches() > 0
    db.close()


def test_load_annotations_gb_with_c3anndb_backend(DATA_DIR):
    """load_annotations with storage_backend='c3anndb' should load GenBank files."""
    path = DATA_DIR / "annotated_seq.gb"
    db = load_annotations(path=path, storage_backend="c3anndb")
    assert isinstance(db, GenbankAnnotationDb)
    assert db.num_matches() > 0
    db.close()


def test_load_annotations_roundtrip_gff(DATA_DIR, tmp_path):
    """Round-trip test: load GFF, save to c3gffdb, reload via load_annotations."""
    # Load from GFF
    orig = load_annotations(path=DATA_DIR / "simple.gff")
    outpath = tmp_path / "roundtrip.c3gffdb"
    orig.write(outpath)

    # Reload via load_annotations (should auto-detect suffix)
    reloaded = load_annotations(path=outpath)
    assert isinstance(reloaded, GffAnnotationDb)
    assert reloaded.num_matches() == orig.num_matches()

    # Verify data integrity
    orig_features = list(orig.get_features_matching())
    reloaded_features = list(reloaded.get_features_matching())
    assert len(orig_features) == len(reloaded_features)

    close_dbs(orig, reloaded)


def test_load_annotations_roundtrip_gb(DATA_DIR, tmp_path):
    """Round-trip test: load GenBank, save to c3gbdb, reload via load_annotations."""
    # Load from GenBank
    orig = load_annotations(path=DATA_DIR / "annotated_seq.gb")
    outpath = tmp_path / "roundtrip.c3gbdb"
    orig.write(outpath)

    # Reload via load_annotations (should auto-detect suffix)
    reloaded = load_annotations(path=outpath)
    assert isinstance(reloaded, GenbankAnnotationDb)
    assert reloaded.num_matches() == orig.num_matches()

    close_dbs(orig, reloaded)


# Tests for LookupTableCache._get_name_by_id coverage


def test_lookup_cache_get_seqid_name_cache_hit():
    """Test get_seqid_name returns cached value when ID is in reverse cache."""
    conn = _make_db_connection(":memory:")
    conn.execute("CREATE TABLE seqids (id INTEGER PRIMARY KEY, name TEXT UNIQUE)")
    conn.execute("INSERT INTO seqids (name) VALUES ('chr1')")
    conn.commit()

    cache = LookupTableCache(conn)
    # First call populates the cache
    seqid_id = cache.get_seqid_id("chr1")
    # Second call should hit the reverse cache
    name = cache.get_seqid_name(seqid_id)
    assert name == "chr1"

    conn.close()


def test_lookup_cache_get_seqid_name_db_lookup():
    """Test get_seqid_name performs DB lookup when ID not in cache."""
    conn = _make_db_connection(":memory:")
    conn.execute("CREATE TABLE seqids (id INTEGER PRIMARY KEY, name TEXT UNIQUE)")
    conn.execute("INSERT INTO seqids (name) VALUES ('chr1')")
    conn.commit()

    # Create cache and manually insert ID without populating reverse cache
    cache = LookupTableCache(conn)
    # Directly query DB to get ID without populating cache
    cursor = conn.execute("SELECT id FROM seqids WHERE name = 'chr1'")
    seqid_id = cursor.fetchone()[0]

    # get_seqid_name should look up from DB
    name = cache.get_seqid_name(seqid_id)
    assert name == "chr1"
    # Verify it was cached
    assert seqid_id in cache._seqid_reverse
    assert "chr1" in cache._seqid_cache

    conn.close()


def test_lookup_cache_get_seqid_name_not_found():
    """Test get_seqid_name returns None when ID doesn't exist."""
    conn = _make_db_connection(":memory:")
    conn.execute("CREATE TABLE seqids (id INTEGER PRIMARY KEY, name TEXT UNIQUE)")
    conn.commit()

    cache = LookupTableCache(conn)
    # Query for non-existent ID
    name = cache.get_seqid_name(999)
    assert name is None

    conn.close()


def test_lookup_cache_get_biotype_name():
    """Test get_biotype_name works correctly."""
    conn = _make_db_connection(":memory:")
    conn.execute("CREATE TABLE biotypes (id INTEGER PRIMARY KEY, name TEXT UNIQUE)")
    conn.execute("INSERT INTO biotypes (name) VALUES ('gene')")
    conn.commit()

    cache = LookupTableCache(conn)
    # First call populates the cache
    biotype_id = cache.get_biotype_id("gene")
    # Get name from ID
    name = cache.get_biotype_name(biotype_id)
    assert name == "gene"

    # Test not found case
    assert cache.get_biotype_name(999) is None

    conn.close()


def test_lookup_cache_get_source_name():
    """Test get_source_name works correctly."""
    conn = _make_db_connection(":memory:")
    conn.execute("CREATE TABLE sources (id INTEGER PRIMARY KEY, name TEXT UNIQUE)")
    conn.execute("INSERT INTO sources (name) VALUES ('GFF3')")
    conn.commit()

    cache = LookupTableCache(conn)
    # First call populates the cache
    source_id = cache.get_source_id("GFF3")
    # Get name from ID
    name = cache.get_source_name(source_id)
    assert name == "GFF3"

    # Test not found case
    assert cache.get_source_name(999) is None

    conn.close()


# Tests for SqliteAnnotationDbMixin._get_feature_rowid coverage


def test_get_feature_rowid_by_name(gff_db):
    """Test _get_feature_rowid finds feature by name."""
    # Get a known feature name from the database
    features = list(gff_db.get_features_matching(biotype="gene"))
    if features:
        feature_name = features[0]["name"]
        rowid = gff_db._get_feature_rowid("gff", feature_name)
        assert rowid is not None
        assert isinstance(rowid, int)


def test_get_feature_rowid_by_name_and_biotype(gff_db):
    """Test _get_feature_rowid finds feature by name and biotype."""
    # Get a known feature
    features = list(gff_db.get_features_matching(biotype="gene"))
    if features:
        feature = features[0]
        rowid = gff_db._get_feature_rowid(
            "gff", feature["name"], biotype=feature["biotype"]
        )
        assert rowid is not None
        assert isinstance(rowid, int)


def test_get_feature_rowid_not_found(gff_db):
    """Test _get_feature_rowid returns None when feature doesn't exist."""
    rowid = gff_db._get_feature_rowid("gff", "nonexistent_feature_xyz")
    assert rowid is None


def test_get_feature_rowid_wrong_biotype(gff_db):
    """Test _get_feature_rowid returns None when biotype doesn't match."""
    # Get a gene feature
    features = list(gff_db.get_features_matching(biotype="gene"))
    if features:
        feature_name = features[0]["name"]
        # Search with wrong biotype
        rowid = gff_db._get_feature_rowid(
            "gff", feature_name, biotype="nonexistent_biotype"
        )
        assert rowid is None


# Tests for upgrade_annotation_db hierarchy building coverage


def test_upgrade_annotation_db_builds_hierarchy(tmp_path):
    """Test that upgrade_annotation_db correctly builds feature hierarchy."""
    import sqlite3

    # Create a database with parent-child relationships
    db_path = tmp_path / "hierarchy_test.c3gffdb"

    # Create database manually with old schema (no hierarchy table)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Create gff table (simplified version of actual schema)
    conn.execute("""
        CREATE TABLE gff (
            seqid TEXT,
            source TEXT,
            biotype TEXT,
            start INTEGER,
            stop INTEGER,
            score TEXT,
            strand TEXT,
            phase TEXT,
            name TEXT,
            parent_id TEXT
        )
    """)

    # Insert parent gene
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'gene', 100, 500, '+', 'gene1', NULL)
    """)

    # Insert child mRNA with parent_id pointing to gene1
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'mRNA', 100, 500, '+', 'mRNA1', 'gene1')
    """)

    # Insert exon with parent_id pointing to mRNA1
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'exon', 100, 200, '+', 'exon1', 'mRNA1')
    """)

    conn.commit()
    conn.close()

    # Upgrade the database
    upgrade_annotation_db(db_path, backup=False)

    # Verify hierarchy was built
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Check feature_hierarchy table exists and has entries
    cursor = conn.execute("SELECT COUNT(*) FROM feature_hierarchy")
    count = cursor.fetchone()[0]
    assert count >= 2  # At least gene->mRNA and mRNA->exon relationships

    # Check specific relationships
    cursor = conn.execute("""
        SELECT h.child_id, h.parent_id
        FROM feature_hierarchy h
        JOIN gff g_child ON g_child.rowid = h.child_id
        JOIN gff g_parent ON g_parent.rowid = h.parent_id
        WHERE g_child.name = 'mRNA1' AND g_parent.name = 'gene1'
    """)
    row = cursor.fetchone()
    assert row is not None

    conn.close()


def test_upgrade_annotation_db_multiple_parents(tmp_path):
    """Test upgrade handles features with comma-separated parent_ids."""
    import sqlite3

    db_path = tmp_path / "multi_parent_test.c3gffdb"

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    conn.execute("""
        CREATE TABLE gff (
            seqid TEXT,
            source TEXT,
            biotype TEXT,
            start INTEGER,
            stop INTEGER,
            score TEXT,
            strand TEXT,
            phase TEXT,
            name TEXT,
            parent_id TEXT
        )
    """)

    # Insert two parent genes
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'gene', 100, 500, '+', 'gene1', NULL)
    """)
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'gene', 600, 1000, '+', 'gene2', NULL)
    """)

    # Insert exon with multiple parents (comma-separated)
    conn.execute("""
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, parent_id)
        VALUES ('chr1', 'test', 'exon', 100, 200, '+', 'shared_exon', 'gene1,gene2')
    """)

    conn.commit()
    conn.close()

    # Upgrade the database
    upgrade_annotation_db(db_path, backup=False)

    # Verify hierarchy has two entries for the shared exon
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    cursor = conn.execute("""
        SELECT COUNT(*)
        FROM feature_hierarchy h
        JOIN gff g ON g.rowid = h.child_id
        WHERE g.name = 'shared_exon'
    """)
    count = cursor.fetchone()[0]
    assert count == 2  # Should have two parent relationships

    conn.close()


def test_upgrade_annotation_db_no_backup(DATA_DIR, tmp_path):
    """Test upgrade_annotation_db with backup=False doesn't create backup file."""
    import sqlite3

    db_path = tmp_path / "no_backup_test.c3gffdb"
    db = load_annotations(path=DATA_DIR / "simple.gff")
    db.write(db_path)
    db.close()

    # Simulate old format
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS schema_version")
    conn.commit()
    conn.close()

    # Upgrade without backup
    upgrade_annotation_db(db_path, backup=False)

    # Verify no backup was created
    backup_path = tmp_path / "no_backup_test.c3gffdb.bak"
    assert not backup_path.exists()

    # Verify schema was still upgraded
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT version FROM schema_version ORDER BY version DESC LIMIT 1"
    )
    version = cursor.fetchone()[0]
    assert version == SCHEMA_VERSION
    conn.close()


def test_upgrade_annotation_db_migrates_columns(tmp_path):
    """Test that upgrade_annotation_db migrates old TEXT columns to new *_id columns.

    Old databases have seqid/biotype/source as TEXT columns.
    New databases have seqid_id/biotype_id/source_id as INTEGER columns
    referencing lookup tables. The upgrade function must migrate the schema.
    """
    import sqlite3

    db_path = tmp_path / "old_schema_test.c3gffdb"

    # Create database with OLD schema (TEXT columns, not *_id INTEGER columns)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Old schema used seqid, biotype, source as TEXT columns (not *_id INTEGER)
    conn.execute("""
        CREATE TABLE gff (
            seqid TEXT,
            source TEXT,
            biotype TEXT,
            start INTEGER,
            stop INTEGER,
            score TEXT,
            strand INTEGER,
            phase TEXT,
            attributes TEXT,
            comments TEXT,
            spans array,
            name TEXT,
            parent_id TEXT
        )
    """)

    # Also need user table with old schema
    conn.execute("""
        CREATE TABLE user (
            seqid TEXT,
            biotype TEXT,
            name TEXT,
            parent_id TEXT,
            start INTEGER,
            stop INTEGER,
            strand INTEGER,
            spans array,
            attributes TEXT,
            on_alignment INT
        )
    """)

    # Insert test data with the old schema
    # Need to serialize spans as numpy array
    from cogent3.core.annotation_db import array_to_sqlite

    spans_data = array_to_sqlite(numpy.array([[100, 500]], dtype=int))
    conn.execute(
        """
        INSERT INTO gff (seqid, source, biotype, start, stop, strand, name, spans)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        ("chr1", "test_source", "gene", 100, 500, 1, "gene1", spans_data),
    )
    conn.commit()
    conn.close()

    # Upgrade the database
    upgrade_annotation_db(db_path, backup=False)

    # Now try to load and use the database - this should work after proper upgrade
    loaded_db = GffAnnotationDb.from_file(db_path)

    # Query features - this will fail if columns weren't migrated properly
    # because the code will try to query seqid_id but the table still has seqid
    features = list(loaded_db.get_features_matching(seqid="chr1"))
    assert len(features) == 1
    assert features[0]["name"] == "gene1"
    assert features[0]["biotype"] == "gene"
    assert features[0]["seqid"] == "chr1"

    # Also verify we can query by biotype
    features = list(loaded_db.get_features_matching(biotype="gene"))
    assert len(features) == 1

    loaded_db.close()


@pytest.fixture
def chunked_vs_single_gff_dbs(DATA_DIR):
    """Yields (single_block_db, chunked_db) for the c_elegans GFF.

    ``CDS:B0019.1`` spans three lines in this file, so a ``num_lines=2`` chunk
    size is guaranteed to straddle the boundary. The reference 'expected' DB
    loads the same file in a single block.
    """
    from cogent3.core.annotation_db import _db_from_gff

    path = DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"
    expected = load_annotations(path=path)
    chunked = _db_from_gff(
        path=path, seqids=None, db=None, write_path=":memory:", num_lines=2
    )
    yield expected, chunked
    expected.close()
    chunked.close()


def test_chunked_load_no_duplicate_rows(chunked_vs_single_gff_dbs):
    # make sure we avoid issues where a feature whose lines straddle
    # a block boundary is re-inserted as a duplicate row.

    expected, chunked = chunked_vs_single_gff_dbs

    assert len(chunked) == len(expected)

    cds_rows = chunked.db.execute(
        "SELECT COUNT(*) FROM gff WHERE name = ?", ("CDS:B0019.1",)
    ).fetchone()[0]
    assert cds_rows == 1, f"CDS:B0019.1 appears {cds_rows} times in gff table"


def test_chunked_load_hierarchy_attached_to_canonical_row(
    chunked_vs_single_gff_dbs,
):
    # make sure name lookup points to the canonical row.

    _, chunked = chunked_vs_single_gff_dbs

    canonical_cds_rowid = chunked.db.execute(
        "SELECT MIN(rowid) FROM gff WHERE name = ?", ("CDS:B0019.1",)
    ).fetchone()[0]

    edge_targets = {
        row[0]
        for row in chunked.db.execute(
            "SELECT child_id FROM feature_hierarchy WHERE child_id IN "
            "(SELECT rowid FROM gff WHERE name = ?)",
            ("CDS:B0019.1",),
        ).fetchall()
    }
    assert canonical_cds_rowid in edge_targets, (
        f"hierarchy edge missing for canonical CDS rowid {canonical_cds_rowid}; "
        f"edges instead attached to: {edge_targets}"
    )


def test_chunked_load_preserves_forward_parent_references(
    chunked_vs_single_gff_dbs,
):
    # make sure we handle the case where a parent occurs after a child
    expected, chunked = chunked_vs_single_gff_dbs

    expected_children = sorted(
        c["name"] for c in expected.get_feature_children(name="Gene:WBGene00000138")
    )
    chunked_children = sorted(
        c["name"] for c in chunked.get_feature_children(name="Gene:WBGene00000138")
    )
    assert chunked_children == expected_children


def test_load_annotations_format_name(tmp_path, DATA_DIR):
    path = DATA_DIR / "annotated_seq.gb"
    outpath = tmp_path / "annotated_seq.gb"
    outpath.write_text(path.read_text())
    db = load_annotations(path=outpath, format_name="genbank")
    assert isinstance(db, GenbankAnnotationDb)


def test_load_annotations_unknown_format_name(tmp_path, DATA_DIR):
    path = DATA_DIR / "annotated_seq.gb"
    outpath = tmp_path / "annotated_seq.gb"
    outpath.write_text(path.read_text())
    with pytest.raises(ValueError):
        load_annotations(path=outpath, format_name="bogus")

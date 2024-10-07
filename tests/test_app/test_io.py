import bz2
import gzip
import json
import os
import pathlib
import pickle
import shutil
import tempfile

import numpy
import pytest
from numpy.testing import assert_allclose

from cogent3 import get_app, get_moltype, open_data_store
from cogent3.app import io as io_app
from cogent3.app.composable import NotCompleted, source_proxy
from cogent3.app.data_store import (
    DataMember,
    DataStoreDirectory,
    Mode,
    ReadOnlyDataStoreZipped,
)
from cogent3.app.io import DEFAULT_DESERIALISER, DEFAULT_SERIALISER
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.maths.util import safe_log
from cogent3.parse.sequence import PARSERS
from cogent3.util.deserialise import deserialise_object
from cogent3.util.table import Table

DNA = get_moltype("dna")


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("io")


@pytest.fixture(scope="function")
def w_dir_dstore(tmp_dir):
    return DataStoreDirectory(tmp_dir, mode="w")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    # this set's the working directory for all tests in this module
    # as a tmp dir
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="function")
def fasta_dir(DATA_DIR, tmp_dir):
    tmp_dir = pathlib.Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    fasta_dir = tmp_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = fasta_dir / fn.name
        dest.write_text(fn.read_text())
    return fasta_dir


@pytest.fixture
def zipped_full(fasta_dir):
    # converts the fasta_dir into a zipped archive
    source = fasta_dir
    path = shutil.make_archive(
        base_name=source, format="zip", base_dir=source, root_dir=source.parent
    )
    return ReadOnlyDataStoreZipped(pathlib.Path(path), suffix="fasta")


def _get_generic_result(source):
    """creates a generic result with a DNA moltype as the single value"""
    from cogent3.app.result import generic_result

    gr = generic_result(source=source)
    gr["dna"] = DNA
    return gr


def test_write_seqs(fasta_dir, tmp_dir):
    """correctly writes sequences out"""
    datastore = DataStoreDirectory(fasta_dir, suffix="fasta")
    datamember = datastore[0]
    data = datamember.read().splitlines()
    data = dict(iter(PARSERS["fasta".lower()](data)))
    seqs = ArrayAlignment(data=data, moltype=None)
    seqs.info.source = datastore.source
    out_data_store = DataStoreDirectory(
        tmp_dir / "test_write_seqs", mode=Mode.w, suffix="fasta"
    )
    writer = io_app.write_seqs(out_data_store, format="fasta")
    wrote = writer(seqs[0], identifier=datamember.unique_id)
    assert isinstance(wrote, DataMember)


def test_source_proxy_simple(fasta_dir):
    """correctly writes sequences out"""
    from cogent3.app.composable import define_app
    from cogent3.app.typing import IdentifierType

    @define_app
    def get_bytes(path: IdentifierType) -> bytes:
        path = pathlib.Path(path)
        return path.read_bytes()

    datastore = DataStoreDirectory(fasta_dir, suffix="fasta")
    datamember = datastore[0]
    reader = get_bytes()
    path = datamember.data_store.source / datamember.unique_id
    data = reader(path)
    # direct call gives you back the annotated type
    assert isinstance(data, (bytes, bytearray))
    # directly calling the intermediate wrap method should work
    got = reader._source_wrapped(source_proxy(path))
    assert isinstance(got, source_proxy)
    # calling with list of data that doesn't have a source should
    # also return source_proxy
    got = list(reader.as_completed([path], show_progress=False))
    assert isinstance(got[0], source_proxy)


@pytest.mark.parametrize("suffix", ("nex", "paml", "fasta"))
def test_load_aligned(DATA_DIR, suffix):
    """should handle nexus too"""
    dstore = DataStoreDirectory(DATA_DIR, suffix=suffix, limit=2)
    loader = io_app.load_aligned(format=suffix)
    results = [loader(m) for m in dstore]
    for result in results:
        assert isinstance(result, ArrayAlignment)


def test_load_unaligned(DATA_DIR):
    """load_unaligned returns degapped sequence collections"""
    fasta_paths = DataStoreDirectory(DATA_DIR, suffix=".fasta", limit=2)
    fasta_loader = io_app.load_unaligned(format="fasta")
    for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
        assert isinstance(seqs, SequenceCollection)
        assert "-" not in "".join(seqs.to_dict().values())
        assert seqs.info.source == fasta_paths[i].unique_id


@pytest.mark.parametrize(
    "loader",
    (io_app.load_aligned, io_app.load_unaligned, io_app.load_tabular, io_app.load_db),
)
def test_load_nonpath(loader):
    # returns NotCompleted when it's given an alignment/sequence
    # collection
    got = loader()({})
    assert isinstance(got, NotCompleted)


def test_load_json(tmp_dir):
    """correctly loads an object from json"""
    import json

    from cogent3.app.data_store import make_record_for_json

    data = make_record_for_json("delme", DNA, True)
    data = json.dumps(data)
    outpath = tmp_dir / "delme.json"
    outpath.write_text(data)
    # straight directory
    reader = io_app.load_json()
    got = reader(outpath)
    assert isinstance(got, DNA.__class__)
    assert got == DNA


def test_load_tabular(tmp_dir):
    """correctly loads tabular data"""
    rows = [[1, 2], [3, 4], [5, 6.5]]
    table = Table(["A", "B"], data=rows)
    load_table = io_app.load_tabular(sep="\t", with_header=True)
    outpath = tmp_dir / "delme.tsv"
    table.write(outpath)
    new = load_table(outpath)
    assert new.title == ""
    assert type(new[0, "B"]) == type(table[0, "B"])
    assert type(new[0, "A"]) == type(table[0, "A"])
    outpath = tmp_dir / "delme2.tsv"
    with open(outpath, "w") as out:
        out.write("\t".join(table.header[:1]) + "\n")
        for row in table.to_list():
            row = "\t".join(map(str, row))
            out.write(row + "\n")
    result = load_table(outpath)
    assert isinstance(result, NotCompleted)


def test_load_tabular_motif_counts_array(w_dir_dstore):
    """correctly loads tabular data for MotifCountsArray"""
    w_dir_dstore.suffix = "tsv"
    data = [[2, 4], [3, 5], [4, 8]]
    mca = MotifCountsArray(data, "AB")
    loader = io_app.load_tabular(sep="\t", as_type="motif_counts")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    outpath = "delme.tsv"
    m = writer.main(data=mca, identifier="delme")
    new = loader(m)
    assert mca.to_dict() == new.to_dict()


def test_load_tabular_motif_freqs_array(w_dir_dstore):
    """correctly loads tabular data for MotifFreqsArray"""
    w_dir_dstore.suffix = "tsv"
    data = [[0.3333, 0.6667], [0.3750, 0.625], [0.3333, 0.6667]]
    mfa = MotifFreqsArray(data, "AB")
    loader = io_app.load_tabular(sep="\t", as_type="motif_freqs")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    outpath = "delme"
    m = writer.main(mfa, identifier="delme")
    new = loader(m)
    assert mfa.to_dict() == new.to_dict()


def test_load_tabular_pssm(w_dir_dstore):
    """correctly loads tabular data for PSSM"""
    w_dir_dstore.suffix = "tsv"
    # data from test_profile
    data = [
        [0.1, 0.3, 0.5, 0.1],
        [0.25, 0.25, 0.25, 0.25],
        [0.05, 0.8, 0.05, 0.1],
        [0.7, 0.1, 0.1, 0.1],
        [0.6, 0.15, 0.05, 0.2],
    ]
    pssm = PSSM(data, "ACTG")
    loader = io_app.load_tabular(sep="\t", as_type="pssm")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(pssm, identifier="delme")
    new = loader(m)
    assert_allclose(pssm.array, new.array, atol=0.0001)


def test_load_tabular_distance_matrix(w_dir_dstore):
    """correctly loads tabular data for DistanceMatrix"""
    w_dir_dstore.suffix = "tsv"
    data = {(0, 0): 0, (0, 1): 4, (1, 0): 4, (1, 1): 0}
    matrix = DistanceMatrix(data)
    loader = io_app.load_tabular(sep="\t", as_type="distances")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(matrix, identifier="delme")
    new = loader(m)
    assert matrix.to_dict() == new.to_dict()


def test_load_tabular_table(w_dir_dstore):
    """correctly loads tabular data"""
    w_dir_dstore.suffix = "tsv"
    rows = [[1, 2], [3, 4], [5, 6.5]]
    table = Table(["A", "B"], data=rows)
    loader = io_app.load_tabular(sep="\t", as_type="table")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(table, identifier="delme")
    new = loader(m)
    assert table.to_dict() == new.to_dict()


def test_write_tabular_motif_counts_array(w_dir_dstore):
    """correctly writes tabular data for MotifCountsArray"""
    w_dir_dstore.suffix = "tsv"
    data = [[2, 4], [3, 5], [4, 8]]
    mca = MotifCountsArray(data, "AB")
    loader = io_app.load_tabular(sep="\t")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(mca, identifier="delme")
    assert isinstance(m, DataMember)
    new = loader(m)
    # when written to file in tabular form
    # the loaded table will have dim-1 dim-2 as column labels
    # and the key-values pairs listed below; in dict form...
    expected = {
        0: {"dim-1": 0, "dim-2": "A", "value": 2},
        1: {"dim-1": 0, "dim-2": "B", "value": 4},
        2: {"dim-1": 1, "dim-2": "A", "value": 3},
        3: {"dim-1": 1, "dim-2": "B", "value": 5},
        4: {"dim-1": 2, "dim-2": "A", "value": 4},
        5: {"dim-1": 2, "dim-2": "B", "value": 8},
    }
    assert expected == new.to_dict()


def test_write_tabular_motif_freqs_array(w_dir_dstore):
    """correctly writes tabular data for MotifFreqsArray"""
    w_dir_dstore.suffix = "tsv"
    data = [[0.3333, 0.6667], [0.3750, 0.625], [0.3333, 0.6667]]
    mfa = MotifFreqsArray(data, "AB")
    loader = io_app.load_tabular(sep="\t")

    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(mfa, identifier="delme")
    new = loader(m)
    # when written to file in tabular form
    # the loaded table will have dim-1 dim-2 as column labels
    # and the key-values pairs listed below; in dict form...
    expected = {
        0: {"dim-1": 0, "dim-2": "A", "value": 0.3333},
        1: {"dim-1": 0, "dim-2": "B", "value": 0.6667},
        2: {"dim-1": 1, "dim-2": "A", "value": 0.3750},
        3: {"dim-1": 1, "dim-2": "B", "value": 0.6250},
        4: {"dim-1": 2, "dim-2": "A", "value": 0.3333},
        5: {"dim-1": 2, "dim-2": "B", "value": 0.6667},
    }
    assert expected == new.to_dict()


def test_write_tabular_pssm(w_dir_dstore):
    """correctly writes tabular data for PSSM"""
    w_dir_dstore.suffix = "tsv"
    # data from test_profile
    data = numpy.array(
        [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
    )
    pssm = PSSM(data, "ACTG")
    loader = io_app.load_tabular(sep="\t")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(pssm, identifier="delme")
    new = loader(m)
    expected = safe_log(data) - safe_log(numpy.array([0.25, 0.25, 0.25, 0.25]))
    for i in range(len(expected)):
        j = i // 4
        assert numpy.isclose(new.array[i][2], expected[j][i - j], atol=0.0001)


def test_write_tabular_distance_matrix(w_dir_dstore):
    """correctly writes tabular data for DistanceMatrix"""
    w_dir_dstore.suffix = "tsv"
    data = {(0, 0): 0, (0, 1): 4, (1, 0): 4, (1, 1): 0}
    matrix = DistanceMatrix(data)
    loader = io_app.load_tabular(sep="\t")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(matrix, identifier="delme")
    new = loader(m)
    # when written to file in tabular form
    # the loaded table will have dim-1 dim-2 as column labels
    # and the key-values pairs listed below; in dict form...
    expected = {
        0: {"dim-1": 0, "dim-2": 1, "value": 4},
        1: {"dim-1": 1, "dim-2": 0, "value": 4},
    }
    assert expected == new.to_dict()


def test_write_tabular_table(w_dir_dstore):
    """correctly writes tabular data"""
    w_dir_dstore.suffix = "tsv"
    rows = [[1, 2], [3, 4], [5, 6.5]]
    table = Table(["A", "B"], data=rows)
    loader = io_app.load_tabular(sep="\t")
    writer = io_app.write_tabular(data_store=w_dir_dstore, format="tsv")
    m = writer.main(table, identifier="delme")
    new = loader(m)
    assert table.to_dict() == new.to_dict()


def test_write_json_with_info(w_dir_dstore):
    """correctly writes an object with source attribute from json"""
    w_dir_dstore.suffix = "json"
    # create a mock object that pretends like it's been derived from
    # something

    obj = _get_generic_result("delme.json")
    writer = io_app.write_json(data_store=w_dir_dstore)
    m = writer(obj)
    assert isinstance(m, DataMember)
    reader = io_app.load_json()
    got = reader(m)
    got.deserialised_values()
    assert got["dna"] == DNA


@pytest.mark.parametrize(
    "serialiser,deserialiser",
    (
        (json.dumps, json.loads),
        (pickle.dumps, pickle.loads),
        (lambda x: x, deserialise_object),
    ),
)
def test_deserialiser(serialiser, deserialiser):
    data = {"1": 1, "abc": [1, 2]}
    deserialised = io_app.from_primitive(deserialiser=deserialiser)
    assert deserialised(serialiser(data)) == data


@pytest.mark.parametrize(
    "data,dser", (([1, 2, 3], None), (DNA, io_app.from_primitive()))
)
def test_pickle_unpickle_apps(data, dser):
    pkld = io_app.to_primitive() + io_app.pickle_it()
    upkld = io_app.unpickle_it() + io_app.from_primitive()
    # need to add custom deserialiser for cogent3 objects
    # upkld = upkld if dser is None else upkld + dser
    assert upkld(pkld(data)) == data


def test_pickle_it_unpickleable():
    def foo():  # can't pickle a local function
        ...

    app = io_app.pickle_it()
    got = app(foo)
    assert isinstance(got, NotCompleted)


@pytest.mark.parametrize(
    "compress,decompress",
    ((bz2.compress, bz2.decompress), (gzip.compress, gzip.decompress)),
)
def test_compress_decompress(compress, decompress):
    data = pickle.dumps({"1": 1, "abc": [1, 2]})

    decompressor = io_app.decompress(decompressor=decompress)
    compressor = io_app.compress(compressor=compress)

    assert decompressor(compressor(data)) == data


@pytest.mark.parametrize("data", ([1, 2, 3], DNA)[1:])
def test_pickled_compress_roundtrip(data):
    serialised = io_app.to_primitive() + io_app.pickle_it() + io_app.compress()
    deserialised = io_app.decompress() + io_app.unpickle_it() + io_app.from_primitive()
    s = serialised(data)
    d = deserialised(s)
    assert d.label == data.label


# todo test objects where there is no unique_id provided, or inferrable,
# should return NotCompletedError


def test_write_db_load_db(fasta_dir, tmp_dir):
    from cogent3.app.sqlite_data_store import DataStoreSqlite
    from cogent3.util.misc import get_object_provenance

    orig_dstore = DataStoreDirectory(fasta_dir, suffix="fasta")
    data_store = DataStoreSqlite(tmp_dir / "test.sqlitedb", mode="w")
    load_seqs = io_app.load_unaligned()
    writer = io_app.write_db(data_store=data_store)
    reader = io_app.load_db()
    for m in orig_dstore:
        orig = load_seqs(m)
        m = writer(orig, identifier=m.unique_id)
        read = reader(m)
        assert orig == read

    assert data_store.record_type == get_object_provenance(orig)


def test_write_read_db_not_completed(tmp_dir):
    from cogent3.app.sqlite_data_store import DataStoreSqlite

    nc = NotCompleted("ERROR", "test", "for tracing", source="blah")
    data_store = DataStoreSqlite(tmp_dir / "test.sqlitedb", mode="w")
    writer = io_app.write_db(data_store=data_store)
    assert len(data_store.not_completed) == 0
    writer.main(nc, identifier="blah")
    assert len(data_store.not_completed) == 1
    reader = io_app.load_db()
    got = reader(data_store.not_completed[0])
    assert got.to_rich_dict() == nc.to_rich_dict()


def test_write_read_db_summary_not_completed(tmp_dir):
    # should not fail, even if not completed content is compressed
    from cogent3.app.sqlite_data_store import DataStoreSqlite

    nc = NotCompleted("ERROR", "test", "for tracing", source="blah")
    data_store = DataStoreSqlite(tmp_dir / "test.sqlitedb", mode="w")
    writer = get_app("write_db", data_store=data_store)
    writer.main(nc, identifier="blah")
    assert isinstance(writer.data_store.summary_not_completed, Table)


def test_write_db_parallel(tmp_dir, fasta_dir):
    """writing with overwrite in parallel should reset db"""
    dstore = open_data_store(fasta_dir, suffix="fasta")
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    members = [m for m in dstore if m.unique_id != "brca1.fasta"]
    reader = get_app("load_unaligned")
    aligner = get_app("align_to_ref")
    writer = get_app("write_seqs", out_dstore)
    process = reader + aligner + writer

    _ = process.apply_to(members, show_progress=False, parallel=True, cleanup=True)
    expect = [str(pathlib.Path(m.data_store.source) / m.unique_id) for m in out_dstore]

    # now get read only and check what's in there
    result = open_data_store(out_dstore.source, suffix="fasta")
    got = [str(pathlib.Path(m.data_store.source) / m.unique_id) for m in result]
    assert got != []
    assert got == expect


def test_define_data_store(fasta_dir):
    """returns an iterable data store"""
    found = open_data_store(fasta_dir, suffix=".fasta")
    assert len(found) > 1
    found = open_data_store(fasta_dir, suffix=".fasta", limit=2)
    assert len(found) == 2

    # and with a suffix
    found = list(open_data_store(fasta_dir, suffix=".fasta*"))
    assert len(found) > 2

    # with a wild-card suffix
    found = list(open_data_store(fasta_dir, suffix="*"))
    assert len(os.listdir(fasta_dir)) == len(found)

    # raises ValueError if suffix not provided or invalid
    with pytest.raises(ValueError):
        _ = open_data_store(fasta_dir)

    with pytest.raises(ValueError):
        _ = open_data_store(fasta_dir, 1)


def seqs():
    from cogent3 import make_unaligned_seqs

    return make_unaligned_seqs(
        data=dict(a="ACGG", b="GGC"),
        moltype="dna",
        info=dict(source="dummy/blah.1.2.fa"),
    )


def table():
    from cogent3 import make_table

    table = make_table(
        data=dict(a=[0, 1, 2], b=[0, 1, 2]),
    )
    table.source = "dummy/blah.1.2.fa"
    return table


def dir_dstore(tmp_dir):
    return open_data_store(tmp_dir / "sample", suffix="fa", mode="w")


def db_dstore(tmp_dir):
    return open_data_store(tmp_dir / "sample.sqlitedb", mode="w")


@pytest.mark.parametrize(
    "writer,data,dstore",
    (
        ("write_seqs", seqs(), dir_dstore),
        ("write_db", seqs(), db_dstore),
        ("write_json", seqs(), dir_dstore),
        ("write_tabular", table(), dir_dstore),
    ),
)
def test_writer_unique_id_arg(tmp_dir, writer, data, dstore):
    def uniqid(source):
        from cogent3.app.data_store import get_data_source

        name = pathlib.Path(get_data_source(source)).name
        return name.split(".")[0]

    writer = get_app(writer, data_store=dstore(tmp_dir), id_from_source=uniqid)
    m = writer(data)
    # directory data stores have a suffix, so we create expected name using that
    suffix = getattr(writer.data_store, "suffix", "")
    suffix = f".{suffix}" if suffix else suffix
    expect = f"blah{suffix}"
    assert m.unique_id == expect


@pytest.mark.parametrize(
    "writer",
    (
        "write_seqs",
        "write_db",
        "write_json",
        "write_tabular",
    ),
)
def test_writer_fails_on_data_store(writer):
    # should raise a type error if not a data store provided
    with pytest.raises(TypeError):
        get_app(writer, data_store="not-a-data_store")


def test_open_suffix_dirname(tmp_dir):
    # open_data_store correctly identifies this as a directory data store
    outpath = tmp_dir / "melsubgroup_aln_flydivas_v1.2"
    outpath.mkdir(exist_ok=True)
    dstore = open_data_store(outpath, suffix="txt")
    assert isinstance(dstore, DataStoreDirectory)


@pytest.mark.parametrize("data", ({"a": [0, 1]}, DNA))
def test_default_serialiser_deserialiser(data):
    # the default deserialiser should successfully reverse the
    # default serialiser
    s = DEFAULT_SERIALISER(data)
    ds = DEFAULT_DESERIALISER(s)
    assert ds == data


def test_to_json():
    to_j = io_app.to_json()
    data = {"a": [0, 1]}
    assert to_j(data) == json.dumps(data)


def test_from_json():
    from_j = io_app.from_json()
    assert from_j('{"a": [0, 1]}') == {"a": [0, 1]}


def test_to_from_json():
    to_j = io_app.to_json()
    from_j = io_app.from_json()
    app = to_j + from_j
    data = {"a": [0, 1]}
    assert app(data) == data
    assert app(data) is not data


def test_to_json_combines():
    app = io_app.to_primitive() + io_app.to_json()
    assert app(DNA) == DNA.to_json()


def test_open_zipped(zipped_full):
    got = open_data_store(zipped_full.source, mode="r", suffix="fasta")
    assert len(got) == len(zipped_full)
    assert isinstance(got, type(zipped_full))


@pytest.fixture(scope="function")
def relpath(DATA_DIR):
    # express the data path as relative to user home
    # have to make a tempdir for this to work in github actions
    data = (DATA_DIR / "brca1_5.paml").read_text()
    user = pathlib.Path("~")
    with tempfile.TemporaryDirectory(dir=user.expanduser()) as dirname:
        out_path = pathlib.Path(dirname) / "brca1_5.paml"
        out_path.write_text(data)
        yield str(user / out_path.parent.name / "brca1_5.paml")


@pytest.mark.parametrize("type_", (str, pathlib.Path))
def test_expand_user(relpath, type_):
    loader = get_app("load_aligned", format="paml")
    # define path using the "~" prefix
    seqs = loader(type_(relpath))
    assert isinstance(seqs, ArrayAlignment)

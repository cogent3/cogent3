import json
import pathlib
import shutil

from itertools import product
from pathlib import Path
from pickle import dumps, loads

import pytest

from scitrack import get_text_hexdigest

import cogent3.app.io as io_app

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store import (
    _MD5_TABLE,
    _NOT_COMPLETED_TABLE,
    APPEND,
    OVERWRITE,
    READONLY,
    DataStoreDirectory,
    ReadOnlyDataStoreZipped,
    convert_directory_datastore,
    convert_tinydb_to_sqlite,
    get_data_source,
    get_unique_id,
    load_record_from_json,
    summary_not_completeds,
)
from cogent3.util.table import Table
from cogent3.util.union_dict import UnionDict


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory):
    return Path(tmp_path_factory.mktemp("datastore"))


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    # this set's the working directory for all tests in this module
    # as a tmp dir
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="function")
def fasta_dir(DATA_DIR, tmp_dir):
    tmp_dir = Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    fasta_dir = tmp_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = fasta_dir / fn.name
        dest.write_text(fn.read_text())
    return fasta_dir


@pytest.fixture(scope="function")
def write_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    write_dir = tmp_dir / "write"
    write_dir.mkdir(parents=True, exist_ok=True)
    yield write_dir
    shutil.rmtree(write_dir)


@pytest.fixture(scope="function")
def nc_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    nc_dir = tmp_dir / "nc_dir"
    nc_dir.mkdir(parents=True, exist_ok=True)
    yield nc_dir
    shutil.rmtree(nc_dir)


@pytest.fixture(scope="function")
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix="fasta", mode=READONLY)


@pytest.fixture(scope="function")
def w_dstore(write_dir):
    return DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)


@pytest.fixture(scope="function")
def nc_dstore(DATA_DIR, nc_dir):
    dstore = DataStoreDirectory(nc_dir, suffix="fasta", mode=OVERWRITE)
    # write one log file
    log_filename = "scitrack.log"
    dstore.write_log(unique_id=log_filename, data=(DATA_DIR / log_filename).read_text())
    # write three not_completed file
    nc = [
        NotCompleted(
            "FAIL", f"dummy{i}", f"dummy_message{i}", source=f"dummy_source{i}"
        )
        for i in range(3)
    ]
    for i, item in enumerate(nc):
        dstore.write_not_completed(unique_id=f"nc{i + 1}", data=item.to_json())
    assert len(dstore.not_completed) == 3
    assert len(list((nc_dir / _MD5_TABLE).glob("*.txt"))) == len(dstore)
    filenames = DATA_DIR.glob("*.fasta")
    # write six fasta file
    for fn in filenames:
        identifier = fn.name
        dstore.write(unique_id=identifier, data=fn.read_text())
    return dstore


@pytest.fixture(scope="function")
def completed_objects(ro_dstore):
    return {f"{Path(m.unique_id).stem}": m.read() for m in ro_dstore}


@pytest.fixture(scope="function")
def nc_objects():
    return {
        f"id_{i}": NotCompleted("ERROR", "location", "message", source=f"id_{i}")
        for i in range(3)
    }


@pytest.fixture(scope="function")
def Sample_oldDirectoryDataStore(DATA_DIR, tmp_dir):
    tmp_dir = Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    old_dir = tmp_dir / "old_dir"
    old_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = old_dir / fn.name
        dest.write_text(fn.read_text())
    return old_dir


@pytest.fixture(scope="session")
def log_data(DATA_DIR):
    path = DATA_DIR / "scitrack.log"
    return path.read_text()


@pytest.fixture(scope="function")
def full_dstore(write_dir, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())

    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)

    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


@pytest.fixture(scope="function")
def tinydbfile_locked(DATA_DIR, tmp_dir):
    path = tmp_dir / "sample_locked.tinydb"
    shutil.copy(DATA_DIR / path.name, path)
    return Path(path)


@pytest.fixture(scope="function")
def tinydbfile_unlocked(DATA_DIR, tmp_dir):
    try:
        from tinydb import Query, TinyDB
        from tinydb.middlewares import CachingMiddleware
        from tinydb.storages import JSONStorage
    except ImportError as e:
        raise ImportError(
            "You need to install tinydb to be able to migrate to new datastore."
        ) from e

    locked_path = tmp_dir / "sample_locked.tinydb"
    shutil.copy(DATA_DIR / locked_path.name, locked_path)
    unlocked_path = locked_path.parent / "sample_unlocked.tinydb"
    locked_path.rename(unlocked_path)

    storage = CachingMiddleware(JSONStorage)
    db = TinyDB(unlocked_path, storage=storage)
    query = Query().identifier.matches("LOCK")
    db.remove(query)
    db.storage.flush()
    db.close()
    return unlocked_path


def test_data_member_eq(ro_dstore, fasta_dir):
    ro_dstore2 = DataStoreDirectory(fasta_dir, mode="r", suffix="fasta")
    name = "brca1.fasta"
    mem1 = [m for m in ro_dstore.completed if m.unique_id == name][0]
    mem2 = [m for m in ro_dstore2.completed if m.unique_id == name][0]
    assert mem1 != mem2


@pytest.mark.parametrize("dest", (None, "mydata1.sqlitedb"))
def test_convert_tinydb_to_sqlite(tmp_dir, dest, tinydbfile_locked):
    dest = dest if dest is None else tmp_dir / dest
    if dest:
        dest.unlink(missing_ok=True)
    else:
        (Path(tinydbfile_locked.parent) / f"{tinydbfile_locked.stem}.sqlitedb").unlink(
            missing_ok=True
        )
    dstore_sqlite = convert_tinydb_to_sqlite(tinydbfile_locked, dest=dest)
    assert len(dstore_sqlite) == 6
    # tinydb has hard-coded value of lock 123
    assert dstore_sqlite._lock_id == 123


@pytest.mark.parametrize("dest", (None, "mydata2.sqlitedb"))
def test_convert_tinydb_unlocked_to_sqlite(tmp_dir, dest, tinydbfile_unlocked):
    dest = dest if dest is None else tmp_dir / dest
    if dest:
        dest.unlink(missing_ok=True)
    else:
        (
            Path(tinydbfile_unlocked.parent) / f"{tinydbfile_unlocked.stem}.sqlitedb"
        ).unlink(missing_ok=True)
    dstore_sqlite = convert_tinydb_to_sqlite(tinydbfile_unlocked, dest=dest)
    assert len(dstore_sqlite) == 6
    assert dstore_sqlite._lock_id == None


def test_convert_tinydb_to_sqlite_error(DATA_DIR, tmp_dir):
    path = tmp_dir / "sample_locked.tinydb"
    shutil.copy(DATA_DIR / path.name, path)
    dest = tmp_dir / "data1.sqlitedb"
    # create a file at dest path to trigger error
    dest.write_text("blah")
    with pytest.raises(IOError):
        _ = convert_tinydb_to_sqlite(path, dest=dest)


@pytest.mark.parametrize(
    "orig,num_logs",
    (
        ("sample_locked.tinydb", 1),
        ("sample_locked_w_log.tinydb", 2),
    ),
)
def test_convert_tinydbs_to_sqlite(DATA_DIR, tmp_dir, orig, num_logs):
    orig = DATA_DIR / orig
    src = tmp_dir / orig.name
    shutil.copy(orig, src)
    dest = tmp_dir / "data1.sqlitedb"
    got = convert_tinydb_to_sqlite(src, dest=dest)
    assert len(got.logs) == num_logs


def test_convert_directory_datastore(Sample_oldDirectoryDataStore, write_dir):
    new_dstore = convert_directory_datastore(
        Sample_oldDirectoryDataStore, write_dir, ".fasta"
    )
    assert len(new_dstore) == 6


def test_fail_try_append(full_dstore, completed_objects):
    full_dstore._mode = APPEND
    id, data = list(completed_objects.items())[0]
    with pytest.raises(IOError):
        full_dstore.write(unique_id=id, data=data)


def test_contains(ro_dstore):
    """correctly identify when a data store contains a member"""
    assert "brca1.fasta" in ro_dstore
    assert "brca1" in ro_dstore


def test_len(ro_dstore):
    """DataStore returns correct len"""
    expect = len(list(ro_dstore.source.glob("*.fasta")))
    assert expect == len(ro_dstore) == len(ro_dstore.members)


def test_getitem(ro_dstore):
    with pytest.raises(IndexError):
        _ = ro_dstore[len(ro_dstore)]

    last = ro_dstore[-1]
    first = ro_dstore[0]
    assert last.unique_id != first.unique_id


def test_iterall(ro_dstore):
    expect = {fn.name for fn in ro_dstore.source.glob("*.fasta")}
    got = {m.unique_id for m in ro_dstore}
    assert expect == got


def test_read(ro_dstore):
    """correctly read content"""
    expect = (ro_dstore.source / "brca1.fasta").read_text()
    got = ro_dstore.read("brca1.fasta")
    assert got == expect


def test_pickleable_roundtrip(ro_dstore):
    """pickling of data stores should be reversible"""
    re_dstore = loads(dumps(ro_dstore))
    assert str(ro_dstore) == str(re_dstore)
    assert ro_dstore[0].read() == re_dstore[0].read()


def test_pickleable_member_roundtrip(ro_dstore):
    """pickling of data store members should be reversible"""
    re_member = loads(dumps(ro_dstore[0]))
    data = re_member.read()
    assert len(data) > 0


def test_empty_directory(fasta_dir):
    dstore = DataStoreDirectory(fasta_dir, suffix=".txt")
    assert 0 == len(dstore)


def test_no_logs(ro_dstore):
    assert len(ro_dstore.logs) == 0


def test_no_not_completed(ro_dstore):
    assert len(ro_dstore.not_completed) == 0


def test_logs(nc_dstore):
    assert len(nc_dstore.logs) == 1
    log = nc_dstore.logs[0].read()
    assert isinstance(log, str)


def test_not_completed(nc_dstore):
    assert len(nc_dstore.not_completed) == 3
    nc = nc_dstore.not_completed[0].read()
    assert isinstance(nc, str)


def test_drop_not_completed(nc_dstore):
    num_completed = len(nc_dstore.completed)
    num_not_completed = len(nc_dstore.not_completed)
    num_md5 = len(list((nc_dstore.source / _MD5_TABLE).glob("*.txt")))
    assert num_not_completed == 3
    assert num_completed == 6
    assert len(nc_dstore) == 9
    assert num_md5 == num_completed + num_not_completed
    nc_dstore.drop_not_completed()
    assert len(nc_dstore.not_completed) == 0
    num_md5 = len(list((nc_dstore.source / _MD5_TABLE).glob("*.txt")))
    assert num_md5 == num_completed


def test_write_read_only_datastore(ro_dstore):
    with pytest.raises(IOError):
        ro_dstore.write(unique_id="brca1.fasta", data="test data")


def test_write(fasta_dir, w_dstore):
    """correctly write content"""
    expect = Path(fasta_dir / "brca1.fasta").read_text()
    identifier = "brca1.fasta"
    w_dstore.write(unique_id=identifier, data=expect)
    got = w_dstore.read(identifier)
    assert got == expect


def test_multi_write(fasta_dir, w_dstore):
    """correctly write multiple files to data store"""
    expect_a = Path(fasta_dir / "brca1.fasta").read_text()
    expect_b = Path(fasta_dir / "primates_brca1.fasta").read_text()
    identifier_a = "brca2.fasta"
    identifier_b = "primates_brca2.fasta"
    w_dstore.write(unique_id=identifier_a, data=expect_a)
    w_dstore.write(unique_id=identifier_b, data=expect_b)
    got_a = w_dstore.read(identifier_a)
    got_b = w_dstore.read(identifier_b)
    # check that both bits of data match
    assert got_a == expect_a
    assert got_b == expect_b


def test_append(w_dstore):
    """correctly write content"""
    identifier = "test1.fasta"
    data = "test data"
    w_dstore.write(unique_id=identifier, data=data)
    got = w_dstore.read(identifier)
    assert got == data


def test_no_not_completed_subdir(nc_dstore):
    expect = f"{len(nc_dstore.completed)+len(nc_dstore.not_completed)}x member"
    assert str(nc_dstore).startswith(expect)
    # first remove not_completed directory
    nc_dstore.drop_not_completed()
    # test repr work without not_completed directory
    assert not Path(nc_dstore.source / _NOT_COMPLETED_TABLE).exists()
    expect = f"{len(nc_dstore.completed)}x member"
    assert str(nc_dstore).startswith(expect)
    expect = f"{len(nc_dstore)}x member"
    assert str(nc_dstore).startswith(expect)
    assert len(nc_dstore) == len(nc_dstore.completed)
    not_dir = nc_dstore.source / _NOT_COMPLETED_TABLE
    not_dir.mkdir(exist_ok=True)


def test_limit_datastore(nc_dstore):  # new changed
    assert len(nc_dstore) == len(nc_dstore.completed) + len(nc_dstore.not_completed)
    nc_dstore._limit = len(nc_dstore.completed) // 2
    nc_dstore._completed = []
    nc_dstore._not_completed = []
    assert len(nc_dstore.completed) == len(nc_dstore.not_completed) == nc_dstore.limit
    assert len(nc_dstore) == len(nc_dstore.completed) + len(nc_dstore.not_completed)
    nc_dstore.drop_not_completed()
    assert len(nc_dstore) == len(nc_dstore.completed)
    assert len(nc_dstore.not_completed) == 0
    nc_dstore._limit = len(nc_dstore.completed) // 2
    nc_dstore._completed = []
    nc_dstore._not_completed = []
    assert len(nc_dstore) == len(nc_dstore.completed) == nc_dstore.limit
    assert len(nc_dstore.not_completed) == 0


def test_md5_sum(nc_dstore):
    for m in nc_dstore.members:
        data = m.read()
        md5 = nc_dstore.md5(m.unique_id)
        assert md5 == get_text_hexdigest(data)


def test_md5_none(fasta_dir):
    dstore = DataStoreDirectory(fasta_dir, suffix="fasta")
    for m in dstore.members:
        assert m.md5 is None


def test_md5_missing(nc_dstore):
    nc_dstore.md5("unknown") is None


def test_summary_logs_missing_field(nc_dstore):
    log_path = Path(nc_dstore.source) / nc_dstore.logs[0].unique_id
    data = [
        l for l in log_path.read_text().splitlines() if "composable function" not in l
    ]
    log_path.write_text("\n".join(data))
    assert isinstance(nc_dstore.summary_logs, Table)


@pytest.mark.parametrize("use_dser", (False, True))
def test_summary_not_completed_func(nc_objects, use_dser):
    dstore = io_app.open_data_store(":memory:", mode="w")
    writer = io_app.write_db(dstore)
    deser = io_app.load_db().deserialiser if use_dser else None
    for nc in nc_objects:
        writer(nc)

    got = summary_not_completeds(dstore, deserialise=deser)
    assert isinstance(got, Table)
    if use_dser:
        assert got.shape[0] >= 1
    else:
        assert got.shape[0] == 0


def test_write_if_member_exists(full_dstore, write_dir):  # new changes
    """correctly write content"""
    expect = Path(write_dir / "brca1.fasta").read_text()
    identifier = "brca1.fasta"
    len_dstore = len(full_dstore)
    full_dstore.write(unique_id=identifier, data=expect)
    assert len_dstore == len(full_dstore)
    got = full_dstore.read(identifier)
    assert got == expect
    full_dstore._mode = OVERWRITE
    full_dstore.write(unique_id=identifier, data=expect)
    assert len_dstore == len(full_dstore)
    got = full_dstore.read(identifier)
    assert got == expect


def test_write_success_replaces_not_completed(full_dstore):
    """correctly write content"""
    nc = full_dstore.not_completed[0].unique_id
    data = full_dstore.completed[0].read()
    new_id = Path(nc.replace(".json", f".{full_dstore.suffix}")).name
    num = len(full_dstore)
    full_dstore.write(unique_id=new_id, data=data)
    assert len(full_dstore) == num


@pytest.mark.parametrize("klass", (str, Path))
def test_get_data_source_attr(klass):
    """handles case where input has source attribute string object or pathlib object"""

    class dummy:
        source = None

    obj = dummy()
    value = klass("some/path.txt")
    obj.source = value
    got = get_data_source(obj)
    assert got == str("path.txt")


_types = tuple(product((dict, UnionDict), (str, Path)))


@pytest.mark.parametrize("container_type,source_stype", _types)
def test_get_data_source_dict(container_type, source_stype):
    """handles case where input is dict (sub)class instance with top level source key"""
    value = source_stype("some/path.txt")
    data = container_type(source=value)
    got = get_data_source(data)
    assert got == "path.txt"


@pytest.mark.parametrize(
    "name", ("path/name.txt", "path/name.gz", "path/name.fasta.gz", "name.fasta.gz")
)
def test_get_unique_id(name):
    got = get_unique_id(name)
    assert got == "name"


@pytest.mark.parametrize("data", (dict(), set(), dict(info=dict())))
def test_get_data_source_none(data):
    assert get_data_source(data) is None


@pytest.mark.parametrize("klass", (str, Path))
def test_get_data_source_seqcoll(klass):
    """handles case where input is sequence collection object"""
    from cogent3 import make_unaligned_seqs

    value = klass("some/path.txt")
    obj = make_unaligned_seqs(
        data=dict(seq1="ACGG"), moltype="dna", info=dict(source=value, random_key=1234)
    )
    got = get_data_source(obj)
    assert got == "path.txt"


def test_load_record_from_json():
    """handle different types of input"""
    orig = {"data": "blah", "identifier": "some.json", "completed": True}
    data = orig.copy()
    data2 = data.copy()
    data2["data"] = json.dumps(data)
    for d in (data, json.dumps(data), data2):
        expected = "blah" if d != data2 else json.loads(data2["data"])
        Id, data_, compl = load_record_from_json(d)
        assert Id == "some.json"
        assert data_ == expected
        assert compl == True


def test_write_read_not_completed(nc_dstore):
    nc_dstore.drop_not_completed()
    assert len(nc_dstore.not_completed) == 0
    nc = NotCompleted("ERROR", "test", "for tracing", source="blah")
    writer = io_app.write_seqs(data_store=nc_dstore)
    writer.main(nc, identifier="blah")
    assert len(nc_dstore.not_completed) == 1
    got = nc_dstore.not_completed[0].read()
    assert nc.to_json() == got


@pytest.fixture
def zipped_basic(fasta_dir):
    # converts the fasta_dir into a zipped archive

    path = shutil.make_archive(
        base_name=fasta_dir.name,
        format="zip",
        base_dir=fasta_dir.name,
        root_dir=fasta_dir.parent,
    )
    return pathlib.Path(path)


@pytest.fixture
def zipped_full(full_dstore):
    # converts the fasta_dir into a zipped archive
    source = pathlib.Path(full_dstore.source)
    path = shutil.make_archive(
        base_name=source.name,
        format="zip",
        base_dir=source.name,
        root_dir=source.parent,
    )
    return ReadOnlyDataStoreZipped(pathlib.Path(path), suffix="fasta")


def test_zipped_ro_fail(zipped_basic):
    with pytest.raises(ValueError):
        ReadOnlyDataStoreZipped(zipped_basic, suffix="fasta", mode="w")


def test_zipped_ro_ioerror():
    with pytest.raises(IOError):
        ReadOnlyDataStoreZipped("blah-1234.zip", suffix="fasta")


def _get_member_data(members) -> dict:
    return {m.unique_id: m.read().splitlines() for m in members}


def test_zipped_ro_basic(zipped_basic, ro_dstore):
    dstore = ReadOnlyDataStoreZipped(zipped_basic, suffix="fasta")
    assert len(dstore.completed) == len(ro_dstore.completed)
    assert len(dstore) == len(ro_dstore)
    expect = _get_member_data(ro_dstore.completed)
    got = _get_member_data(dstore.completed)
    assert expect == got
    # there are no not-completeds in a basic dstore
    expect = _get_member_data(ro_dstore.not_completed)
    got = _get_member_data(dstore.not_completed)
    assert expect == got


def test_zipped_ro_full(zipped_full, full_dstore):
    got_ids = {m.unique_id for m in zipped_full.completed}
    expect_ids = {m.unique_id for m in full_dstore.completed}
    assert got_ids == expect_ids

    got_ids = {m.unique_id for m in zipped_full.not_completed}
    expect_ids = {m.unique_id for m in full_dstore.not_completed}
    assert got_ids == expect_ids
    assert len(zipped_full) == len(full_dstore)

    expect = _get_member_data(full_dstore.completed)
    got = _get_member_data(zipped_full.completed)
    assert expect == got
    # there are not-completeds in the full dstore
    expect = _get_member_data(full_dstore.not_completed)
    got = _get_member_data(zipped_full.not_completed)
    assert expect == got


def test_zipped_logs(zipped_full, full_dstore):
    assert len(zipped_full.logs) == len(full_dstore.logs)
    expect = _get_member_data(full_dstore.logs)
    got = _get_member_data(zipped_full.logs)
    assert expect == got


def test_zipped_md5(zipped_full, full_dstore):
    expect = {m.unique_id: full_dstore.md5(m.unique_id) for m in full_dstore.completed}
    got = {m.unique_id: zipped_full.md5(m.unique_id) for m in zipped_full.completed}
    assert got == expect


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_zipped_contains(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    assert "formattest.fasta" in full_dstore


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_members(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    completed = full_dstore.completed
    not_completed = full_dstore.not_completed
    assert full_dstore.members == completed + not_completed


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_describe(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    got = full_dstore.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_iter(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    members = list(full_dstore)
    assert members == full_dstore.members


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_validate(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    validation_table = full_dstore.validate()
    assert isinstance(validation_table, Table)
    assert validation_table["Num md5sum correct", "Value"] == len(full_dstore)
    assert validation_table["Num md5sum incorrect", "Value"] == 0
    assert validation_table["Has log", "Value"] == True


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_summary_logs(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    # log summary has a row per log file and a column for each property
    got = full_dstore.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


@pytest.mark.parametrize("dstore", ("full_dstore", "zipped_full"))
def test_summary_not_completed(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    got = full_dstore.summary_not_completed
    assert got.shape >= (1, 1)
    assert isinstance(got, Table)


@pytest.fixture(scope="function")
def app_dstore_in(tmp_path):
    from cogent3 import get_app, open_data_store

    in_path = tmp_path / "in_data"
    in_path.mkdir(parents=True)
    fasta_content = ">seq\nACGT"
    with open(in_path / "one.fa", "w") as file:
        file.write(fasta_content)

    dstore_in = open_data_store(in_path, suffix=".fa", mode="r")
    dstore_out = open_data_store(tmp_path / "data_out", suffix="fa", mode="w")
    loader = get_app("load_unaligned")
    writer = get_app("write_seqs", dstore_out)

    pipe = loader + writer
    return pipe, dstore_in


def test_write_multiple_times_apply_to(app_dstore_in):
    app, dstore_in = app_dstore_in
    app.apply_to(dstore_in)
    orig_length = len(app.data_store)
    app.apply_to(dstore_in)
    assert len(app.data_store) == orig_length


def test_directory_data_store_write_compressed(tmp_path):
    from cogent3 import get_app, make_aligned_seqs, open_data_store

    out = open_data_store(base_path=tmp_path / "demo", suffix="fa.gz", mode="w")
    writer = get_app("write_seqs", data_store=out)
    seqs = make_aligned_seqs(
        dict(s1="CG--T", s2="CGTTT"), moltype="dna", info=dict(source="test")
    )
    got = writer(seqs)  # pylint: disable=not-callable
    assert got, got

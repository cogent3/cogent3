import shutil

from pathlib import Path
from pickle import dumps, loads

import pytest

from scitrack import get_text_hexdigest

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store import (
    ReadOnlyDirectoryDataStore,
    ReadOnlyTinyDbDataStore,
    ReadOnlyZippedDataStore,
    SingleReadDataStore,
    WritableDirectoryDataStore,
    WritableTinyDbDataStore,
)
from cogent3.app.data_store_new import (
    _MD5_TABLE,
    _NOT_COMPLETED_TABLE,
    APPEND,
    OVERWRITE,
    READONLY,
    DataMember,
    DataStoreDirectory,
    _legacy_tiny_reader,
    upgrade_data_store,
)
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

DATA_DIR = Path(__file__).parent.parent / "data"


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("datastore")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    # this set's the working directory for all tests in this module
    # as a tmp dir
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="session")
def fasta_dir(tmp_dir):
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


@pytest.fixture(scope="session")
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix="fasta", mode=READONLY)


@pytest.fixture(scope="function")
def w_dstore(write_dir):
    return DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)


@pytest.fixture(scope="function")
def nc_dstore(nc_dir):
    dstore = DataStoreDirectory(nc_dir, suffix="fasta", mode=OVERWRITE)
    # write one log file
    log_filename = "scitrack.log"
    dstore.write_log(log_filename, (DATA_DIR / log_filename).read_text())
    # write three not_completed file
    nc = [
        NotCompleted(
            "FAIL", f"dummy{i}", f"dummy_message{i}", source=f"dummy_source{i}"
        )
        for i in range(3)
    ]
    for i, item in enumerate(nc):
        dstore.write_not_completed(f"nc{i + 1}", item.to_json())
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
def log_data():
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


def test_convert_tinydb_to_Sqlitedb(fasta_dir):
    path = DATA_DIR / "data1.tinydb"
    dstore_sqlite = _legacy_tiny_reader(path, ".tinydb")
    assert len(dstore_sqlite) == 6


def test_upgrade_dstore(fasta_dir, write_dir):
    dstore = ReadOnlyDirectoryDataStore(fasta_dir, suffix=".fasta")
    new_dstore = upgrade_data_store(dstore.source, write_dir, ".fasta", ".fasta")
    assert len(dstore) == len(new_dstore)


def test_fail_try_append(full_dstore, completed_objects):
    full_dstore.mode = APPEND
    id, data = list(completed_objects.items())[0]
    with pytest.raises(IOError):
        full_dstore.write(unique_id=id, data=data)


def test_contains(ro_dstore):
    """correctly identify when a data store contains a member"""
    assert "brca1.fasta" in ro_dstore
    assert "brca1" in ro_dstore


def test_iter(ro_dstore):
    """DataStore objects allow iteration over members"""
    members = list(ro_dstore)
    assert members == ro_dstore.members


def test_members(ro_dstore):
    for m in ro_dstore:
        assert isinstance(m, DataMember)
    assert len(ro_dstore) == len(ro_dstore.completed) + len(ro_dstore.not_completed)


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


def test_write_read_only_member(ro_dstore):
    with pytest.raises(IOError):
        ro_dstore[0].write("test data")


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
    assert repr(nc_dstore).startswith(expect)
    # first remove not_completed directory
    nc_dstore.drop_not_completed()
    # test repr work without not_completed directory
    assert not Path(nc_dstore.source / _NOT_COMPLETED_TABLE).exists()
    expect = f"{len(nc_dstore.completed)}x member"
    assert repr(nc_dstore).startswith(expect)
    expect = f"{len(nc_dstore)}x member"
    assert repr(nc_dstore).startswith(expect)
    assert len(nc_dstore) == len(nc_dstore.completed)
    not_dir = nc_dstore.source / _NOT_COMPLETED_TABLE
    not_dir.mkdir(exist_ok=True)


def test_limit_datastore(nc_dstore):  # new changed
    assert len(nc_dstore) == len(nc_dstore.completed) + len(nc_dstore.not_completed)
    nc_dstore.limit = len(nc_dstore.completed) // 2
    nc_dstore._completed = []
    nc_dstore._not_completed = []
    assert len(nc_dstore.completed) == len(nc_dstore.not_completed) == nc_dstore.limit
    assert len(nc_dstore) == len(nc_dstore.completed) + len(nc_dstore.not_completed)
    nc_dstore.drop_not_completed()
    assert len(nc_dstore) == len(nc_dstore.completed)
    assert len(nc_dstore.not_completed) == 0
    nc_dstore.limit = len(nc_dstore.completed) // 2
    nc_dstore._completed = []
    nc_dstore._not_completed = []
    assert len(nc_dstore) == len(nc_dstore.completed) == nc_dstore.limit
    assert len(nc_dstore.not_completed) == 0


def test_md5_sum(nc_dstore):
    for m in nc_dstore.members:
        data = m.read()
        md5 = nc_dstore.md5(m.unique_id)
        assert md5 == get_text_hexdigest(data)


def test_summary_logs(full_dstore):
    # log summary has a row per log file and a column for each property
    got = full_dstore.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


def test_summary_not_completed(full_dstore):
    got = full_dstore.summary_not_completed
    assert got.shape >= (1, 1)
    assert isinstance(got, Table)


def test_describe(full_dstore):
    got = full_dstore.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


def test_iter(full_dstore):
    members = list(full_dstore)
    assert members == full_dstore.members


def test_members(full_dstore):
    completed = full_dstore.completed
    not_completed = full_dstore.not_completed
    assert full_dstore.members == completed + not_completed


def test_validate(full_dstore):
    validation_table = full_dstore.validate()
    assert isinstance(validation_table, Table)
    assert validation_table["Num md5sum correct", "Value"] == len(full_dstore)
    assert validation_table["Num md5sum incorrect", "Value"] == 0
    assert validation_table["Has log", "Value"] == True


def test_write_if_member_exists(full_dstore, write_dir): #new changes
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

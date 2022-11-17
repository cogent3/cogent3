from datetime import datetime
from pathlib import Path
from pickle import dumps, loads

import pytest

from scitrack import get_text_hexdigest

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store_new import (
    OVERWRITE,
    READONLY,
    DataMemberABC,
    DataStoreDirectory,
)
from cogent3.app.sqlite_data_store import (
    _LOG_TABLE,
    _RESULT_TABLE,
    SqliteDataStore,
    open_sqlite_db_rw,
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


@pytest.fixture(scope="function")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("sqldb")


@pytest.fixture(scope="function")
def db_dir(tmp_dir):
    db_dir = tmp_dir / "data.sqlitedb"
    return db_dir


@pytest.fixture(scope="function")
def ro_dir_dstore():
    return DataStoreDirectory(DATA_DIR, suffix="fasta")


@pytest.fixture(scope="function")
def sql_dstore(ro_dir_dstore, db_dir):
    # we now need to write these out to a path
    sql_dstore = SqliteDataStore(db_dir, mode=OVERWRITE)
    for m in ro_dir_dstore:
        sql_dstore.write(data=m.read(), unique_id=m.unique_id)
    return sql_dstore


@pytest.fixture(scope="function")
def ro_sql_dstore(sql_dstore):
    # we now need to write these out to a path
    sql_dstore._mode = READONLY
    return sql_dstore


@pytest.fixture(scope="function")
def completed_objects(ro_dir_dstore):
    return {f"{Path(m.unique_id).stem}": m.read() for m in ro_dir_dstore}


@pytest.fixture(scope="function")
def nc_objects():
    return {
        f"id_{i}": NotCompleted("ERROR", "location", "message", source=f"id_{i}")
        for i in range(3)
    }


@pytest.fixture(scope="function")
def log_data():
    path = Path(__file__).parent.parent / "data" / "scitrack.log"
    return path.read_text()


@pytest.fixture(scope="function")
def full_dstore_directory(db_dir, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(db_dir, suffix="fasta", mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())

    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)

    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


@pytest.fixture(scope="function")
def full_dstore_sqlite(db_dir, nc_objects, completed_objects, log_data):
    dstore = SqliteDataStore(db_dir, mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())
    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)
    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


@pytest.fixture(scope="function")
def dstore_on_disk(full_dstore_sqlite):
    path = full_dstore_sqlite.source
    full_dstore_sqlite.close()
    return path


def test_open_existing(dstore_on_disk):
    ro = SqliteDataStore(dstore_on_disk, mode=READONLY)
    assert len(ro) > 0
    assert len(ro.completed) > 0
    assert len(ro.not_completed) > 0
    assert len(ro.logs) > 0
    for attr in ("summary_logs", "summary_not_completed", "describe"):
        assert isinstance(getattr(ro, attr), Table)


def test_open_to_append(dstore_on_disk):
    ro = SqliteDataStore(dstore_on_disk, mode=OVERWRITE)


def test_db_creation():
    db = SqliteDataStore(":memory:", mode=OVERWRITE)
    db = db.db
    result = db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    assert len(result) == 3
    created_names = {r["name"] for r in result}
    assert created_names == {
        _LOG_TABLE,
        _RESULT_TABLE,
        "state",
    }
    rows = db.execute(f"Select * from {_LOG_TABLE}").fetchall()
    assert len(rows) == 0


def test_db_init_log():
    dstore = SqliteDataStore(":memory:", mode=OVERWRITE)
    dstore._init_log()
    rows = dstore.db.execute(f"Select * from {_LOG_TABLE}").fetchall()
    assert len(rows) == 1
    assert rows[0]["date"].date() == datetime.today().date()


def test_open_sqlite_db_rw():
    db = open_sqlite_db_rw(":memory:")
    # should make tables for completed, not completed, md5, log and state
    result = db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    assert len(result) == 3
    created_names = {r["name"] for r in result}
    assert created_names == {
        _LOG_TABLE,
        _RESULT_TABLE,
        "state",
    }


def test_rw_sql_dstore_mem(completed_objects):
    """in memory dstore"""
    db = SqliteDataStore(":memory:", mode=OVERWRITE)
    for unique_id, obj in completed_objects.items():
        db.write(data=obj, unique_id=unique_id)
    expect = len(completed_objects)
    query = f"SELECT count(*) as c FROM {_RESULT_TABLE} WHERE not_completed=?"
    got = db.db.execute(query, (0,)).fetchone()["c"]
    assert got == expect, f"Failed for {_RESULT_TABLE} number of rows"
    assert len(db.completed) == expect


def test_not_completed(nc_objects):
    db = SqliteDataStore(":memory:", mode=OVERWRITE)
    for unique_id, obj in nc_objects.items():
        db.write_not_completed(data=obj.to_json(), unique_id=unique_id)
    expect = len(nc_objects)
    query = f"SELECT count(*) as c FROM {_RESULT_TABLE} WHERE not_completed=?"
    got = db.db.execute(query, (1,)).fetchone()["c"]
    assert got == expect, f"Failed for {_RESULT_TABLE} number of rows"
    assert len(db.not_completed) == expect


def test_logdata(log_data):
    db = SqliteDataStore(":memory:", mode=OVERWRITE)
    db.write_log(data=log_data, unique_id="scitrack.log")
    query = f"select count(*) as c from {_LOG_TABLE}"
    got = db.db.execute(query).fetchone()["c"]
    assert got == 1


def test_drop_not_completed(nc_objects):
    db = SqliteDataStore(":memory:", mode=OVERWRITE)
    for unique_id, obj in nc_objects.items():
        db.write_not_completed(data=obj.to_json(), unique_id=unique_id)
    expect = len(nc_objects)
    assert len(db.not_completed) == expect
    db.drop_not_completed()
    assert len(db.not_completed) == 0


def test_contains(sql_dstore):
    """correctly identify when a data store contains a member"""
    assert "brca1.fasta" in sql_dstore


def test_iter(sql_dstore):
    """DataStore objects allow iteration over members"""
    members = list(sql_dstore)
    assert members == sql_dstore.members


def test_members(sql_dstore):
    for m in sql_dstore:
        assert isinstance(m, DataMemberABC)
    assert len(sql_dstore) == len(sql_dstore.completed) + len(sql_dstore.not_completed)


def test_len(sql_dstore, ro_dir_dstore):
    """DataStore returns correct len"""
    expect = len(list(ro_dir_dstore.source.glob("*.fasta")))
    assert expect == len(sql_dstore) == len(sql_dstore.members)


def test_md5_sum(sql_dstore):
    for m in sql_dstore.members:
        data = m.read()
        md5 = sql_dstore.md5(m.unique_id)
        assert md5 == get_text_hexdigest(data)


def test_iterall(sql_dstore, ro_dir_dstore):
    expect = {fn.name for fn in ro_dir_dstore.source.glob("*.fasta")}
    got = {Path(m.unique_id).name for m in sql_dstore}
    assert expect == got


def test_read(full_dstore_sqlite):
    """correctly read content"""
    c = full_dstore_sqlite.completed[0]
    nc = full_dstore_sqlite.not_completed[0]
    log = full_dstore_sqlite.logs[0]
    for r in (c, nc, log):
        assert isinstance(r.read(), str)


def test_read_log(sql_dstore, log_data):
    """correctly read content"""
    sql_dstore.write_log(unique_id="brca1.fasta", data=log_data)

    got = sql_dstore.read(f"{_LOG_TABLE}/brca1.fasta")
    assert got == log_data


def test_write(sql_dstore, ro_dir_dstore):
    """correctly write content"""
    expect = Path(ro_dir_dstore.source / "brca1.fasta").read_text()
    identifier = "result_table/brca1.fasta"
    sql_dstore.write(unique_id=identifier, data=expect)
    got = sql_dstore.read(identifier)
    assert got == expect


def test_write_if_member_exists(sql_dstore, ro_dir_dstore):
    """correctly write content"""
    expect = Path(ro_dir_dstore.source / "brca1.fasta").read_text()
    identifier = "result_table/brca1.fasta"
    sql_dstore.write(unique_id=identifier, data=expect)
    got = sql_dstore.read(identifier)
    assert got == expect
    sql_dstore._mode = OVERWRITE
    identifier = "result_table/brca1.fasta"
    sql_dstore.write(unique_id=identifier, data=expect)
    got = sql_dstore.read(identifier)
    assert got == expect


def test_summary_logs(full_dstore_sqlite):
    # log summary has a row per log file and a column for each property
    got = full_dstore_sqlite.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


def test_summary_not_completed(full_dstore_sqlite):
    got = full_dstore_sqlite.summary_not_completed
    assert got.shape >= (1, 1)
    assert isinstance(got, Table)


def test_no_not_completed_subdir(full_dstore_sqlite):
    expect = f"{len(full_dstore_sqlite.completed)+len(full_dstore_sqlite.not_completed)}x member"
    assert repr(full_dstore_sqlite).startswith(expect)
    # first remove not_completed directory
    full_dstore_sqlite.drop_not_completed()
    # test repr work without not_completed directory
    expect = f"{len(full_dstore_sqlite.completed)}x member"
    assert repr(full_dstore_sqlite).startswith(expect)
    expect = f"{len(full_dstore_sqlite)}x member"
    assert repr(full_dstore_sqlite).startswith(expect)
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed)


def test_describe(full_dstore_sqlite):
    got = full_dstore_sqlite.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


def test_pickleable_roundtrip(full_dstore_sqlite):
    """pickling of data stores should be reversible"""
    re_dstore = loads(dumps(full_dstore_sqlite))
    assert str(re_dstore) == str(full_dstore_sqlite)
    assert re_dstore[0].read() == full_dstore_sqlite[0].read()


def test_pickleable_member_roundtrip(full_dstore_sqlite):
    """pickling of data store members should be reversible"""
    re_member = loads(dumps(full_dstore_sqlite[0]))
    data = re_member.read()
    assert len(data) > 0


def test_getitem(full_dstore_sqlite):
    with pytest.raises(IndexError):
        _ = full_dstore_sqlite[len(full_dstore_sqlite)]

    last = full_dstore_sqlite[-1]
    first = full_dstore_sqlite[0]
    assert last.unique_id != first.unique_id


def test_empty_data_store(db_dir):
    dstore = SqliteDataStore(db_dir, mode=OVERWRITE)
    assert 0 == len(dstore)


def test_no_logs(db_dir):
    dstore = SqliteDataStore(db_dir, mode=OVERWRITE)
    assert len(dstore.logs) == 0


def test_limit_datastore(full_dstore_sqlite):
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed) + len(
        full_dstore_sqlite.not_completed
    )
    full_dstore_sqlite._limit = len(full_dstore_sqlite.completed)
    full_dstore_sqlite.drop_not_completed()
    assert len(full_dstore_sqlite) == full_dstore_sqlite._limit
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed)


def test_validate(full_dstore_sqlite):
    r = full_dstore_sqlite.validate()
    assert r.shape == (3, 2)


def test_no_not_completed(sql_dstore):
    assert len(sql_dstore.not_completed) == 0


def test_write_read_only_datastore(ro_sql_dstore):
    with pytest.raises(IOError):
        ro_sql_dstore.write(unique_id="brca1.fasta", data="test data")


def test_write_read_only_member(ro_sql_dstore):
    with pytest.raises(IOError):
        ro_sql_dstore[0].write("test data")


def test_new_write_read(full_dstore_sqlite):
    """correctly write content"""
    identifier = "test1.fasta"
    data = "test data"
    m = full_dstore_sqlite.write(unique_id=identifier, data=data)
    got = full_dstore_sqlite.read(m.unique_id)
    assert got == data


'''
def test_append to reset log_id

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
'''


def test_summary_logs(full_dstore_sqlite):
    # log summary has a row per log file and a column for each property
    got = full_dstore_sqlite.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


def test_summary_not_completed(full_dstore_sqlite):
    got = full_dstore_sqlite.summary_not_completed
    assert got.shape >= (1, 1)
    assert isinstance(got, Table)

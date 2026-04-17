from datetime import UTC, datetime
from pathlib import Path

import pytest
from scinexus.composable import NotCompleted
from scinexus.data_store import (
    APPEND,
    OVERWRITE,
    READONLY,
    DataStoreDirectory,
)
from scinexus.io import open_data_store
from scinexus.sqlite_data_store import (
    LOG_TABLE,
    RESULT_TABLE,
    DataStoreSqlite,
    has_valid_schema,
)

from cogent3 import get_app
from cogent3.core.table import Table


@pytest.fixture
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("sqldb")


@pytest.fixture
def db_dir(tmp_dir):
    return tmp_dir / "data.sqlitedb"


@pytest.fixture
def ro_dir_dstore(DATA_DIR):
    return DataStoreDirectory(DATA_DIR, suffix="fasta")


@pytest.fixture
def completed_objects(ro_dir_dstore):
    return {f"{Path(m.unique_id).stem}": m.read() for m in ro_dir_dstore}


@pytest.fixture
def nc_objects():
    return {
        f"id_{i}": NotCompleted("ERROR", "location", "message", source=f"id_{i}")
        for i in range(3)
    }


@pytest.fixture
def log_data():
    path = Path(__file__).parent.parent / "data" / "scitrack.log"
    return path.read_text()


@pytest.fixture
def full_dstore_sqlite(db_dir, nc_objects, completed_objects, log_data):
    dstore = DataStoreSqlite(db_dir, mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())
    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)
    dstore.write_log(unique_id="scitrack.log", data=log_data)
    yield dstore
    dstore.close()


@pytest.fixture
def dstore_on_disk(full_dstore_sqlite):
    path = full_dstore_sqlite.source
    full_dstore_sqlite.close()
    return path


def test_invalid_sqlite(full_dstore_sqlite):
    query = "CREATE TABLE newtable (log_id INTEGER PRIMARY KEY)"
    full_dstore_sqlite.db.execute(query)
    assert not has_valid_schema(full_dstore_sqlite.db)


def test_open_existing(dstore_on_disk):
    ro = DataStoreSqlite(dstore_on_disk, mode=READONLY)
    assert len(ro) > 0
    assert len(ro.completed) > 0
    assert len(ro.not_completed) > 0
    assert len(ro.logs) > 0
    for attr in ("summary_logs", "summary_not_completed", "describe"):
        assert isinstance(getattr(ro, attr), Table)


def test_open_to_append(dstore_on_disk):
    DataStoreSqlite(dstore_on_disk, mode=APPEND)


def test_open_to_write(dstore_on_disk):
    DataStoreSqlite(dstore_on_disk, mode=OVERWRITE)


def test_db_creation():
    db = DataStoreSqlite(":memory:", mode=OVERWRITE)
    db = db.db
    result = db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    assert len(result) == 4
    created_names = {r["name"] for r in result}
    assert created_names == {
        LOG_TABLE,
        RESULT_TABLE,
        "state",
        "citations",
    }
    rows = db.execute(f"Select * from {LOG_TABLE}").fetchall()
    assert len(rows) == 0


def test_db_init_log():
    dstore = DataStoreSqlite(":memory:", mode=OVERWRITE)
    dstore._init_log()
    rows = dstore.db.execute(f"Select * from {LOG_TABLE}").fetchall()
    assert len(rows) == 1
    assert rows[0]["date"].date() == datetime.now(tz=UTC).date()


def test_not_completed(nc_objects):
    db = DataStoreSqlite(":memory:", mode=OVERWRITE)
    for unique_id, obj in nc_objects.items():
        db.write_not_completed(data=obj.to_json(), unique_id=unique_id)
    expect = len(nc_objects)
    query = f"SELECT count(*) as c FROM {RESULT_TABLE} WHERE is_completed=?"
    got = db.db.execute(query, (0,)).fetchone()["c"]
    assert got == expect, f"Failed for {RESULT_TABLE} number of rows"
    assert len(db.not_completed) == expect


def test_limit_datastore(full_dstore_sqlite):  # new
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed) + len(
        full_dstore_sqlite.not_completed,
    )
    full_dstore_sqlite._limit = len(full_dstore_sqlite.completed) // 2
    full_dstore_sqlite._completed = []
    full_dstore_sqlite._not_completed = []
    assert (
        len(full_dstore_sqlite.completed)
        == len(full_dstore_sqlite.not_completed)
        == full_dstore_sqlite._limit
    )
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed) + len(
        full_dstore_sqlite.not_completed,
    )
    full_dstore_sqlite.drop_not_completed()
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed)
    assert len(full_dstore_sqlite.not_completed) == 0
    full_dstore_sqlite._limit = len(full_dstore_sqlite.completed) // 2
    full_dstore_sqlite._completed = []
    full_dstore_sqlite._not_completed = []
    assert (
        len(full_dstore_sqlite)
        == len(full_dstore_sqlite.completed)
        == full_dstore_sqlite._limit
    )
    assert len(full_dstore_sqlite.not_completed) == 0


def test_summary_logs(full_dstore_sqlite):
    # log summary has a row per log file and a column for each property
    got = full_dstore_sqlite.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


def test_no_not_completed_subdir(full_dstore_sqlite):
    expect = f"{len(full_dstore_sqlite.completed) + len(full_dstore_sqlite.not_completed)}x member"
    assert str(full_dstore_sqlite).startswith(expect)
    # first remove not_completed directory
    full_dstore_sqlite.drop_not_completed()
    # test repr work without not_completed directory
    expect = f"{len(full_dstore_sqlite.completed)}x member"
    assert str(full_dstore_sqlite).startswith(expect)
    expect = f"{len(full_dstore_sqlite)}x member"
    assert str(full_dstore_sqlite).startswith(expect)
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed)


def test_describe(full_dstore_sqlite):
    got = full_dstore_sqlite.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


def test_limit_datastore(full_dstore_sqlite):
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed) + len(
        full_dstore_sqlite.not_completed,
    )
    full_dstore_sqlite._limit = len(full_dstore_sqlite.completed)
    full_dstore_sqlite.drop_not_completed()
    assert len(full_dstore_sqlite) == full_dstore_sqlite._limit
    assert len(full_dstore_sqlite) == len(full_dstore_sqlite.completed)


def test_validate(full_dstore_sqlite):
    r = full_dstore_sqlite.validate()
    assert r.shape == (4, 2)


def test_set_record_type(full_dstore_sqlite):
    from scinexus.misc import get_object_provenance

    from cogent3 import make_table

    assert full_dstore_sqlite.record_type is None
    t = make_table(data={"a": [0, 2]})
    full_dstore_sqlite.record_type = t
    assert full_dstore_sqlite.record_type == get_object_provenance(t)


def test_lock_firsttime(full_dstore_sqlite):
    full_dstore_sqlite.db.execute("DELETE FROM state WHERE state_id=1")
    full_dstore_sqlite.lock()
    assert full_dstore_sqlite.locked
    full_dstore_sqlite.unlock()
    assert not full_dstore_sqlite.locked


@pytest.fixture
def md5_none(full_dstore_sqlite):
    """create a data store with empty md5 fields"""
    full_dstore_sqlite.db.execute(
        "UPDATE results SET md5=? WHERE record_id LIKE '%'",
        (None,),
    )
    return full_dstore_sqlite


def test_validate_missing_md5(md5_none):
    t = md5_none.validate()
    assert t["md5_missing", "Value"] == 9
    for c in ("md5_correct", "md5_incorrect"):
        assert t[c, "Value"] == 0


def _make_appendable_dstore(path, suffix):
    return open_data_store(path, suffix=suffix, mode="a")


def _make_and_run_proc(out_path, suffix, members):
    out_dstore = _make_appendable_dstore(out_path, suffix)
    loader = get_app("load_unaligned", moltype="dna", format_name="fasta")
    mlength = get_app("sample.min_length", 400)

    if suffix:
        writer = get_app("write_seqs", out_dstore, format_name="fasta")
    else:
        writer = get_app("write_db", out_dstore)

    app = loader + mlength + writer
    return app.apply_to(members, cleanup=True)


@pytest.mark.parametrize(
    ("name", "suffix"),
    [("appended", "fa"), ("appended.sqlitedb", None)],
)
def test_append_makes_logs(tmp_dir, ro_dir_dstore, name, suffix):
    # do half the records in the first call
    num = len(ro_dir_dstore.completed) // 2
    # make a path for writeable dstore
    out_path = tmp_dir / name
    got1 = _make_and_run_proc(out_path, suffix, ro_dir_dstore[:num])
    assert len(got1.logs) == 1

    # creating a separate instance should result in a
    # new log file
    got2 = _make_and_run_proc(out_path, suffix, ro_dir_dstore[num:])
    assert len(got2.logs) == 2
    # should be a row for each log
    summary = got2.summary_logs
    assert summary.shape[0] == 2


def test_summary_not_completed(nc_objects):
    dstore = open_data_store(":memory:", mode="w")
    writer = get_app("write_db", dstore)
    for nc in nc_objects.values():
        writer(nc)

    # relying on the fact that all nc_objects have same origin
    # and message, so those columns can be readily interrogated
    summary = dstore.summary_not_completed
    vals = summary.to_list(columns=["origin", "message", "num"])
    assert len(vals) == 1
    assert vals[0] == ["location", "'message'", 3]


@pytest.fixture
def sample_citations():
    from citeable import Software

    c1 = Software(author=["A Author"], title="Tool One", year=2024)
    c1.app = "app_one"
    c2 = Software(author=["B Author"], title="Tool Two", year=2025)
    c2.app = "app_two"
    return (c1, c2)


def test_write_citations_overwrites_sqlite(sample_citations):
    from citeable import Software

    dstore = DataStoreSqlite(":memory:", mode=OVERWRITE)
    dstore.write_citations(data=sample_citations)
    c3 = Software(author=["C Author"], title="Tool Three", year=2026)
    c3.app = "app_three"
    dstore.write_citations(data=(c3,))
    loaded = dstore._load_citations()
    assert len(loaded) == 1
    assert loaded[0].title == "Tool Three"


def test_summary_citations_sqlite(sample_citations):
    dstore = DataStoreSqlite(":memory:", mode=OVERWRITE)
    dstore.write_citations(data=sample_citations)
    table = dstore.summary_citations
    assert isinstance(table, Table)
    assert table.shape[0] == 2
    assert "app" in table.header
    assert "citation" in table.header


def test_old_schema_without_citations_table(tmp_dir):
    """Opening an old-style database without a citations table should not fail."""
    import sqlite3

    db_path = Path(str(tmp_dir)) / "old.sqlitedb"
    # Create a database with the old schema (no citations table)
    db = sqlite3.connect(str(db_path))
    db.execute(
        "CREATE TABLE state(state_id INTEGER PRIMARY KEY, record_type TEXT, lock_pid INTEGER)",
    )
    db.execute(
        f"CREATE TABLE {LOG_TABLE}(log_id INTEGER PRIMARY KEY, log_name TEXT, date timestamp, data BLOB)",
    )
    db.execute(
        f"CREATE TABLE {RESULT_TABLE}(record_id TEXT PRIMARY KEY, log_id INTEGER, md5 BLOB, is_completed INTEGER, data BLOB)",
    )
    db.close()

    # Should open without error in read-only mode
    dstore = DataStoreSqlite(db_path, mode=READONLY)
    assert dstore._load_citations() == []

    # summary_citations should return empty table
    table = dstore.summary_citations
    assert table.shape[0] == 0
    dstore.close()


def test_old_schema_write_citations_creates_table(tmp_dir, sample_citations):
    """Writing citations to an old-style database should create the table."""
    import sqlite3

    db_path = Path(str(tmp_dir)) / "old_rw.sqlitedb"
    # Create a database with the old schema (no citations table)
    db = sqlite3.connect(str(db_path))
    db.execute(
        "CREATE TABLE state(state_id INTEGER PRIMARY KEY, record_type TEXT, lock_pid INTEGER)",
    )
    db.execute(
        f"CREATE TABLE {LOG_TABLE}(log_id INTEGER PRIMARY KEY, log_name TEXT, date timestamp, data BLOB)",
    )
    db.execute(
        f"CREATE TABLE {RESULT_TABLE}(record_id TEXT PRIMARY KEY, log_id INTEGER, md5 BLOB, is_completed INTEGER, data BLOB)",
    )
    db.close()

    # Open in append mode and write citations
    dstore = DataStoreSqlite(db_path, mode=APPEND)
    dstore.write_citations(data=sample_citations)
    loaded = dstore._load_citations()
    assert len(loaded) == 2
    dstore.close()

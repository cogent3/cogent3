from pathlib import Path
from pickle import dumps, loads

import pytest

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store_new import (
    _LOG_TABLE,
    _MD5_TABLE,
    _NOT_COMPLETED_TABLE,
    OVERWRITE,
    READONLY,
    SKIP,
    DataStoreDirectory,
)


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


@pytest.fixture(scope="session")
def write_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    write_dir = tmp_dir / "write"
    write_dir.mkdir(parents=True, exist_ok=True)
    return write_dir


@pytest.fixture(scope="function")
def nc_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    nc_dir = tmp_dir / "nc_dir"
    logs_dir = nc_dir / _LOG_TABLE
    not_dir = nc_dir / _NOT_COMPLETED_TABLE
    md5_dir = nc_dir / _MD5_TABLE
    if (logs_dir / "scitrack.log").exists():
        (logs_dir / "scitrack.log").unlink()
    for fn in list(not_dir.glob("*.json")):
        fn.unlink()
    for fn in list(md5_dir.glob("*.txt")):
        fn.unlink()
    for fn in list(nc_dir.glob("*.*")):
        fn.unlink()
    if Path(logs_dir).exists():
        Path(logs_dir).rmdir()
    if Path(not_dir).exists():
        Path(not_dir).rmdir()
    if Path(md5_dir).exists():
        Path(md5_dir).rmdir()
    if Path(nc_dir).exists():
        Path(nc_dir).rmdir()
    nc_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    not_dir.mkdir(exist_ok=True)
    (logs_dir / "scitrack.log").write_text((DATA_DIR / "scitrack.log").read_text())
    return nc_dir


@pytest.fixture(scope="session")
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix="fasta", if_dest_exists=READONLY)


@pytest.fixture(scope="session")
def w_dstore(write_dir):
    return DataStoreDirectory(write_dir, suffix="fasta", if_dest_exists=OVERWRITE)


@pytest.fixture(scope="function")
def nc_dstore(nc_dir):
    dstore = DataStoreDirectory(nc_dir, suffix="fasta", if_dest_exists=SKIP)
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
    for fn in filenames:
        identifier = fn.name
        dstore.write(identifier, fn.read_text())
    return dstore


def test_contains(ro_dstore):
    """correctly identify when a data store contains a member"""
    assert "brca1.fasta" in ro_dstore


def test_iter(ro_dstore):
    """DataStore objects allow iteration over members"""
    members = list(ro_dstore)
    assert members == ro_dstore.members


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


def test_write_read_only_datastore(ro_dstore):
    with pytest.raises(IOError):
        ro_dstore.write("brca1.fasta", "test data")


def test_write_read_only_member(ro_dstore):
    with pytest.raises(IOError):
        ro_dstore[0].write("test data")


def test_write(fasta_dir, w_dstore):
    """correctly write content"""
    expect = Path(fasta_dir / "brca1.fasta").read_text()
    identifier = "brca1.fasta"
    w_dstore.write(identifier, expect)
    got = w_dstore.read(identifier)
    assert got == expect


def test_multi_write(fasta_dir, w_dstore):
    """correctly write multiple files to data store"""
    expect_a = Path(fasta_dir / "brca1.fasta").read_text()
    expect_b = Path(fasta_dir / "primates_brca1.fasta").read_text()
    identifier_a = "brca2.fasta"
    identifier_b = "primates_brca2.fasta"
    w_dstore.write(identifier_a, expect_a)
    w_dstore.write(identifier_b, expect_b)
    got_a = w_dstore.read(identifier_a)
    got_b = w_dstore.read(identifier_b)
    # check that both bits of data match
    assert got_a == expect_a
    assert got_b == expect_b


def test_append(w_dstore):
    """correctly write content"""
    identifier = "test1.fasta"
    data = "test data"
    w_dstore.write(identifier, data)
    got = w_dstore.read(identifier)
    assert got == data


def test_notcompleted(nc_dstore):
    assert len(nc_dstore.not_completed) == 3
    nc_dstore.drop_not_completed()
    assert len(nc_dstore.not_completed) == 0


def test_no_not_completed_subdir(nc_dstore):
    # first remove not_completed directory
    nc_dstore.drop_not_completed()
    Path(nc_dstore.source / _NOT_COMPLETED_TABLE).rmdir()
    # test repr work without not_completed directory
    assert not Path(nc_dstore.source / _NOT_COMPLETED_TABLE).exists()
    expect = "9x member"
    assert repr(nc_dstore).startswith(expect)
    not_dir = nc_dstore.source / _NOT_COMPLETED_TABLE
    not_dir.mkdir(exist_ok=True)


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

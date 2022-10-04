from pathlib import Path
from pickle import dumps, loads

import pytest

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store_new import (
    _LOG_TABLE,
    _NOT_COMPLETED_TABLE,
    OVERWRITE,
    READONLY,
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


@pytest.fixture(scope="session")
def ic_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    ic_dir = tmp_dir / "ic_dir"
    ic_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = ic_dir / fn.name
        dest.write_text(fn.read_text())
    logs_dir = ic_dir / _LOG_TABLE
    nc_dir = ic_dir / _NOT_COMPLETED_TABLE
    logs_dir.mkdir(exist_ok=True)
    nc_dir.mkdir(exist_ok=True)
    (logs_dir / "scitrack.log").write_text((DATA_DIR / "scitrack.log").read_text())
    nc = NotCompleted("FAIL", "dummy", "dummy_message", source="dummy_source")
    (nc_dir / "nc.json").write_text(nc.to_json())
    return ic_dir


@pytest.fixture(scope="session")
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix=".fasta", if_dest_exists=READONLY)


@pytest.fixture(scope="session")
def w_dstore(write_dir):
    return DataStoreDirectory(write_dir, if_dest_exists=OVERWRITE)


@pytest.fixture(scope="session")
def ic_dstore(ic_dir):
    return DataStoreDirectory(ic_dir, suffix=".fasta")


def test_contains(ro_dstore):
    """correctly identify when a data store contains a member"""
    assert "brca1.fasta" in ro_dstore


def test_iter(ro_dstore):
    """DataStore objects allow iteration over members"""
    members = list(ro_dstore)
    assert members == ro_dstore.members


def test_len(fasta_dir, ro_dstore):
    """DataStore returns correct len"""
    expect = len(list(fasta_dir.glob("*.fasta")))
    assert expect == len(ro_dstore) == len(ro_dstore.members)


def test_getitem(ro_dstore):
    with pytest.raises(IndexError):
        _ = ro_dstore[len(ro_dstore)]

    last = ro_dstore[-1]
    first = ro_dstore[0]
    assert last.unique_id != first.unique_id


def test_iterall(fasta_dir, ro_dstore):
    expect = {fn.name for fn in fasta_dir.glob("*.fasta")}
    got = {m.unique_id for m in ro_dstore}
    assert expect == got


def test_read(fasta_dir, ro_dstore):
    """correctly read content"""
    expect = (fasta_dir / "brca1.fasta").read_text()
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


def test_logs(ic_dstore):
    assert len(ic_dstore.logs) == 1
    log = ic_dstore.logs[0].read()
    assert isinstance(log, str)


def test_not_completed(ic_dstore):
    assert len(ic_dstore.not_completed) == 1
    nc = ic_dstore.not_completed[0].read()
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
    identifier_a = "brca1.fasta"
    identifier_b = "primates_brca1.fasta"
    w_dstore.write(identifier_a, expect_a)
    w_dstore.write(identifier_b, expect_b)
    got_a = w_dstore.read(identifier_a)
    got_b = w_dstore.read(identifier_b)
    # check that both bits of data match
    assert got_a == expect_a
    assert got_b == expect_b


'''
def test_write_wout_suffix(w_dstore):
    """appends suffix expected to records"""
    with pytest.raises(ValueError):
        w_dstore.write("1", str(dict(a=24, b="some text")))

    w_dstore.write("1.fasta", str(dict(a=24, b="some text")))
    assert len(w_dstore)== 1

@skipIf(sys.platform.lower() != "darwin", "broken on linux")
def test_md5_write(write__dir, w_dstore):
    """tracks md5 sums of written data"""
    expect = Path( write_dir / "brca1.fasta").read_text()

    identifier = make_identifier(w_dstore, "brca1.fasta", absolute=True)
    abs_id = w_dstore.write(identifier, expect)
    md5 = "05a7302479c55c0b5890b50f617c5642"
    assert w_dstore.md5(abs_id) == md5
    assert w_dstore[0].md5 == md5

    # does not have md5 if not set
    identifier = make_identifier(w_dstore, "brca1.fasta", absolute=True)
    abs_id = w_dstore.write(identifier, expect)
    got = w_dstore.md5(abs_id, force=False)
    assert got is None
    # but if you set force=True, you get it
    md5 = "05a7302479c55c0b5890b50f617c5642"
    got = w_dstore.md5(abs_id, force=True)
    assert got == md5


def test_add_file():
    """correctly add an arbitrarily named file"""
    data = Path(f"data{os.sep}brca1.fasta").read_text()

    log_path = os.path.join(dirname, "some.log")
    with open(log_path, "w") as out:
        out.write("some text")

    path = os.path.join(dirname, basedir)
    dstore = DataStoreDirectory(path, suffix=".fa", create=True)
    _ = dstore.write("brca1.fa", data)
    dstore.add_file(log_path)
    assert "some.log" in dstore
    assert os.path.exists(log_path)

    log_path = os.path.join(dirname, "some.log")
    with open(log_path, "w") as out:
        out.write("some text")

    path = os.path.join(dirname, basedir)
    dstore = DataStoreDirectory(path, suffix=".fa", create=True)
    _ = dstore.write("brca1.fa", data)
    dstore.add_file(log_path, cleanup=True)
    assert "some.log" in dstore
    assert not os.path.exists(log_path)

def test_make_identifier():
    """correctly construct an identifier for a new member"""

    if dirname.startswith(f".{os.sep}"):
        dirname = dirname[2:]

    path = os.path.join(dirname, basedir)
    base_path = path.replace(".zip", "")
    dstore = DataStoreDirectory(path, suffix=".json", create=True)
    name = "brca1.fasta"
    got = make_identifier(dstore, name, absolute = True)
    expect = os.path.join(base_path, name.replace("fasta", "json"))
    assert got == expect

    # now using a DataMember
    member = DataMember(
        os.path.join(f"blah{os.sep}blah", f"2-{name}"), None
    )
    got = make_identifier(dstore, member, absolute=True)
    expect = os.path.join(base_path, member.name.replace("fasta", "json"))
    assert got == expect


def test_summary_logs():
    ...

def test_summary_not_completed():
   ...

'''

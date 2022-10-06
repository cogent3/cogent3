import shutil

from pathlib import Path
from pickle import dumps, loads

import pytest

from scitrack import get_text_hexdigest

from cogent3.app.composable import NotCompleted
from cogent3.app.data_store_new import (
    _LOG_TABLE,
    _MD5_TABLE,
    _NOT_COMPLETED_TABLE,
    OVERWRITE,
    RAISE,
    READONLY,
    SKIP,
    DataStoreDirectory,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


DATA_DIR = Path(__file__).parent.parent / "data"


def retain_nc_dstore(dstore):
    # retain nc_dstore to the previous state
    ncdir = dstore.source
    nc = [
        NotCompleted(
            "FAIL", f"dummy{i}", f"dummy_message{i}", source=f"dummy_source{i}"
        )
        for i in range(3)
    ]
    for i, item in enumerate(nc):
        dstore.write_not_completed(f"nc{i+1}", item.to_json())
    assert 3 == len(dstore.not_completed)
    assert 3 == len(list((ncdir / _MD5_TABLE).glob("*.txt")))


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
    return DataStoreDirectory(fasta_dir, suffix="fasta", if_dest_exists=READONLY)


@pytest.fixture(scope="function")
def w_dstore(write_dir):
    return DataStoreDirectory(write_dir, suffix="fasta", if_dest_exists=OVERWRITE)


@pytest.fixture(scope="function")
def nc_dstore(nc_dir):
    dstore = DataStoreDirectory(nc_dir, suffix="fasta", if_dest_exists=SKIP)
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


def test_limit_datastore(nc_dstore):
    assert len(nc_dstore) == len(nc_dstore.completed) + len(nc_dstore.not_completed)
    nc_dstore._limit = len(nc_dstore.completed)
    nc_dstore.drop_not_completed()
    assert len(nc_dstore) == nc_dstore._limit
    assert len(nc_dstore) == len(nc_dstore.completed)


def test_md5_sum(nc_dstore):
    for m in nc_dstore.members:
        data = m.read()
        md5 = nc_dstore.md5(m.unique_id)
        assert md5 == get_text_hexdigest(data)

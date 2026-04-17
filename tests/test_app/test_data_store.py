import pathlib
import shutil
from pathlib import Path

import pytest
from scinexus.composable import NotCompleted
from scinexus.data_store import (
    MD5_TABLE,
    OVERWRITE,
    READONLY,
    DataStoreDirectory,
    ReadOnlyDataStoreZipped,
)
from scinexus.io import open_data_store

from cogent3.app import io as io_app
from cogent3.app.data_store import convert_directory_datastore
from cogent3.core.table import Table


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return Path(tmp_path_factory.mktemp("datastore"))


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    # this set's the working directory for all tests in this module
    # as a tmp dir
    monkeypatch.chdir(tmp_dir)


@pytest.fixture
def fasta_dir(DATA_DIR, tmp_dir):
    tmp_dir = Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    fasta_dir = tmp_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = fasta_dir / fn.name
        dest.write_text(fn.read_text())
    return fasta_dir


@pytest.fixture
def write_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    write_dir = tmp_dir / "write"
    write_dir.mkdir(parents=True, exist_ok=True)
    yield write_dir
    shutil.rmtree(write_dir)


@pytest.fixture
def nc_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    nc_dir = tmp_dir / "nc_dir"
    nc_dir.mkdir(parents=True, exist_ok=True)
    yield nc_dir
    shutil.rmtree(nc_dir)


@pytest.fixture
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix="fasta", mode=READONLY)


@pytest.fixture
def nc_dstore(DATA_DIR, nc_dir):
    dstore = DataStoreDirectory(nc_dir, suffix="fasta", mode=OVERWRITE)
    # write one log file
    log_filename = "scitrack.log"
    dstore.write_log(unique_id=log_filename, data=(DATA_DIR / log_filename).read_text())
    # write three not_completed file
    nc = [
        NotCompleted(
            "FAIL",
            f"dummy{i}",
            f"dummy_message{i}",
            source=f"dummy_source{i}",
        )
        for i in range(3)
    ]
    for i, item in enumerate(nc):
        dstore.write_not_completed(unique_id=f"nc{i + 1}", data=item.to_json())
    assert len(dstore.not_completed) == 3
    assert len(list((nc_dir / MD5_TABLE).glob("*.txt"))) == len(dstore)
    filenames = DATA_DIR.glob("*.fasta")
    # write six fasta file
    for fn in filenames:
        identifier = fn.name
        dstore.write(unique_id=identifier, data=fn.read_text())
    return dstore


@pytest.fixture
def completed_objects(ro_dstore):
    return {f"{Path(m.unique_id).stem}": m.read() for m in ro_dstore}


@pytest.fixture
def nc_objects():
    return {
        f"id_{i}": NotCompleted("ERROR", "location", "message", source=f"id_{i}")
        for i in range(3)
    }


@pytest.fixture
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


@pytest.fixture
def full_dstore(write_dir, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())

    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)

    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


def test_convert_directory_datastore(Sample_oldDirectoryDataStore, write_dir):
    new_dstore = convert_directory_datastore(
        Sample_oldDirectoryDataStore,
        write_dir,
        ".fasta",
    )
    assert len(new_dstore) == 6


def test_summary_not_completed_func(nc_objects):
    dstore = open_data_store(":memory:", mode="w")
    writer = io_app.write_db(dstore)
    for nc in nc_objects.values():
        writer(nc)

    got = dstore.summary_not_completed
    assert isinstance(got, Table)
    assert got.shape[0] >= 1


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


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_describe(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    got = full_dstore.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_validate(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    validation_table = full_dstore.validate()
    assert isinstance(validation_table, Table)
    assert validation_table["md5_correct", "Value"] == len(full_dstore)
    assert validation_table["md5_incorrect", "Value"] == 0
    assert validation_table["has_log", "Value"] is True


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_summary_logs(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    # log summary has a row per log file and a column for each property
    got = full_dstore.summary_logs
    assert got.shape == (1, 6)
    assert isinstance(got, Table)


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_summary_not_completed(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    got = full_dstore.summary_not_completed
    assert got.shape >= (1, 1)
    assert isinstance(got, Table)

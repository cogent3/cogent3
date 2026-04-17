import pathlib
import shutil
from itertools import product
from pathlib import Path

import pytest
from scinexus.data_store import (
    _MD5_TABLE,
    OVERWRITE,
    READONLY,
    DataStoreDirectory,
    ReadOnlyDataStoreZipped,
    get_data_source,
)

from cogent3.app import io as io_app
from cogent3.app import sample as sample_app
from cogent3.app.composable import NotCompleted
from cogent3.app.data_store import convert_directory_datastore, summary_not_completeds
from cogent3.core.table import Table
from cogent3.util.union_dict import UnionDict


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
    assert len(list((nc_dir / _MD5_TABLE).glob("*.txt"))) == len(dstore)
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


def test_summary_logs_missing_field(nc_dstore):
    log_path = Path(nc_dstore.source) / nc_dstore.logs[0].unique_id
    data = [
        l for l in log_path.read_text().splitlines() if "composable function" not in l
    ]
    log_path.write_text("\n".join(data))
    assert isinstance(nc_dstore.summary_logs, Table)


@pytest.mark.parametrize("use_dser", [False, True])
def test_summary_not_completed_func(nc_objects, use_dser):
    dstore = io_app.open_data_store(":memory:", mode="w")
    writer = io_app.write_db(dstore)
    deser = io_app.load_db().deserialiser if use_dser else None
    for nc in nc_objects.values():
        writer(nc)

    got = summary_not_completeds(dstore, deserialise=deser)
    assert isinstance(got, Table)
    if use_dser:
        assert got.shape[0] >= 1
    else:
        assert got.shape[0] == 0


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


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_zipped_contains(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    assert "formattest.fasta" in full_dstore


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_members(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    completed = full_dstore.completed
    not_completed = full_dstore.not_completed
    assert full_dstore.members == completed + not_completed


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_describe(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    got = full_dstore.describe
    assert got.shape >= (3, 2)
    assert isinstance(got, Table)


@pytest.mark.parametrize("dstore", ["full_dstore", "zipped_full"])
def test_iter(dstore, request):
    full_dstore = request.getfixturevalue(dstore)
    members = list(full_dstore)
    assert members == full_dstore.members


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


@pytest.fixture
def app_dstore_in(tmp_path):
    from scinexus.io import open_data_store

    from cogent3 import get_app

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
    from scinexus.io import open_data_store

    from cogent3 import get_app, make_aligned_seqs

    out = open_data_store(base_path=tmp_path / "demo", suffix="fa.gz", mode="w")
    writer = get_app("write_seqs", data_store=out)
    seqs = make_aligned_seqs(
        {"s1": "CG--T", "s2": "CGTTT"},
        moltype="dna",
        info={"source": "test"},
    )
    got = writer(seqs)  # pylint: disable=not-callable
    assert got, got


def test_apply_to_not_completed(nc_dstore, tmp_path):
    loader = io_app.load_unaligned()
    num_seqs = sample_app.take_n_seqs(number=3, fixed_choice=False)
    out_dstore = io_app.open_data_store(tmp_path / "output", suffix="fa", mode="w")
    writer = io_app.write_seqs(data_store=out_dstore, format_name="fasta")
    app = loader + num_seqs + writer
    fini = app.apply_to(nc_dstore)
    assert 0 < len(fini.completed) <= len(nc_dstore.completed)


@pytest.fixture
def sample_citations():
    from citeable import Software

    c1 = Software(author=["A Author"], title="Tool One", year=2024)
    c1.app = "app_one"
    c2 = Software(author=["B Author"], title="Tool Two", year=2025)
    c2.app = "app_two"
    return (c1, c2)


def test_summary_citations_directory(write_dir, sample_citations):
    dstore = DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)
    dstore.write_citations(data=sample_citations)
    table = dstore.summary_citations
    assert isinstance(table, Table)
    assert table.shape[0] == 2
    assert "app" in table.header
    assert "citation" in table.header


def test_write_bib_tilde_path(write_dir, sample_citations, HOME_TMP_DIR):
    dstore = DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)
    dstore.write_citations(data=sample_citations)
    bib_path = f"~/{HOME_TMP_DIR.name}/refs.bib"
    dstore.write_bib(bib_path)
    expected = pathlib.Path(bib_path).expanduser()
    assert expected.exists()
    content = expected.read_text()
    assert "Tool One" in content
    assert "Tool Two" in content


def test_summary_citations_zipped(write_dir, sample_citations):
    dstore = DataStoreDirectory(write_dir, suffix="fasta", mode=OVERWRITE)
    dstore.write_citations(data=sample_citations)
    source = pathlib.Path(dstore.source)
    path = shutil.make_archive(
        base_name=source.name,
        format="zip",
        base_dir=source.name,
        root_dir=source.parent,
    )
    zipped = ReadOnlyDataStoreZipped(pathlib.Path(path), suffix="fasta")
    table = zipped.summary_citations
    assert isinstance(table, Table)
    assert table.shape[0] == 2


def test_write_citations_zipped_raises(zipped_basic):
    zipped = ReadOnlyDataStoreZipped(zipped_basic, suffix="fasta")
    with pytest.raises(TypeError, match="read only"):
        zipped.write_citations(data=(None,))


def test_old_directory_store_without_citations(fasta_dir):
    """Opening a directory store created before citations were added works."""
    # fasta_dir has .fasta files but no .citations file
    dstore = DataStoreDirectory(fasta_dir, suffix="fasta", mode=READONLY)
    assert dstore._load_citations() == []
    table = dstore.summary_citations
    assert isinstance(table, Table)
    assert table.shape[0] == 0


_types = tuple(product((dict, UnionDict), (str, Path)))


@pytest.mark.parametrize(("container_type", "source_stype"), _types)
def test_get_data_source_dict(container_type, source_stype):
    """handles case where input is dict (sub)class instance with top level source key"""
    value = source_stype("some/path.txt")
    data = container_type(source=value)
    got = get_data_source(data)
    assert got == "path.txt"


@pytest.mark.parametrize("klass", [str, Path])
def test_get_data_source_seqcoll(klass):
    """handles case where input is sequence collection object"""
    from cogent3 import make_unaligned_seqs

    value = klass("some/path.txt")
    obj = make_unaligned_seqs(
        {"seq1": "ACGG"},
        moltype="dna",
        info={"random_key": 1234},
        source=value,
    )
    got = get_data_source(obj)
    assert got == "path.txt"

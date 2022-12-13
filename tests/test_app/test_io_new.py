from pathlib import Path

import numpy
import pytest

from numpy.testing import assert_allclose

from cogent3 import DNA
from cogent3.app import io_new as io_app
from cogent3.app.composable import NotCompleted, source_proxy
from cogent3.app.data_store_new import DataMember, DataStoreDirectory, Mode
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.maths.util import safe_log
from cogent3.parse.sequence import PARSERS
from cogent3.util.deserialise import deserialise_object
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
    return Path(tmpdir_factory.mktemp("io"))


@pytest.fixture(scope="function")
def w_dir_dstore(tmp_dir):
    return DataStoreDirectory(tmp_dir, mode="w")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    # this set's the working directory for all tests in this module
    # as a tmp dir
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="function")
def fasta_dir(tmp_dir):
    tmp_dir = Path(tmp_dir)
    filenames = DATA_DIR.glob("*.fasta")
    fasta_dir = tmp_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = fasta_dir / fn.name
        dest.write_text(fn.read_text())
    return fasta_dir


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
    wrote = writer(seqs[0], datamember.unique_id)
    assert isinstance(wrote, DataMember)


def test_source_proxy_simple(fasta_dir):
    """correctly writes sequences out"""
    from cogent3.app.composable import define_app
    from cogent3.app.typing import IdentifierType

    @define_app
    def get_bytes(path: IdentifierType) -> bytes:
        path = Path(path)
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
    got = list(reader.as_completed([path]))
    assert isinstance(got[0], source_proxy)


@pytest.mark.parametrize("suffix", ("nex", "paml", "fasta"))
def test_load_aligned(suffix):
    """should handle nexus too"""
    nexus_paths = DataStoreDirectory(DATA_DIR, suffix=suffix, limit=2)
    loader = io_app.load_aligned(format=suffix)
    results = [loader(m) for m in nexus_paths]
    for result in results:
        assert isinstance(result, ArrayAlignment)


def test_load_unaligned():
    """load_unaligned returns degapped sequence collections"""
    fasta_paths = DataStoreDirectory(DATA_DIR, suffix=".fasta", limit=2)
    fasta_loader = io_app.load_unaligned(format="fasta")
    for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
        assert isinstance(seqs, SequenceCollection)
        assert "-" not in "".join(seqs.to_dict().values())
        assert seqs.info.source == fasta_paths[i].unique_id


@pytest.mark.parametrize(
    "loader", (io_app.load_aligned, io_app.load_unaligned, io_app.load_tabular)
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
        for row in table.tolist():
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


# todo test objects where there is no unique_id provided, or inferrable,
# should return NotCompletedError

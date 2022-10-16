from pathlib import Path

import pytest

from cogent3.app import io_new as io_app_new
from cogent3.app.composable_new import (
    _source_wrapped,
    source_proxy,
)
from cogent3.app.data_store_new import (
    SKIP,
    DataMember,
    DataStoreDirectory,
)
from cogent3.core.alignment import ArrayAlignment
from cogent3.parse.sequence import PARSERS


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


def test_write_seqs(fasta_dir, tmp_dir):
    """correctly writes sequences out"""
    datastore = DataStoreDirectory(fasta_dir, suffix="fasta")
    datamember = datastore[0]
    data = datamember.read().splitlines()
    data = dict(iter(PARSERS["fasta".lower()](data)))
    seqs = ArrayAlignment(data=data, moltype=None)
    seqs.info.source = datastore.source
    writer = io_app_new.WriteSeqs(tmp_dir / "write", if_dest_exists=SKIP)
    wrote = writer(seqs[0], datamember.unique_id)
    assert isinstance(wrote, DataMember)


def test_source_proxy_simple(fasta_dir):
    """correctly writes sequences out"""
    datastore = DataStoreDirectory(fasta_dir, suffix="fasta")
    datamember = datastore[0]
    reader = io_app_new.get_bytes()
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

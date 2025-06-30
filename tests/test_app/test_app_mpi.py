from pathlib import Path

import pytest

from cogent3 import get_app, open_data_store
from cogent3.util import parallel


@pytest.fixture
def tmp_dir(tmpdir_factory):
    return Path(tmpdir_factory.mktemp("mpirun"))


@pytest.mark.skipif(not parallel.USING_MPI, reason="Not using MPI")
def test_write_db(tmp_dir):
    """writing with overwrite in MPI should reset db"""
    dstore = open_data_store("data", suffix="fasta")
    members = [m for m in dstore if m.unique_id != "brca1.fasta"]
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    reader = get_app("load_unaligned")
    aligner = get_app("align_to_ref")
    writer = get_app("write_db", out_dstore)
    process = reader + aligner + writer

    process.apply_to(
        members,
        show_progress=False,
        parallel=True,
        par_kw={"use_mpi": True},
    )

    expect = [str(m) for m in process.data_store]

    # now get read only and check what's in there
    result = open_data_store(out_dstore.source)
    got = [str(m) for m in result]

    assert got == expect


def is_master(n):
    return parallel.is_master_process()


@pytest.mark.skipif(not parallel.USING_MPI, reason="Not using MPI")
def test_is_master_process():
    """
    is_master_process() should return False
    for all child processes
    """
    assert parallel.is_master_process()  # this should be master!
    index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # but workers should not
    master_processes = sum(
        bool(result) for result in parallel.imap(is_master, index, use_mpi=True)
    )
    assert master_processes == 0

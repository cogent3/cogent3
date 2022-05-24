from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main, skipUnless

from cogent3.app import align as align_app
from cogent3.app import io as io_app
from cogent3.app.io import write_db
from cogent3.util import parallel


__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Sheng Han Moses Koh"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class MPITests(TestCase):
    basedir = "data"

    @skipUnless(parallel.USING_MPI, reason="Not using MPI")
    def test_write_db(self):
        """writing with overwrite in MPI should reset db"""
        dstore = io_app.get_data_store("data", suffix="fasta")
        members = dstore.filtered(callback=lambda x: "brca1.fasta" not in x.split("/"))
        with TemporaryDirectory(dir=".") as dirname:
            path = Path(dirname) / "delme.tinydb"
            reader = io_app.load_unaligned()
            aligner = align_app.align_to_ref()
            writer = write_db(path, create=True, if_exists="overwrite")
            process = reader + aligner + writer

            r = process.apply_to(
                members,
                show_progress=False,
                parallel=True,
                par_kw=dict(use_mpi=True),
            )

            expect = [str(m) for m in process.data_store]
            process.data_store.close()

            # now get read only and check what's in there
            result = io_app.get_data_store(path)
            got = [str(m) for m in result]

            assert got == expect


if __name__ == "__main__":
    main()

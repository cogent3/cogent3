import os
import pathlib

from os.path import join
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3 import DNA, open_data_store
from cogent3.app import align as align_app
from cogent3.app import io as io_app
from cogent3.app.composable import NotCompleted
from cogent3.app.io import write_db
from cogent3.app.result import generic_result


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _get_generic_result(source):
    """creates a generic result with a DNA moltype as the single value"""
    gr = generic_result(source=source)
    gr["dna"] = DNA
    return gr


class TestIo(TestCase):
    basedir = "data"

    def tearDown(self) -> None:
        path = pathlib.Path("delme.tinydb")
        if path.exists():
            path.unlink()

    def test_define_data_store(self):
        """returns an iterable data store"""
        found = open_data_store(self.basedir, suffix=".fasta")
        self.assertTrue(len(found) > 1)
        found = open_data_store(self.basedir, suffix=".fasta", limit=2)
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(open_data_store(self.basedir, suffix=".fasta*"))
        self.assertTrue(len(found) > 2)

        # with a wild-card suffix
        found = list(open_data_store(self.basedir, suffix="*"))
        self.assertEqual(len(os.listdir(self.basedir)), len(found))

        # raises ValueError if suffix not provided or invalid
        with self.assertRaises(ValueError):
            _ = open_data_store(self.basedir)

        with self.assertRaises(ValueError):
            _ = open_data_store(self.basedir, 1)

    def test_load_db_failure_json_file(self):
        """informative load_db error message when given a json file path"""
        # todo this test has a trapped exception about being unable to delete
        # a file
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme")
            writer = write_db(outpath, create=True, if_exists="ignore")
            gr = _get_generic_result(join("blah", "delme.json"))
            got = writer(gr)
            writer.data_store.db.close()
            reader = io_app.load_db()
            outpath = join(dirname, "dummy.json")
            with open(outpath, mode="w") as outfile:
                outfile.write("\n\n")

            got = reader(outpath)
            self.assertIsInstance(got, NotCompleted)
            self.assertTrue("json" in got.message)

    def test_write_json_no_info(self):
        """correctly writes an object with out an info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme")
            gr = _get_generic_result(join("blah", "delme.json"))
            writer = io_app.write_json(outdir, create=True)
            _ = writer(gr)
            reader = io_app.load_json()
            got = reader(writer.data_store[0])
            got.deserialised_values()
            self.assertEqual(got["dna"], DNA)

    def test_restricted_usage_of_tinydb_suffix(self):
        """can only use tinydb in a load_db, write_db context"""
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme.tinydb")
            for writer_class in (
                io_app.write_seqs,
                io_app.write_json,
                io_app.write_tabular,
            ):
                with self.assertRaises(ValueError):
                    writer_class(outdir, create=True, if_exists="skip")
            # but OK for write_db
            w = io_app.write_db(outdir, create=True, if_exists="skip")
            w.data_store.close()

    def test_write_db_parallel(self):
        """writing with overwrite in parallel should reset db"""
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme.tinydb")
            dstore = open_data_store(self.basedir, suffix="fasta")
            members = dstore.filtered(
                callback=lambda x: "brca1.fasta" not in x.split("/")
            )
            reader = io_app.load_unaligned()
            aligner = align_app.align_to_ref()
            writer = write_db(outdir, create=True, if_exists="overwrite")
            process = reader + aligner + writer

            _ = process.apply_to(
                members, show_progress=False, parallel=True, cleanup=True
            )
            expect = [str(m) for m in process.data_store]
            process.data_store.close()

            # now get read only and check what's in there
            result = open_data_store(outdir)
            got = [str(m) for m in result]
            self.assertNotEqual(got, [])
            result.close()
            self.assertEqual(got, expect)


if __name__ == "__main__":
    main()

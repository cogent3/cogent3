import json
import os
import shutil
import zipfile

from os.path import basename, join
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from unittest.mock import Mock, patch

from cogent3 import DNA
from cogent3.app import io as io_app
from cogent3.app.composable import NotCompleted
from cogent3.app.data_store import WritableZippedDataStore
from cogent3.app.io import write_db
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestIo(TestCase):
    basedir = "data"

    def test_findall(self):
        """find all files recursively"""
        found = list(io_app.findall(self.basedir, suffix=".fasta"))
        self.assertTrue(len(found) > 1)
        found = list(io_app.findall(self.basedir, suffix=".fasta", limit=2))
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.findall(self.basedir, suffix=".fasta*"))
        self.assertTrue(len(found) > 2)

    def test_findall_zip(self):
        """find all files recursively in a zip archive"""
        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, "new")
            shutil.make_archive(zip_path, "zip", self.basedir)
            zip_path = zip_path + ".zip"  # because shutil adds the suffix
            found = list(io_app.findall(zip_path, suffix="fasta"))
            self.assertTrue(len(found) > 1)
            found = list(io_app.findall(zip_path, suffix="fasta", limit=2))
            self.assertTrue(len(found) == 2)

            # and with a suffix
            found = list(io_app.findall(zip_path, suffix=".fasta*"))
            self.assertTrue(len(found) > 2)

    def test_define_data_store(self):
        """returns an iterable data store"""
        found = io_app.get_data_store(self.basedir, suffix=".fasta")
        self.assertTrue(len(found) > 1)
        found = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.get_data_store(self.basedir, suffix=".fasta*"))
        self.assertTrue(len(found) > 2)

        # with a wild-card suffix
        found = list(io_app.get_data_store(self.basedir, suffix="*"))
        self.assertEqual(len(os.listdir(self.basedir)), len(found))

        # raises ValueError if suffix not provided or invalid
        with self.assertRaises(ValueError):
            _ = io_app.get_data_store(self.basedir)

        with self.assertRaises(ValueError):
            _ = io_app.get_data_store(self.basedir, 1)

    def test_load_aligned(self):
        """correctly loads aligned seqs"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                self.assertEqual(aln.info.source, paths[i])

        fasta_paths = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        fasta_loader = io_app.load_aligned(format="fasta")
        validate(fasta_paths, fasta_loader)

    def test_load_aligned_from_zip(self):
        """correctly loads aligned seqs from a zip archive"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                # paths is only the basename when workjing with zip archives
                # whereas the inpath will have full path of zip archive
                self.assertEqual(aln.info.source, paths[i])
                self.assertEqual(aln.info.source, paths[i])

        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, self.basedir.replace(".zip", ""))
            shutil.make_archive(
                base_name=zip_path, root_dir=".", format="zip", base_dir=self.basedir
            )
            zip_path = zip_path + ".zip"  # because shutil adds the suffix
            fasta_paths = list(io_app.findall(zip_path, suffix=".fasta", limit=2))
            fasta_loader = io_app.load_aligned(format="fasta")
            validate(fasta_paths, fasta_loader)

    def test_load_unaligned(self):
        """load_unaligned returns degapped sequence collections"""
        fasta_paths = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        fasta_loader = io_app.load_unaligned(format="fasta")
        for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
            self.assertIsInstance(seqs, SequenceCollection)
            self.assertTrue("-" not in "".join(seqs.todict().values()))
            self.assertEqual(seqs.info.source, fasta_paths[i])

        # should also handle case where it's given an alignment/sequence
        # collection
        got = fasta_loader(seqs)
        self.assertEqual(got, seqs)

    def test_write_seqs(self):
        """correctly writes sequences out"""
        fasta_paths = list(io_app.findall(self.basedir, suffix=".fasta", limit=2))
        fasta_loader = io_app.load_aligned(format="fasta")
        alns = list(map(fasta_loader, fasta_paths))
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_seqs(dirname, if_exists="ignore")
            wrote = list(map(writer, alns))
            written = list(io_app.findall(dirname, suffix="fasta"))
            for i, wrote in enumerate(written):
                self.assertEqual(alns[i].info.stored, join(dirname, wrote))

    def test_load_json(self):
        """correctly loads an object from json"""
        from cogent3.app.data_store import make_record_for_json

        data = make_record_for_json("delme", DNA, True)
        data = json.dumps(data)
        # straight directory
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme.json")
            with open(outpath, "w") as outfile:
                outfile.write(data)
            reader = io_app.load_json()
            got = reader(outpath)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

        # zipped directory
        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, "delme.zip")
            outpath = "delme/delme.json"
            with zipfile.ZipFile(zip_path, "a") as out:
                out.writestr(outpath, data)

            dstore = io_app.get_data_store(zip_path, suffix="json")
            member = dstore.get_member("delme.json")
            reader = io_app.load_json()
            got = reader(member)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

    def test_write_db_load_db(self):
        """correctly write/load from tinydb"""
        data = DNA.to_json()
        # straight directory
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme")
            writer = write_db(outpath, create=True, if_exists="ignore")
            mock = patch("data.source", autospec=True)
            mock.to_json = DNA.to_json
            mock.source = join("blah", "delme.json")
            got = writer(mock)
            writer.data_store.db.close()
            dstore = io_app.get_data_store(f"{outpath}.tinydb", suffix="json")
            reader = io_app.load_db()
            got = reader(dstore[0])
            dstore.close()
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

    def test_load_tabular(self):
        """correctly loads tabular data"""
        rows = [["1", "2"], ["3", "4"], ["5", "6.5"]]
        table = Table(["A", "B"], rows=rows)
        load_table = io_app.load_tabular(sep="\t", with_header=True)
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme.tsv")
            table.write(outpath)
            new = load_table(outpath)
            self.assertEqual(new.title, "")
            self.assertEqual(type(new[0, "B"]), float)
            self.assertEqual(type(new[0, "A"]), int)
            outpath = join(dirname, "delme2.tsv")
            with open(outpath, "w") as out:
                out.write("\t".join(table.header[:1]) + "\n")
                for row in table.tolist():
                    row = "\t".join(map(str, row))
                    out.write(row + "\n")
            result = load_table(outpath)
            self.assertIsInstance(result, NotCompleted)

        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme.zip")
            dstore = WritableZippedDataStore(outpath, suffix="tsv", create=True)
            dstore.write("sample1.tsv", table.tostring("tsv"))
            new = load_table(dstore[0])
            self.assertEqual(type(new[0, "B"]), float)
            self.assertEqual(type(new[0, "A"]), int)

    def test_write_json_with_info(self):
        """correctly writes an object with info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme")
            mock = Mock()
            mock.to_rich_dict = DNA.to_rich_dict
            mock.info.source = join("blah", "delme.json")
            writer = io_app.write_json(outdir, create=True)
            _ = writer(mock)
            reader = io_app.load_json()
            got = reader(join(outdir, "delme.json"))
            self.assertEqual(got, DNA)

        # now with a zipped archive
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme.zip")
            mock = Mock()
            mock.to_rich_dict = DNA.to_rich_dict
            mock.info.source = join("blah", "delme.json")
            writer = io_app.write_json(outdir, create=True)
            identifier = writer(mock)
            reader = io_app.load_json()
            got = reader(writer.data_store[0])
            self.assertEqual(got, DNA)
            expect = join(outdir.replace(".zip", ""), "delme.json")
            if expect.startswith("." + os.sep):
                expect = expect[2:]
            self.assertEqual(identifier, expect)

    def test_write_json_no_info(self):
        """correctly writes an object with out an info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme")
            mock = patch("data.source", autospec=True)
            mock.to_rich_dict = DNA.to_rich_dict
            mock.source = join("blah", "delme.json")
            writer = io_app.write_json(outdir, create=True)
            _ = writer(mock)
            reader = io_app.load_json()
            got = reader(writer.data_store[0])
            self.assertEqual(got, DNA)

        # now with a zipped archive
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme.zip")
            mock = patch("data.source", autospec=True)
            mock.to_rich_dict = DNA.to_rich_dict
            mock.source = join("blah", "delme.json")
            writer = io_app.write_json(outdir, create=True)
            identifier = writer(mock)
            reader = io_app.load_json()
            # checking loadable from a data store member too
            got = reader(writer.data_store[0])
            self.assertEqual(got, DNA)
            expect = join(outdir.replace(".zip", ""), "delme.json")
            if expect.startswith("." + os.sep):
                expect = expect[2:]
            self.assertEqual(identifier, expect)


if __name__ == "__main__":
    main()

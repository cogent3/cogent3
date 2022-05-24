import bz2
import gzip
import os
import pathlib
import tempfile
import zipfile

from os import remove, rmdir
from os.path import exists
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3.util.io import (
    _path_relative_to_zip_parent,
    atomic_write,
    get_format_suffixes,
    open_,
    path_exists,
    remove_files,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class AtomicWriteTests(TestCase):
    """testing the atomic_write class."""

    def test_does_not_write_if_exception(self):
        """file does not exist if an exception raised before closing"""
        # create temp file directory
        with tempfile.TemporaryDirectory(".") as dirname:
            dirname = pathlib.Path(dirname)
            test_filepath = dirname / "Atomic_write_test"
            with self.assertRaises(AssertionError):
                with atomic_write(test_filepath, mode="w") as f:
                    f.write("abc")
                    raise AssertionError
            self.assertFalse(test_filepath.exists())

    def test_writes_compressed_formats(self):
        """correctly writes / reads different compression formats"""
        fpath = pathlib.Path("data/sample.tsv")
        with open(fpath) as infile:
            expect = infile.read()

        with tempfile.TemporaryDirectory(".") as dirname:
            dirname = pathlib.Path(dirname)
            for suffix in ["gz", "bz2", "zip"]:
                outpath = dirname / f"{fpath.name}.{suffix}"
                with atomic_write(outpath, mode="wt") as f:
                    f.write(expect)

                with open_(outpath) as infile:
                    got = infile.read()

                self.assertEqual(got, expect, msg=f"write failed for {suffix}")

    def test_rename(self):
        """Renames file as expected"""
        # create temp file directory
        with tempfile.TemporaryDirectory(".") as dirname:
            # create temp filepath
            dirname = pathlib.Path(dirname)
            test_filepath = dirname / "Atomic_write_test"
            # touch the filepath so it exists
            f = open(test_filepath, "w").close()
            self.assertTrue(exists(test_filepath))
            # file should overwrite file if file already exists
            with atomic_write(test_filepath, mode="w") as f:
                f.write("abc")

    def test_atomic_write_noncontext(self):
        """atomic write works as more regular file object"""
        with TemporaryDirectory(dir=".") as dirname:
            path = pathlib.Path(dirname) / "foo.txt"
            zip_path = path.parent / f"{path.name}.zip"
            aw = atomic_write(path, in_zip=zip_path, mode="w")
            aw.write("some data")
            aw.close()
            with open_(zip_path) as ifile:
                got = ifile.read()
            self.assertEqual(got, "some data")

    def test_open_handles_bom(self):
        """handle files with a byte order mark"""
        with TemporaryDirectory(dir=".") as dirname:
            # create the different file types
            dirname = pathlib.Path(dirname)

            text = "some text"

            # plain text
            textfile = dirname / "sample.txt"
            textfile.write_text(text, encoding="utf-8-sig")

            # gzipped
            gzip_file = dirname / "sample.txt.gz"
            with gzip.open(gzip_file, "wt", encoding="utf-8-sig") as outfile:
                outfile.write(text)

            # bzipped
            bzip_file = dirname / "sample.txt.bz2"
            with bz2.open(bzip_file, "wt", encoding="utf-8-sig") as outfile:
                outfile.write(text)

            # zipped
            zip_file = dirname / "sample.zip"
            with zipfile.ZipFile(zip_file, "w") as outfile:
                outfile.write(textfile, "sample.txt")

            for path in (bzip_file, gzip_file, textfile, zip_file):
                with open_(path) as infile:
                    got = infile.read()
                    self.assertEqual(got, text, msg=f"failed reading {path}")

    def test_aw_zip_from_path(self):
        """supports inferring zip archive name from path"""
        with TemporaryDirectory(dir=".") as dirname:
            path = pathlib.Path(dirname) / "foo.txt"
            zip_path = path.parent / f"{path.name}.zip"
            aw = atomic_write(zip_path, in_zip=True, mode="w")
            aw.write("some data")
            aw.close()
            with open_(zip_path) as ifile:
                got = ifile.read()
                self.assertEqual(got, "some data")

            path = pathlib.Path(dirname) / "foo2.txt"
            zip_path = path.parent / f"{path.name}.zip"
            aw = atomic_write(path, in_zip=zip_path, mode="w")
            aw.write("some data")
            aw.close()
            with open_(zip_path) as ifile:
                got = ifile.read()
                self.assertEqual(got, "some data")

    def test_expanduser(self):
        """expands user correctly"""
        # create temp file directory
        home = pathlib.Path("~").expanduser()
        with tempfile.TemporaryDirectory(dir=home) as dirname:
            # create temp filepath
            dirname = pathlib.Path(dirname)
            test_filepath = dirname / "Atomic_write_test"
            test_filepath = str(test_filepath).replace(str(home), "~")
            with atomic_write(test_filepath, mode="w") as f:
                f.write("abc")

    def test_path_relative_to_zip_parent(self):
        """correctly generates member paths for a zip archive"""
        zip_path = pathlib.Path("some/path/to/a/data.zip")
        for member in ("data/member.txt", "member.txt", "a/b/c/member.txt"):
            got = _path_relative_to_zip_parent(zip_path, pathlib.Path(member))
            self.assertEqual(got.parts[0], "data")


class UtilsTests(TestCase):
    """Tests of individual functions in utils"""

    def setUp(self):
        """ """
        self.files_to_remove = []
        self.dirs_to_remove = []

    def tearDown(self):
        """ """
        list(map(remove, self.files_to_remove))
        list(map(rmdir, self.dirs_to_remove))

    def test_remove_files(self):
        """Remove files functions as expected"""
        # create list of temp file paths
        test_filepaths = [
            tempfile.NamedTemporaryFile(prefix="remove_files_test").name
            for i in range(5)
        ]

        # try to remove them with remove_files and verify that an IOError is
        # raises
        self.assertRaises(OSError, remove_files, test_filepaths)
        # now get no error when error_on_missing=False
        remove_files(test_filepaths, error_on_missing=False)

        # touch one of the filepaths so it exists
        open(test_filepaths[2], "w").close()
        # check that an error is raised on trying to remove the files...
        self.assertRaises(OSError, remove_files, test_filepaths)
        # ... but that the existing file was still removed
        self.assertFalse(exists(test_filepaths[2]))

        # touch one of the filepaths so it exists
        open(test_filepaths[2], "w").close()
        # no error is raised on trying to remove the files
        # (although 4 don't exist)...
        remove_files(test_filepaths, error_on_missing=False)
        # ... and the existing file was removed
        self.assertFalse(exists(test_filepaths[2]))

    def test_get_format_suffixes_returns_lower_case(self):
        """should always return lower case"""
        a, b = get_format_suffixes("suffixes.GZ")
        self.assertTrue(a == None and b == "gz")
        a, b = get_format_suffixes("suffixes.ABCD")
        self.assertTrue(a == "abcd" and b == None)
        a, b = get_format_suffixes("suffixes.ABCD.BZ2")
        self.assertTrue(a == "abcd" and b == "bz2")
        a, b = get_format_suffixes("suffixes.abcd.BZ2")
        self.assertTrue(a == "abcd" and b == "bz2")
        a, b = get_format_suffixes("suffixes.ABCD.bz2")
        self.assertTrue(a == "abcd" and b == "bz2")

    def test_get_format_suffixes(self):
        """correctly return suffixes for compressed etc.. formats"""
        a, b = get_format_suffixes("no_suffixes")
        self.assertTrue(a == b == None)
        a, b = get_format_suffixes("suffixes.gz")
        self.assertTrue(a == None and b == "gz")
        a, b = get_format_suffixes("suffixes.abcd")
        self.assertTrue(a == "abcd" and b == None)
        a, b = get_format_suffixes("suffixes.abcd.bz2")
        self.assertTrue(a == "abcd" and b == "bz2")
        a, b = get_format_suffixes("suffixes.zip")
        self.assertTrue(a == None and b == "zip")

    def test_get_format_suffixes_pathlib(self):
        """correctly return suffixes for compressed etc.. formats from pathlib"""
        Path = pathlib.Path
        a, b = get_format_suffixes(Path("no_suffixes"))
        self.assertTrue(a == b == None)
        a, b = get_format_suffixes(Path("suffixes.gz"))
        self.assertTrue(a == None and b == "gz")
        a, b = get_format_suffixes(Path("suffixes.abcd"))
        self.assertTrue(a == "abcd" and b == None)
        a, b = get_format_suffixes(Path("suffixes.abcd.bz2"))
        self.assertTrue(a == "abcd" and b == "bz2")
        a, b = get_format_suffixes(Path("suffixes.zip"))
        self.assertTrue(a == None and b == "zip")

    def test_path_exists(self):
        """robustly identifies whether an object is a valid path and exists"""
        self.assertFalse(path_exists({}))
        self.assertFalse(path_exists("not an existing path"))
        self.assertFalse(path_exists("(a,b,(c,d))"))
        self.assertFalse(path_exists("(a:0.1,b:0.1,(c:0.1,d:0.1):0.1)"))
        # works for a Path instance
        p = pathlib.Path(__file__)
        self.assertTrue(path_exists(p))
        # or string instance
        self.assertTrue(path_exists(__file__))

    def test_open_reads_zip(self):
        """correctly reads a zip compressed file"""
        with TemporaryDirectory(dir=".") as dirname:
            text_path = os.path.join(dirname, "foo.txt")
            with open(text_path, "w") as f:
                f.write("any str")

            zip_path = os.path.join(dirname, "foo.zip")
            with zipfile.ZipFile(zip_path, "w") as zip:
                zip.write(text_path)

            with open_(zip_path) as got:
                self.assertEqual(got.readline(), "any str")

    def test_open_writes_zip(self):
        """correctly writes a zip compressed file"""
        with TemporaryDirectory(dir=".") as dirname:
            zip_path = pathlib.Path(dirname) / "foo.txt.zip"

            with open_(zip_path, "w") as f:
                f.write("any str")

            with zipfile.ZipFile(zip_path, "r") as zip:
                name = zip.namelist()[0]
                got = zip.open(name).read()
                self.assertEqual(got, b"any str")

    def test_open_zip_multi(self):
        """zip with multiple records cannot be opened using open_"""
        with TemporaryDirectory(dir=".") as dirname:
            text_path1 = os.path.join(dirname, "foo.txt")
            with open(text_path1, "w") as f:
                f.write("any str")

            text_path2 = os.path.join(dirname, "bar.txt")
            with open(text_path2, "w") as f:
                f.write("any str")

            zip_path = os.path.join(dirname, "foo.zip")
            with zipfile.ZipFile(zip_path, "w") as zip:
                zip.write(text_path1)
                zip.write(text_path2)

            with self.assertRaises(ValueError):
                open_(zip_path)


if __name__ == "__main__":
    main()

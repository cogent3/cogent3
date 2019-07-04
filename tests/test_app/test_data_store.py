import os
import shutil
import sys

from tempfile import TemporaryDirectory
from unittest import TestCase, main, skipIf

from cogent3.app.data_store import (
    DataStoreMember,
    ReadOnlyDirectoryDataStore,
    ReadOnlyZippedDataStore,
    SingleReadDataStore,
    WritableDirectoryDataStore,
    WritableZippedDataStore,
)
from cogent3.parse.fasta import MinimalFastaParser


class DataStoreBaseTests:
    basedir = "data"
    ReadClass = None
    WriteClass = None

    def test_findall(self):
        """correctly identify all files with a suffix"""
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        num = len(dstore.members)
        self.assertEqual(num, 5)

        dstore = self.ReadClass(self.basedir, suffix=".fasta", limit=2)

        num = len(dstore.members)
        self.assertEqual(num, 2)

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        # member path based on uncompressed
        basedir = self.basedir.split(".")[0]
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        for member in dstore.members:
            self.assertTrue(member.startswith(basedir))
            self.assertNotEqual(str(member), member.name)

    def test_absolute_identifier(self):
        """correctly returns the absolute identifier"""
        # member path based on uncompressed
        basedir = self.basedir.split(".")[0]
        expect = {
            os.path.join(basedir, n) for n in ("brca1.fasta", "primates_brca1.fasta")
        }
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        got = {dstore.get_absolute_identifier(p) for p in expect}
        self.assertEqual(got, expect)

    def test_contains(self):
        """correctly identify when a data store contains a member"""
        # member path based on uncompressed
        basedir = self.basedir.split(".")[0]
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        self.assertTrue("brca1.fasta" in dstore)
        self.assertTrue(f"{basedir}/brca1.fasta" in dstore)

    def test_get_member(self):
        """returns a matching member"""
        basedir = self.basedir.split(".")[0]
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        member = dstore.get_member("brca1.fasta")
        self.assertNotEqual(member, None)

    def test_iter(self):
        """DataStore objects allow iteration over members"""
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        members = [m for m in dstore]
        self.assertEqual(members, dstore.members)

    def test_len(self):
        """DataStore returns correct len"""
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        self.assertEqual(len(dstore), len(dstore.members))

    def test_make_identifier(self):
        """correctly construct an identifier for a new member"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            base_path = path.replace(".zip", "")
            dstore = self.WriteClass(path, suffix=".json", create=True)
            name = "brca1.fasta"
            got = dstore.make_absolute_identifier(name)
            expect = os.path.join(base_path, name.replace("fasta", "json"))
            self.assertEqual(got, expect)

            # now using a DataStoreMember
            member = DataStoreMember(os.path.join("blah/blah", f"2-{name}"), None)
            got = dstore.make_absolute_identifier(member)
            expect = os.path.join(base_path, member.name.replace("fasta", "json"))
            self.assertEqual(got, expect)

    def test_read(self):
        """correctly read content"""
        with open("data/brca1.fasta") as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        basedir = self.basedir.replace(".zip", "")
        data = dstore.read(os.path.join(basedir, "brca1.fasta"))
        data = data.splitlines()
        got = {l: s for l, s in MinimalFastaParser(data)}
        self.assertEqual(got, expect)

    # todo not really bnroken, but something to do with line-feeds I
    #  suspect. This means scitrack needs a more platform robust approach...
    @skipIf(sys.platform.lower() != "darwin", "broken on linux")
    def test_md5_read(self):
        """tracks md5 checksums of read data"""
        dstore = self.ReadClass(self.basedir, suffix=".fasta")
        basedir = self.basedir.replace(".zip", "")
        identifier = os.path.join(basedir, "brca1.fasta")
        md5 = "05a7302479c55c0b5890b50f617c5642"
        self.assertEqual(dstore.md5(identifier, force=True), md5)
        # this property also directly available on the member
        member = dstore.get_member(identifier)
        self.assertEqual(member.md5, md5)

    def test_write(self):
        """correctly write content"""
        with open("data/brca1.fasta") as infile:
            expect = infile.read()

        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            identifier = dstore.make_absolute_identifier("brca1.fasta")
            abs_id = dstore.write(identifier, expect)
            got = dstore.read(abs_id)
            self.assertEqual(got, expect)

    @skipIf(sys.platform.lower() != "darwin", "broken on linux")
    def test_md5_write(self):
        """tracks md5 sums of written data"""
        with open("data/brca1.fasta") as infile:
            expect = infile.read()

        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            identifier = dstore.make_absolute_identifier("brca1.fasta")
            abs_id = dstore.write(identifier, expect)
            md5 = "05a7302479c55c0b5890b50f617c5642"
            self.assertEqual(dstore.md5(abs_id), md5)
            # this property also directly available on the member
            self.assertEqual(dstore[0].md5, md5)

        # does not have md5 if not set
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True, md5=False)
            identifier = dstore.make_absolute_identifier("brca1.fasta")
            abs_id = dstore.write(identifier, expect)
            got = dstore.md5(abs_id, force=False)
            self.assertEqual(got, None)
            # but if you set force=True, you get it
            md5 = "05a7302479c55c0b5890b50f617c5642"
            got = dstore.md5(abs_id, force=True)
            self.assertEqual(got, md5)

    def test_multi_write(self):
        """correctly write multiple files to data store"""
        with open("data/brca1.fasta") as infile:
            expect_a = infile.read()

        with open("data/primates_brca1.fasta") as infile:
            expect_b = infile.read()

        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            identifier_a = dstore.make_absolute_identifier("brca1.fasta")
            identifier_b = dstore.make_absolute_identifier("primates_brca1.fasta")
            abs_id_a = dstore.write(identifier_a, expect_a)
            abs_id_b = dstore.write(identifier_b, expect_b)
            got_a = dstore.read(abs_id_a)
            got_b = dstore.read(abs_id_b)
            # check that both bits of data match
            self.assertEqual(got_a, expect_a)
            self.assertEqual(got_b, expect_b)

    def test_filter(self):
        """filter method should return correctly matching members"""
        dstore = self.ReadClass(self.basedir, suffix="*")
        got = [m.name for m in dstore.filtered(callback=lambda x: "brca1" in str(x))]
        self.assertTrue(len(set(got)), 2)
        got = dstore.filtered(pattern="*brca1*")
        expect = [
            path
            for path in os.listdir(self.basedir.replace(".zip", ""))
            if "brca1" in path
        ]
        self.assertEqual(len(got), len(expect))

    def test_pickleable_roundtrip(self):
        """pickling of data stores should be reversible"""
        from pickle import dumps, loads

        dstore = self.ReadClass(self.basedir, suffix="*")
        re_dstore = loads(dumps(dstore))
        got = re_dstore[0].read()
        self.assertEqual(str(dstore), str(re_dstore))
        self.assertEqual(dstore[0].read(), re_dstore[0].read())

    def test_pickleable_member_roundtrip(self):
        """pickling of data store members should be reversible"""
        from pickle import dumps, loads

        dstore = self.ReadClass(self.basedir, suffix="*")
        re_member = loads(dumps(dstore[0]))
        data = re_member.read()
        self.assertTrue(len(data) > 0)

    def test_add_file(self):
        """correctly add an arbitrarily named file"""
        with open("data/brca1.fasta") as infile:
            data = infile.read()

        with TemporaryDirectory(dir=".") as dirname:
            log_path = os.path.join(dirname, "some.log")
            with open(log_path, "w") as out:
                out.write("some text")

            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            _ = dstore.write("brca1.fa", data)
            dstore.add_file(log_path)
            self.assertTrue("some.log" in dstore)
            self.assertTrue(os.path.exists(log_path))

        with TemporaryDirectory(dir=".") as dirname:
            log_path = os.path.join(dirname, "some.log")
            with open(log_path, "w") as out:
                out.write("some text")

            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            _ = dstore.write("brca1.fa", data)
            dstore.add_file(log_path, cleanup=True)
            self.assertTrue("some.log" in dstore)
            self.assertFalse(os.path.exists(log_path))


class DirectoryDataStoreTests(TestCase, DataStoreBaseTests):
    basedir = "data"
    ReadClass = ReadOnlyDirectoryDataStore
    WriteClass = WritableDirectoryDataStore


class ZippedDataStoreTests(TestCase, DataStoreBaseTests):
    basedir = "data.zip"
    ReadClass = ReadOnlyZippedDataStore
    WriteClass = WritableZippedDataStore

    def setUp(self):
        basedir = self.basedir.split(".")[0]
        shutil.make_archive(
            base_name=basedir, format="zip", base_dir=basedir, root_dir="."
        )

    def tearDown(self):
        os.remove(self.basedir)

    def test_write_no_parent(self):
        """zipped data store handles archive with no parent dir"""
        self.WriteClass("delme.zip", create=True, suffix="fa")

    def test_store_suffix(self):
        """data store adds file suffix if not provided"""
        source = self.basedir.split(".")[0]
        dstore = self.ReadClass(source, suffix="*")
        self.assertEqual(dstore.source, self.basedir)
        self.assertTrue(len(dstore) > 1)


class SingleReadStoreTests(TestCase):
    basedir = "data/brca1.fasta"
    Class = SingleReadDataStore

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        dstore = self.Class(self.basedir)
        self.assertEqual(dstore.members, ["data/brca1.fasta"])

    def test_absolute_identifier(self):
        """correctly returns the absolute identifier"""
        dstore = self.Class(self.basedir)
        got = dstore.get_absolute_identifier("brca1.fasta")
        self.assertEqual(got, self.basedir)

    def test_contains(self):
        """correctly identify when a data store contains a member"""
        dstore = self.Class(self.basedir)
        self.assertTrue("brca1.fasta" in dstore)
        self.assertTrue(self.basedir in dstore)

    def test_read(self):
        """correctly read content"""
        with open("data/brca1.fasta") as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.Class(self.basedir, suffix=".fasta")
        data = dstore.read(self.basedir)
        data = data.splitlines()
        got = {l: s for l, s in MinimalFastaParser(data)}
        self.assertEqual(got, expect)


if __name__ == "__main__":
    main()

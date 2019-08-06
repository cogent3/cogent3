import os
import shutil
import sys

from tempfile import TemporaryDirectory
from unittest import TestCase, main, skipIf

from cogent3.app.data_store import (
    DataStoreMember,
    ReadOnlyDirectoryDataStore,
    ReadOnlyTinyDbDataStore,
    ReadOnlyZippedDataStore,
    SingleReadDataStore,
    WritableDirectoryDataStore,
    WritableTinyDbDataStore,
    WritableZippedDataStore,
)
from cogent3.parse.fasta import MinimalFastaParser


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


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
        self.assertTrue(f"{basedir}{os.sep}brca1.fasta" in dstore)

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
            if dirname.startswith("." + os.sep):
                dirname = dirname[2:]

            path = os.path.join(dirname, self.basedir)
            base_path = path.replace(".zip", "")
            dstore = self.WriteClass(path, suffix=".json", create=True)
            name = "brca1.fasta"
            got = dstore.make_absolute_identifier(name)
            expect = os.path.join(base_path, name.replace("fasta", "json"))
            self.assertEqual(got, expect)

            # now using a DataStoreMember
            member = DataStoreMember(
                os.path.join("blah" + os.sep + "blah", f"2-{name}"), None
            )
            got = dstore.make_absolute_identifier(member)
            expect = os.path.join(base_path, member.name.replace("fasta", "json"))
            self.assertEqual(got, expect)
            dstore.close()

    def test_read(self):
        """correctly read content"""
        with open("data" + os.sep + "brca1.fasta") as infile:
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
        with open("data" + os.sep + "brca1.fasta") as infile:
            expect = infile.read()

        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix=".fa", create=True)
            identifier = dstore.make_absolute_identifier("brca1.fasta")
            abs_id = dstore.write(identifier, expect)
            got = dstore.read(abs_id)
            self.assertEqual(got, expect)
            dstore.close()

    @skipIf(sys.platform.lower() != "darwin", "broken on linux")
    def test_md5_write(self):
        """tracks md5 sums of written data"""
        with open("data" + os.sep + "brca1.fasta") as infile:
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
            dstore.close()

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
            dstore.close()

    def test_multi_write(self):
        """correctly write multiple files to data store"""
        with open("data" + os.sep + "brca1.fasta") as infile:
            expect_a = infile.read()

        with open("data" + os.sep + "primates_brca1.fasta") as infile:
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
            dstore.close()

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
        with open("data" + os.sep + "brca1.fasta") as infile:
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
            dstore.close()

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
            dstore.close()


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


class TinyDBDataStoreTests(TestCase):
    basedir = "data"
    ReadClass = ReadOnlyTinyDbDataStore
    WriteClass = WritableTinyDbDataStore

    def setUp(self):
        dstore = ReadOnlyDirectoryDataStore(self.basedir, suffix="fasta")
        data = {m.name: m.read() for m in dstore}
        self.data = data

    def test_len(self):
        """len tindydb data store correct"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            for id_, data in self.data.items():
                identifier = dstore.make_relative_identifier(id_)
                dstore.write(identifier, data)
            self.assertTrue(len(dstore), len(self.data))
            dstore.close()

    def test_tiny_contains(self):
        """contains operation works for tinydb data store"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            for id_, data in self.data.items():
                identifier = dstore.make_relative_identifier(id_)
                got = dstore.write(identifier, data)
                # just string value matching element should return True
                self.assertTrue(got in dstore)
                # but if database member parent is not the same as the
                # data store, it will be False
                got.parent = "abcd"
                self.assertTrue(got not in dstore)
            self.assertTrue("brca1.json" in dstore)
            dstore.close()

    def test_add_file(self):
        """adding file to tinydb should work"""
        with TemporaryDirectory(dir=".") as dirname:
            log_path = os.path.join(dirname, "some.log")
            with open(log_path, "w") as out:
                out.write("some text")

            keys = [k for k in self.data.keys()]
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            identifier = dstore.make_relative_identifier(keys[0])
            got = dstore.write(identifier, self.data[keys[0]])
            path = dstore.add_file(log_path, keep_suffix=True, cleanup=False)
            self.assertTrue("some.log" in dstore)
            dstore.close()

    def test_tiny_get_member(self):
        """get member works on TinyDbDataStore"""
        with TemporaryDirectory(dir=".") as dirname:
            keys = [k for k in self.data.keys()]
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            identifier = dstore.make_relative_identifier(keys[0])
            inserted = dstore.write(identifier, self.data[keys[0]])
            got = dstore.get_member(identifier)
            self.assertEqual(got.name, identifier)
            self.assertEqual(got.id, inserted.id)
            self.assertEqual(got.read(), self.data[keys[0]])
            dstore.close()

    def test_tinydb_iter(self):
        """tinydb iter works"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            for id_, data in self.data.items():
                identifier = dstore.make_relative_identifier(id_)
                dstore.write(identifier, data)

            members = [m for m in dstore]
            self.assertEqual(members, dstore.members)
            dstore.close()

    def test_tinydb_filter(self):
        """filtering members by name"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            for id_, data in self.data.items():
                identifier = dstore.make_relative_identifier(id_)
                dstore.write(identifier, data)
            matches = dstore.filtered("*brca1*")
            self.assertGreater(len(matches), 2)
            dstore.close()

    def test_pickleable_roundtrip(self):
        """pickling of data stores should be reversible"""
        from pickle import dumps, loads

        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="ignore")
            for id_, data in self.data.items():
                identifier = dstore.make_relative_identifier(id_)
                dstore.write(identifier, data)
            dstore.db.storage.flush()  # make sure written to disk
            re_dstore = loads(dumps(dstore))
            got = re_dstore[0].read()
            re_dstore.close()
            self.assertEqual(str(dstore), str(re_dstore))
            self.assertEqual(got, dstore[0].read())
            dstore.close()

    def test_tiny_write_incomplete(self):
        """write an incomplete result to tinydb"""
        from cogent3.app.composable import NotCompleted

        keys = list(self.data)
        incomplete = [
            keys.pop(0),
            NotCompleted("FAIL", "somefunc", "checking", source="testing.txt"),
        ]
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            id_ = dstore.make_relative_identifier(incomplete[0])
            dstore.write_incomplete(id_, incomplete[1])
            for k in keys:
                id_ = dstore.make_relative_identifier(k)
                dstore.write(id_, self.data[k])
            dstore.close()

            # all records are contained
            dstore = self.ReadClass(path)
            for k in self.data:
                id_ = f"{k.split('.')[0]}.json"
                self.assertTrue(id_ in dstore)

            # but len(dstore) reflects only members with completed==True
            len_store = len(dstore)
            # but len(dstore.db) reflects all members
            self.assertEqual(len(dstore.db), len_store + 1)
            # the incomplete property contains the incomplete ones
            got = dstore.incomplete[0].read()
            self.assertTrue("notcompleted" in got["type"].lower())
            dstore.close()

    def test_summary_methods(self):
        """produce a table"""
        from cogent3.app.composable import NotCompleted

        keys = list(self.data)
        incomplete = [
            keys.pop(0),
            NotCompleted("FAIL", "somefunc", "checking", source="testing.txt"),
        ]
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            id_ = dstore.make_relative_identifier(incomplete[0])
            dstore.write_incomplete(id_, incomplete[1])
            for k in keys:
                id_ = dstore.make_relative_identifier(k)
                dstore.write(id_, self.data[k])
            # now add a log file
            dstore.add_file("data" + os.sep + "scitrack.log", cleanup=False)
            got = dstore.describe
            # table has rows for completed, incomplete and log
            self.assertEqual(got.shape, (3, 2))
            # log summary has a row per log file and a column for each property
            got = dstore.summary_logs
            self.assertEqual(got.shape, (1, 6))
            # incomplete summmary has a row per type and 5 columns
            got = dstore.summary_incomplete
            self.assertEqual(got.shape, (1, 5))
            dstore.close()

    def test_dblock(self):
        """locking/unlocking of db"""
        from cogent3.app.data_store import _db_lockid
        from pathlib import Path

        keys = list(self.data)
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, if_exists="overwrite")
            for k in keys:
                id_ = dstore.make_relative_identifier(k)
                dstore.write(id_, self.data[k])
            self.assertTrue(dstore.locked)
            dstore.unlock(force=True)
            # now introduce an artificial lock, making sure lock flushed to disk
            dstore.db.insert({"identifier": "LOCK", "pid": 123})
            dstore.db.storage.flush()
            self.assertTrue(dstore.locked)
            self.assertEqual(_db_lockid(dstore.source), 123)
            # now calling _source_create_delete with overwrite should have no
            # effect
            dstore._source_create_delete("overwrite", False)
            path = Path(dstore.source)
            self.assertTrue(path.exists())
            # unlocking with wrong pid has no effect
            dstore.unlock()
            self.assertTrue(dstore.locked)
            # but we can force it
            dstore.unlock(force=True)
            self.assertFalse(dstore.locked)
            dstore.close()
            # and now a call to _source_create_delete will delete
            dstore._source_create_delete("overwrite", False)
            self.assertFalse(path.exists())
            dstore.close()


class SingleReadStoreTests(TestCase):
    basedir = f"data{os.sep}brca1.fasta"
    Class = SingleReadDataStore

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        dstore = self.Class(self.basedir)
        self.assertEqual(dstore.members, [f"data{os.sep}brca1.fasta"])

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
        with open("data" + os.sep + "brca1.fasta") as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.Class(self.basedir, suffix=".fasta")
        data = dstore.read(self.basedir)
        data = data.splitlines()
        got = {l: s for l, s in MinimalFastaParser(data)}
        self.assertEqual(got, expect)


if __name__ == "__main__":
    main()

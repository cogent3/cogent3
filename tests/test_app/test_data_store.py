import os
import shutil
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3.app.data_store import (SingleReadDataStore,
                                    ReadOnlyDirectoryDataStore,
                                    WritableDirectoryDataStore,
                                    WritableZippedDataStore,
                                    ReadOnlyZippedDataStore, DataStoreMember, )
from cogent3.parse.fasta import MinimalFastaParser


class DataStoreBaseTests:
    basedir = 'data'
    ReadClass = None
    WriteClass = None

    def test_findall(self):
        """correctly identify all files with a suffix"""
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        num = len(dstore.members)
        self.assertEqual(num, 3)

        dstore = self.ReadClass(self.basedir, suffix='.fasta', limit=2)

        num = len(dstore.members)
        self.assertEqual(num, 2)

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        # member path based on uncompressed
        basedir = self.basedir.split('.')[0]
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        for member in dstore.members:
            self.assertTrue(member.startswith(basedir))
            self.assertNotEqual(str(member), member.name)

    def test_absolute_identifier(self):
        """correctly returns the absolute identifier"""
        # member path based on uncompressed
        basedir = self.basedir.split('.')[0]
        expect = {os.path.join(basedir, n) for n in ('brca1.fasta',
                                                     'primates_brca1.fasta')}
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        got = {dstore.get_absolute_identifier(p) for p in expect}
        self.assertEqual(got, expect)

    def test_contains(self):
        """correctly identify when a data store contains a member"""
        # member path based on uncompressed
        basedir = self.basedir.split('.')[0]
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        self.assertTrue('brca1.fasta' in dstore)
        self.assertTrue(f'{basedir}/brca1.fasta' in dstore)

    def test_iter(self):
        """DataStore objects allow iteration over members"""
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        members = [m for m in dstore]
        self.assertEqual(members, dstore.members)

    def test_len(self):
        """DataStore returns correct len"""
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        self.assertEqual(len(dstore), len(dstore.members))

    def test_make_identifier(self):
        """correctly construct an identifier for a new member"""
        with TemporaryDirectory(dir='.') as dirname:
            path = os.path.join(dirname, self.basedir)
            base_path = path.replace('.zip', '')
            dstore = self.WriteClass(path, suffix='.json', create=True)
            name = 'brca1.fasta'
            got = dstore.make_absolute_identifier(name)
            expect = os.path.join(base_path, name.replace('fasta', 'json'))
            self.assertEqual(got, expect)

            # now using a DataStoreMember
            member = DataStoreMember(os.path.join('blah/blah', f'2-{name}'),
                                     None)
            got = dstore.make_absolute_identifier(member)
            expect = os.path.join(
                base_path, member.name.replace('fasta', 'json'))
            self.assertEqual(got, expect)

    def test_read(self):
        """correctly read content"""
        with open('data/brca1.fasta') as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        basedir = self.basedir.replace('.zip', '')
        data = dstore.read(os.path.join(basedir, 'brca1.fasta'))
        data = data.splitlines()
        got = {l: s for l, s in MinimalFastaParser(data)}
        self.assertEqual(got, expect)

    def test_write(self):
        """correctly write content"""
        with open('data/brca1.fasta') as infile:
            expect = infile.read()

        with TemporaryDirectory(dir='.') as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix='.fa', create=True)
            identifier = dstore.make_absolute_identifier('brca1.fasta')
            abs_id = dstore.write(identifier, expect)
            got = dstore.read(abs_id)
            self.assertEqual(got, expect)

    def test_multi_write(self):
        """correctly write multiple files to data store"""
        with open('data/brca1.fasta') as infile:
            expect_a = infile.read()

        with open('data/primates_brca1.fasta') as infile:
            expect_b = infile.read()

        with TemporaryDirectory(dir='.') as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix='.fa', create=True)
            identifier_a = dstore.make_absolute_identifier('brca1.fasta')
            identifier_b = dstore.make_absolute_identifier(
                'primates_brca1.fasta')
            abs_id_a = dstore.write(identifier_a, expect_a)
            abs_id_b = dstore.write(identifier_b, expect_b)
            got_a = dstore.read(abs_id_a)
            got_b = dstore.read(abs_id_b)
            # check that both bits of data match
            self.assertEqual(got_a, expect_a)
            self.assertEqual(got_b, expect_b)

    def test_filter(self):
        """filter method should return correctly matching members"""
        dstore = self.ReadClass(self.basedir, suffix='*')
        got = [m.name for m in
               dstore.filtered(callback=lambda x: 'brca1' in str(x))]
        self.assertTrue(len(set(got)), 2)
        got = dstore.filtered(pattern='*brca1*')
        expect = [path for path in os.listdir(self.basedir.replace('.zip', ''))
                  if 'brca1' in path]
        self.assertEqual(len(got), len(expect))

    def test_pickleable(self):
        """data store members should be pickleable"""
        from pickle import dumps
        dstore = self.ReadClass(self.basedir, suffix='*')
        r = dumps(dstore[0])
        dumps(dstore)


class DirectoryDataStoreTests(TestCase, DataStoreBaseTests):
    basedir = 'data'
    ReadClass = ReadOnlyDirectoryDataStore
    WriteClass = WritableDirectoryDataStore


class ZippedDataStoreTests(TestCase, DataStoreBaseTests):
    basedir = 'data.zip'
    ReadClass = ReadOnlyZippedDataStore
    WriteClass = WritableZippedDataStore

    def setUp(self):
        basedir = self.basedir.split('.')[0]
        shutil.make_archive(base_name=basedir,
                            format='zip',
                            base_dir=basedir,
                            root_dir='.')

    def tearDown(self):
        os.remove(self.basedir)

    def test_write_no_parent(self):
        """zipped data store handles archive with no parent dir"""
        self.WriteClass('delme.zip', create=True, suffix='fa')


class SingleReadStoreTests(TestCase):
    basedir = 'data/brca1.fasta'
    Class = SingleReadDataStore

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        dstore = self.Class(self.basedir)
        self.assertEqual(dstore.members, ['data/brca1.fasta'])

    def test_absolute_identifier(self):
        """correctly returns the absolute identifier"""
        dstore = self.Class(self.basedir)
        got = dstore.get_absolute_identifier('brca1.fasta')
        self.assertEqual(got, self.basedir)

    def test_contains(self):
        """correctly identify when a data store contains a member"""
        dstore = self.Class(self.basedir)
        self.assertTrue('brca1.fasta' in dstore)
        self.assertTrue(self.basedir in dstore)

    def test_read(self):
        """correctly read content"""
        with open('data/brca1.fasta') as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.Class(self.basedir, suffix='.fasta')
        data = dstore.read(self.basedir)
        data = data.splitlines()
        got = {l: s for l, s in MinimalFastaParser(data)}
        self.assertEqual(got, expect)


if __name__ == '__main__':
    main()

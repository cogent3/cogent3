import os
import shutil
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3.app.data_store import (SingleReadDataStore,
                                    ReadOnlyDirectoryDataStore,
                                    WritableDirectoryDataStore,
                                    WritableZippedDataStore,
                                    ReadOnlyZippedDataStore, )
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
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        got = set(dstore.members)
        self.assertTrue({'brca1.fasta', 'primates_brca1.fasta'} < got)

    def test_absolute_identifier(self):
        """correctly returns the absolute identifier"""
        expect = {os.path.join(self.basedir, n) for n in ('brca1.fasta',
                                                          'primates_brca1.fasta')}
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        got = {dstore.get_absolute_identifier(p) for p in expect}
        self.assertEqual(got, expect)

    def test_contains(self):
        """correctly identify when a data store contains a member"""
        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        self.assertTrue('brca1.fasta' in dstore)
        self.assertTrue(f'{self.basedir}/brca1.fasta' in dstore)

    def test_make_identifier(self):
        """correctly construct an identifier for a new member"""
        with TemporaryDirectory(dir='.') as dirname:
            path = os.path.join(dirname, self.basedir)
            dstore = self.WriteClass(path, suffix='.json', create=True)
            name = 'brca1.fasta'
            got = dstore.make_absolute_identifier(name)
            expect = os.path.join(path, name.replace('fasta', 'json'))
            self.assertEqual(got, expect)

    def test_read(self):
        """correctly read content"""
        with open('data/brca1.fasta') as infile:
            expect = {l: s for l, s in MinimalFastaParser(infile)}

        dstore = self.ReadClass(self.basedir, suffix='.fasta')
        data = dstore.read(os.path.join(self.basedir, 'brca1.fasta'))
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


class WriteableDirectoryDataStore(object):
    pass


class DirectoryDataStoreTests(TestCase, DataStoreBaseTests):
    ReadClass = ReadOnlyDirectoryDataStore
    WriteClass = WritableDirectoryDataStore


class ZippedDataStoreTests(TestCase, DataStoreBaseTests):
    basedir = 'data.zip'
    ReadClass = ReadOnlyZippedDataStore
    WriteClass = WritableZippedDataStore

    def setUp(self):
        basedir = self.basedir.split('.')[0]
        shutil.make_archive(basedir, 'zip', basedir)

    def tearDown(self):
        os.remove(self.basedir)


class SingleReadStoreTests(TestCase):
    basedir = 'data/brca1.fasta'
    Class = SingleReadDataStore

    def test_get_relative_identifier(self):
        """correctly returns the relative identifier"""
        dstore = self.Class(self.basedir)
        self.assertEqual(dstore.members, ['brca1.fasta'])

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

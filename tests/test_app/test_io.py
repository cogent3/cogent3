from os.path import join, basename
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
from unittest import TestCase, main
from cogent3.app import io as io_app
from cogent3 import DNA
from cogent3.core.alignment import SequenceCollection, ArrayAlignment

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestIo(TestCase):
    basedir = 'data'

    def test_findall(self):
        """find all files recursively"""
        found = list(io_app.findall(self.basedir, suffix='.fasta'))
        self.assertTrue(len(found) > 1)
        found = list(io_app.findall(self.basedir, suffix='.fasta', limit=2))
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.findall(self.basedir, suffix='.fasta*'))
        self.assertTrue(len(found) > 2)

    def test_findall_zip(self):
        """find all files recursively in a zip archive"""
        with TemporaryDirectory(dir='.') as dirname:
            zip_path = join(dirname, 'new')
            shutil.make_archive(zip_path, 'zip', self.basedir)
            zip_path = zip_path + '.zip'  # because shutil adds the suffix
            found = list(io_app.findall(zip_path, suffix='fasta'))
            self.assertTrue(len(found) > 1)
            found = list(io_app.findall(zip_path, suffix='fasta', limit=2))
            self.assertTrue(len(found) == 2)

            # and with a suffix
            found = list(io_app.findall(zip_path, suffix='.fasta*'))
            self.assertTrue(len(found) > 2)

    def test_load_aligned(self):
        """correctly loads aligned seqs"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                expect_source = join(loader.data_store.source, paths[i])
                self.assertEqual(aln.info.source, expect_source)

        fasta_paths = list(io_app.findall(self.basedir, suffix='.fasta',
                                          limit=2))
        fasta_loader = io_app.load_aligned(self.basedir, format='fasta')
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
                expect_source = join(loader.data_store.source, paths[i])
                self.assertEqual(basename(aln.info.source), paths[i])
                self.assertEqual(aln.info.source,
                                 expect_source)

        with TemporaryDirectory(dir='.') as dirname:
            zip_path = join(dirname, 'new')
            shutil.make_archive(zip_path, 'zip', self.basedir)
            zip_path = zip_path + '.zip'  # because shutil adds the suffix
            fasta_paths = list(io_app.findall(zip_path, suffix='.fasta',
                                              limit=2))
            fasta_loader = io_app.load_aligned(format='fasta',
                                               data_path=zip_path)
            validate(fasta_paths, fasta_loader)

    def test_load_unaligned(self):
        """load_unaligned returns degapped sequence collections"""
        fasta_paths = list(io_app.findall(
            self.basedir, suffix='.fasta', limit=2))
        fasta_loader = io_app.load_unaligned(self.basedir, format='fasta')
        for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
            self.assertIsInstance(seqs, SequenceCollection)
            self.assertTrue('-' not in ''.join(seqs.todict().values()))
            expect_source = join(self.basedir, fasta_paths[i])
            self.assertEqual(seqs.info.source, expect_source)

    def test_write_seqs(self):
        """correctly writes sequences out"""
        fasta_paths = list(io_app.findall(self.basedir, suffix='.fasta',
                                          limit=2))
        fasta_loader = io_app.load_aligned(self.basedir, format='fasta',
                                           suffix='fasta')
        alns = list(map(fasta_loader, fasta_paths))
        with TemporaryDirectory(dir='.') as dirname:
            writer = io_app.write_seqs(dirname, if_exists='ignore')
            wrote = list(map(writer, alns))
            written = list(io_app.findall(dirname, suffix='fasta'))
            for i, wrote in enumerate(written):
                self.assertEqual(alns[i].info.stored, join(dirname, wrote))


if __name__ == '__main__':
    main()

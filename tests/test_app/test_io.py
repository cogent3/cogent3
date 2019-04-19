import zipfile
from os.path import join, basename
from tempfile import TemporaryDirectory
import shutil
from unittest import TestCase, main
from unittest.mock import Mock, patch
from cogent3 import DNA
from cogent3.app import io as io_app
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

    def test_define_data_store(self):
        """returns an iterable data store"""
        found = io_app.get_data_store(self.basedir, suffix='.fasta')
        self.assertTrue(len(found) > 1)
        found = io_app.get_data_store(
            self.basedir, suffix='.fasta', limit=2)
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.get_data_store(
            self.basedir, suffix='.fasta*'))
        self.assertTrue(len(found) > 2)

    def test_load_aligned(self):
        """correctly loads aligned seqs"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                self.assertEqual(aln.info.source, paths[i])

        fasta_paths = io_app.get_data_store(self.basedir, suffix='.fasta',
                                            limit=2)
        fasta_loader = io_app.load_aligned(format='fasta')
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

        with TemporaryDirectory(dir='.') as dirname:
            zip_path = join(dirname, self.basedir.replace('.zip', ''))
            shutil.make_archive(base_name=zip_path,
                                root_dir='.',
                                format='zip',
                                base_dir=self.basedir)
            zip_path = zip_path + '.zip'  # because shutil adds the suffix
            fasta_paths = list(io_app.findall(zip_path, suffix='.fasta',
                                              limit=2))
            fasta_loader = io_app.load_aligned(format='fasta')
            validate(fasta_paths, fasta_loader)

    def test_load_unaligned(self):
        """load_unaligned returns degapped sequence collections"""
        fasta_paths = io_app.get_data_store(self.basedir, suffix='.fasta',
                                            limit=2)
        fasta_loader = io_app.load_unaligned(format='fasta')
        for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
            self.assertIsInstance(seqs, SequenceCollection)
            self.assertTrue('-' not in ''.join(seqs.todict().values()))
            self.assertEqual(seqs.info.source, fasta_paths[i])

        # should also handle case where it's given an alignment/sequence
        # collection
        got = fasta_loader(seqs)
        self.assertEqual(got, seqs)

    def test_write_seqs(self):
        """correctly writes sequences out"""
        fasta_paths = list(io_app.findall(self.basedir, suffix='.fasta',
                                          limit=2))
        fasta_loader = io_app.load_aligned(format='fasta', suffix='fasta')
        alns = list(map(fasta_loader, fasta_paths))
        with TemporaryDirectory(dir='.') as dirname:
            writer = io_app.write_seqs(dirname, if_exists='ignore')
            wrote = list(map(writer, alns))
            written = list(io_app.findall(dirname, suffix='fasta'))
            for i, wrote in enumerate(written):
                self.assertEqual(alns[i].info.stored, join(dirname, wrote))

    def test_load_json(self):
        """correctly loads an object from json"""
        data = DNA.to_json()
        # straight directory
        with TemporaryDirectory(dir='.') as dirname:
            outpath = join(dirname, 'delme.json')
            with open(outpath, 'w') as outfile:
                outfile.write(data)
            reader = io_app.load_json(dirname)
            got = reader(outpath)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

        # zipped directory
        with TemporaryDirectory(dir='.') as dirname:
            zip_path = join(dirname, 'delme.zip')
            outpath = 'delme/delme.json'
            with zipfile.ZipFile(zip_path, 'a') as out:
                out.writestr(outpath, data)

            reader = io_app.load_json(zip_path)
            got = reader(outpath)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

    def test_write_json_with_info(self):
        """correctly writes an object with info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir='.') as dirname:
            outdir = join(dirname, 'delme')
            mock = Mock()
            mock.to_json = DNA.to_json
            mock.info.source = join('blah', 'delme.json')
            writer = io_app.write_json(outdir, create=True)
            _ = writer(mock)
            reader = io_app.load_json(outdir)
            got = reader(join(outdir, 'delme.json'))
            self.assertEqual(got, DNA)

        # now with a zipped archive
        with TemporaryDirectory(dir='.') as dirname:
            outdir = join(dirname, 'delme.zip')
            mock = Mock()
            mock.to_json = DNA.to_json
            mock.info.source = join('blah', 'delme.json')
            writer = io_app.write_json(outdir, create=True)
            identifier = writer(mock)
            reader = io_app.load_json(outdir)
            got = reader(join(outdir.replace('.zip', ''), 'delme.json'))
            self.assertEqual(got, DNA)
            self.assertEqual(identifier, join(outdir.replace('.zip', ''),
                                              'delme.json'))

    def test_write_json_no_info(self):
        """correctly writes an object with out an info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir='.') as dirname:
            outdir = join(dirname, 'delme')
            mock = patch('data.source', autospec=True)
            mock.to_json = DNA.to_json
            mock.source = join('blah', 'delme.json')
            writer = io_app.write_json(outdir, create=True)
            _ = writer(mock)
            reader = io_app.load_json(outdir)
            got = reader(join(outdir, 'delme.json'))
            self.assertEqual(got, DNA)

        # now with a zipped archive
        with TemporaryDirectory(dir='.') as dirname:
            outdir = join(dirname, 'delme.zip')
            mock = patch('data.source', autospec=True)
            mock.to_json = DNA.to_json
            mock.source = join('blah', 'delme.json')
            writer = io_app.write_json(outdir, create=True)
            identifier = writer(mock)
            reader = io_app.load_json(outdir)
            got = reader(join(outdir.replace('.zip', ''), 'delme.json'))
            self.assertEqual(got, DNA)
            self.assertEqual(identifier, join(outdir.replace('.zip', ''),
                                              'delme.json'))


if __name__ == '__main__':
    main()

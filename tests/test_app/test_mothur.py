#!/usr/bin/env python
# test_mothur.py


from __future__ import with_statement
from cStringIO import StringIO
from os import remove, rmdir
from tempfile import mkdtemp, mkstemp
from cogent.util.unit_test import TestCase, main
from cogent.app.mothur import Mothur, mothur_from_file


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class MothurTests(TestCase):
    def setUp(self):
        self.small_fasta = (
            '>aaaaaa\nTAGGCTCTGATATAATAGCTCTC---------\n'
            '>cccccc\n------------TGACTACGCAT---------\n'
            '>bbbbbb\n----TATCGCTTCGACGATTCTCTGATAGAGA\n'
            )
        self.small_otus = (
            'unique\t3\tcccccc\tbbbbbb\taaaaaa\t\n'
            '0.62\t2\tbbbbbb,cccccc\taaaaaa\t\n'
            '0.67\t1\taaaaaa,bbbbbb,cccccc\t\n'
            )
        self.small_otus_parsed = [
            (float('0'), [['cccccc'], ['bbbbbb'], ['aaaaaa']]),
            (float('0.62'), [['bbbbbb', 'cccccc'], ['aaaaaa']]),
            (float('0.67'), [['aaaaaa', 'bbbbbb', 'cccccc']]),
            ]
        self.complement_fasta = (
            '>a\n--AGGGGTAATAA--\n'
            '>b\n--TTATTACCCCT--\n'
            '>c\n-------AAAAAA--\n'
            )
        self.complement_otus = (
            'unique\t3\tc\ta\tb\t\n'
            '0.50\t2\ta,c\tb\t\n'
            '1.00\t1\tb,a,c\t\n'
            )

    def test_get_help(self):
        """Mothur.getHelp() should return help string"""
        expected_help = (
            'See manual, available on the MOTHUR wiki:\n'
            'http://schloss.micro.umass.edu/mothur/'
            )
        self.assertEqual(Mothur.getHelp(), expected_help)

    def test_compile_mothur_script(self):
        """Mothur._compile_mothur_script() should return valid Mothur script"""
        app = Mothur()
        app._input_filename = 'test.fasta'
        observed_script = app._compile_mothur_script()
        expected_script = (
            '#unique.seqs(fasta=test.fasta); '
            'dist.seqs(fasta=test.unique.fasta); '
            'read.dist(column=test.unique.dist, name=test.names); '
            'cluster(method=furthest)')
        self.assertEqual(observed_script, expected_script)

    def test_get_result_paths(self):
        """Mothur._get_result_paths() should guess correct output paths"""
        # Implicitly tests all of the _derive_<something>_path() methods.
        app = Mothur()
        app._input_filename = 'test.fasta'
        paths = app._get_result_paths()
        observed_paths = dict([(k, v.Path) for (k, v) in paths.items()])
        expected_paths = {
            'log': 'mothur.logFile',
            'otu list': 'test.unique.fn.list',
            'distance matrix': 'test.unique.dist',
            'unique seqs': 'test.unique.fasta',
            'rank abundance': 'test.unique.fn.rabund',
            'unique names': 'test.names',
            'species abundance': 'test.unique.fn.sabund',
            }
        self.assertEqual(observed_paths, expected_paths)

    def test_working_directory(self):
        """Mothur.WorkingDir attribute should not be cast to FilePath object"""
        app = Mothur(WorkingDir='/tmp')
        self.assertEquals(str(app.WorkingDir), '/tmp')

    def test_call_with_multiline_string(self):
        """Mothur.__call__() should return correct otu's for input as single string"""
        app = Mothur()
        result = app(self.small_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()

    def test_call_with_lines(self):
        """Mothur.__call__() should return correct otu's for input as lines"""
        lines = self.small_fasta.split('\n')
        app = Mothur(InputHandler='_input_as_lines')
        result = app(lines)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()

    def test_call_with_path(self):
        """Mothur.__call__() should return correct otu's for input as path"""
        working_dir = mkdtemp()
        _, filename = mkstemp(dir=working_dir, suffix='.fasta')
        with open(filename, 'w') as f:
            f.write(self.small_fasta)
        app = Mothur(InputHandler='_input_as_path', WorkingDir=working_dir)
        result = app(filename)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        remove(filename)
        result.cleanUp()
        rmdir(working_dir)

    def test_call_with_working_dir(self):
        """Mothur.__call__() should return correct otu's when input dir is changed"""
        working_dir = mkdtemp()
        app = Mothur(WorkingDir=working_dir)
        result = app(self.small_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()
        rmdir(working_dir)

    def test_call_with_complement(self):
        """Mothur.__call__() should return correct otu's for input sequences which are reverse complements"""
        app = Mothur()
        result = app(self.complement_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.complement_otus)
        result.cleanUp()

    def test_mothur_from_file(self):
        """mothur_from_file() should return parsed otus"""
        f = StringIO(self.small_fasta)
        f.seek(0)
        parsed_otus = mothur_from_file(f)
        self.assertEquals(parsed_otus, self.small_otus_parsed)

if __name__ == '__main__':
    main()

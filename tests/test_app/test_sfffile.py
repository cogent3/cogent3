#!/usr/bin/env python
# test_sfffile.py

import os
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.parse.binary_sff import parse_binary_sff
from cogent.app.util import ApplicationError
from cogent.app.sfffile import Sfffile



__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.5.2-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class SfffileTests(TestCase):
    """Test the Sfffile application controller.
    """
    def setUp(self):
        test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.sff_fp = os.path.join(test_dir, 'data', 'test.sff')

    def _check_unmodified_sff_contents(self, sff_file):
        """Extracting repeated code from sfffile tests"""
        sff_file.seek(0)
        header, reads_gen = parse_binary_sff(sff_file)
        reads = list(reads_gen)

        self.assertEqual(header["number_of_reads"], 1)
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0]['Name'], 'FA6P1OK01CGMHQ')

    def test_exit_status(self):
        """Sfffile should raise ApplicationError if exit status is nonzero."""
        a = Sfffile()
        self.assertRaises(ApplicationError, a, 'an_sff_file_that_does_not_exist.sff')

    def test_call(self):
        """Simple sfffile call should produce expected output."""
        a = Sfffile()
        app_results = a(self.sff_fp)
        self._check_unmodified_sff_contents(app_results['sff'])
        app_results.cleanUp()

    def test_call_with_output_path(self):
        """Sfffile should store output to specified filepath."""
        _, output_fp = tempfile.mkstemp()
        a = Sfffile()
        a.Parameters['-o'].on(output_fp)
        app_results = a(self.sff_fp)

        self._check_unmodified_sff_contents(open(output_fp))
        self._check_unmodified_sff_contents(app_results['sff'])
        app_results.cleanUp()

    def test_call_with_included_accession_numbers(self):
        """Sfffile should include specified accession numbers in output."""
        accno_file = tempfile.NamedTemporaryFile()
        accno_file.write('FA6P1OK01CGMHQ\n')
        accno_file.seek(0)

        a = Sfffile()
        a.Parameters['-i'].on(accno_file.name)
        app_results = a(self.sff_fp)

        self._check_unmodified_sff_contents(app_results['sff'])
        app_results.cleanUp()

    def test_call_with_excluded_accession_numbers(self):
        """Sfffile should exclude specified accession numbers in output."""
        accno_file = tempfile.NamedTemporaryFile()
        accno_file.write('FA6P1OK01CGMHQ\n')
        accno_file.seek(0)

        a = Sfffile()
        a.Parameters['-e'].on(accno_file.name)
        app_results = a(self.sff_fp)

        header, reads_gen = parse_binary_sff(app_results['sff'])
        reads = list(reads_gen)

        self.assertEqual(header["number_of_reads"], 0)
        self.assertEqual(len(reads), 0)
        app_results.cleanUp()


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# test_sfffile.py

import hashlib
import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError
from cogent.app.sfffile import Sfffile


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2010, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class SfffileTests(TestCase):
    """Test the Sfffile application controller.
    """

    def setUp(self):
        test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.sff_fp = os.path.join(test_dir, 'data', 'test.sff')

    def test_exit_status(self):
        """Sfffile should raise ApplicationError if exit status is nonzero."""
        a = Sfffile()
        self.assertRaises(ApplicationError, a, 'an_sff_file_that_does_not_exist.sff')

    def test_call(self):
        """Simple sfffile call should produce expected output."""
        a = Sfffile()
        app_results = a(self.sff_fp)

        # compare md5 digests
        expected = hashlib.md5(open(self.sff_fp, 'rb').read()).hexdigest()
        observed = hashlib.md5(app_results['sff'].read()).hexdigest()
        self.assertEqual(observed, expected)

        app_results.cleanUp()

    def test_call_with_output_path(self):
        """Sfffile should store output to specified filepath."""
        _, output_fp = tempfile.mkstemp()
        a = Sfffile()
        a.Parameters['-o'].on(output_fp)
        app_results = a(self.sff_fp)

        expected = hashlib.md5(open(self.sff_fp, 'rb').read()).hexdigest()
        observed = hashlib.md5(open(output_fp, 'rb').read()).hexdigest()
        self.assertEqual(observed, expected)

        # Double check that app_results are still valid
        observed = hashlib.md5(app_results['sff'].read()).hexdigest()
        self.assertEqual(observed, expected)

        app_results.cleanUp()

    def test_call_with_included_accession_numbers(self):
        """Sfffile should include specified accession numbers in output."""
        accno_file = tempfile.NamedTemporaryFile()
        accno_file.write('FA6P1OK01CGMHQ\n')
        accno_file.seek(0)

        a = Sfffile()
        a.Parameters['-i'].on(accno_file.name)
        app_results = a(self.sff_fp)

        expected = hashlib.md5(open(self.sff_fp, 'rb').read()).hexdigest()
        observed = hashlib.md5(app_results['sff'].read()).hexdigest()
        self.assertEqual(observed, expected)

        app_results.cleanUp()

    def test_call_with_excluded_accession_numbers(self):
        """Sfffile should exclude specified accession numbers in output."""
        accno_file = tempfile.NamedTemporaryFile()
        accno_file.write('FA6P1OK01CGMHQ\n')
        accno_file.seek(0)

        a = Sfffile()
        a.Parameters['-e'].on(accno_file.name)
        app_results = a(self.sff_fp)

        expected = hashlib.md5(open(self.sff_fp, 'rb').read()).hexdigest()
        observed = hashlib.md5(app_results['sff'].read()).hexdigest()
        self.assertNotEqual(observed, expected)

        app_results.cleanUp()


if __name__ == '__main__':
    main()

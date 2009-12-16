#!/usr/bin/env python

import os, tempfile

try:
    from cogent.util.unit_test import TestCase, main
    from cogent.parse.pdb import PDBParser
    from cogent.app.stride import Stride
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.app.stride import Stride

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class StrideTest(TestCase):
    """Tests for Stride application controller"""

    def setUp(self):

        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))

    def test_stdout_input_from_entity(self):
        """Test Stride when input is an entity"""

        s = Stride()
        res = s(self.input_structure)
        self.assertEqual(res['ExitStatus'], 0)
        assert res['StdOut'] is not None
        self.assertTrue(res['StdOut'].readline().endswith('---------  2E12\n'))
        self.assertEquals(len(res['StdOut'].readlines()), 267)
        res.cleanUp()

    def test_stdout_input_from_path(self):
        """Test Stride when input is an entity"""

        s = Stride(InputHandler='_input_as_path')
        res = s(self.input_file)
        self.assertEqual(res['ExitStatus'], 0)
        assert res['StdOut'] is not None
        self.assertTrue(res['StdOut'].readline().endswith('---------  2E12\n'))
        self.assertEquals(len(res['StdOut'].readlines()), 267)
        res.cleanUp()

    def test_get_result_path(self):
        """Tests stride result path"""
        s = Stride()
        fd, name = tempfile.mkstemp()
        os.close(fd)
        s.Parameters['-f'].on(name)
        res = s(self.input_structure)
        self.assertEqual(res['ExitStatus'], 0)
        self.assertEqualItems(res.keys(), ['StdOut', 'StdErr', 'ExitStatus', 'File'])
        self.assertTrue(res['File'].readline().endswith('---------  2E12\n'))
        self.assertEquals(len(res['File'].readlines()), 267)
        res.cleanUp()
        self.assertFalse(os.path.exists(name))


if __name__ == '__main__':
    main()

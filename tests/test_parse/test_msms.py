#!/usr/bin/env python

import os, tempfile
from cogent.util.unit_test import TestCase, main
from cogent.parse.msms import parse_VertFile
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class MsmsTest(TestCase):
    """Tests for Msms application output parsers"""

    def setUp(self):
        vs = "1. 2. 3.\n" + \
             "4. 5. 6.\n" + \
             "7. 8. 9.\n"
        self.vertfile = StringIO(vs)
        
    def test_parseVertFile(self):
        out_arr = parse_VertFile(self.vertfile)
        assert out_arr.dtype == 'float64'
        assert out_arr.shape == (3,3)
        assert out_arr[0][0] == 1.
        
if __name__ == '__main__':
    main()


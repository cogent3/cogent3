#!/usr/bin/env python
#file cogent.parse.mothur.py

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2010, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"

from cStringIO import StringIO

from cogent.util.unit_test import TestCase, main
from cogent.parse.mothur import parse_otu_list

class FunctionTests(TestCase):
    def test_parse_otu_list(self):
        observed = list(parse_otu_list(StringIO(mothur_output)))
        expected = [
            (0.0, [['cccccc'], ['bbbbbb'], ['aaaaaa']]),
            (0.62, [['bbbbbb', 'cccccc'], ['aaaaaa']]),
            (0.67000000000000004, [['aaaaaa', 'bbbbbb', 'cccccc']])
            ]
        self.assertEqual(observed, expected)

mothur_output = """\
unique	3	cccccc	bbbbbb	aaaaaa
0.62	2	bbbbbb,cccccc	aaaaaa
0.67	1	aaaaaa,bbbbbb,cccccc
"""

if __name__ == '__main__':
    main()

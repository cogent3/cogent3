#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.rnaalifold import rnaalifold_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RnaalifoldParserTest(TestCase):
    """Provides tests for RNAALIFOLD RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function """
        
        #output
        self.rnaalifold_out = RNAALIFOLD
        #expected
        self.rnaalifold_exp = [['GGCUAGGUAAAUCC',[(0,13),(1,12),(3,10)],-26.50]]

        
    def test_rnaalifold_output(self):
        """Test for rnaalifold format parser"""
        
        obs = rnaalifold_parser(self.rnaalifold_out)
        self.assertEqual(obs,self.rnaalifold_exp)

RNAALIFOLD = ['GGCUAGGUAAAUCC\n', 
'((.(......).)) (-26.50 = -26.50 +   0.00) \n']

if __name__ == '__main__':
    main()

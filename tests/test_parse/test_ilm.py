#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.ilm import ilm_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class IlmParserTest(TestCase):
    """Provides tests for ILM RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function"""
        
        #output
        self.ilm_out = ILM
        #expected
        self.ilm_exp = [[(0,13),(1,12),(2,11),(6,7)]]

        
    def test_ilm_output(self):
        """Test for ilm format"""
        
        obs = ilm_parser(self.ilm_out)
        self.assertEqual(obs,self.ilm_exp)

ILM = ['\n', 'Final Matching:\n', '1 14\n', '2 13\n', '3 12\n', '4 0\n', 
'5 0\n', '6 0\n', '7 8\n', '8 7\n', '9 0\n', '10 0\n', '11 0\n', '12 3\n', 
'13 2\n', '14 1\n']

if __name__ == '__main__':
    main()

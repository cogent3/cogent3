#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.cove import coves_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class CovesParserTest(TestCase):
    """Provides tests for Coves RNA secondary structure parsers"""

    def setUp(self):
        """Setup function """
        
        #output
        self.cove_out = COVE
        #expected
        self.cove_exp = [['UAGAUGGCUCUCAU',[(0,13),(2,11),(3,10)]]]

        
    def test_cove_output(self):
        """Test coves format parser"""
        
        obs = coves_parser(self.cove_out)
        self.assertEqual(obs,self.cove_exp)

COVE =  ['coves - scoring and structure prediction of RNA sequences\n', '        using a covariance model\n', 
'     version 2.4.4, January 1996\n', '\n', 
'---------------------------------------------------\n', 
'Database to search/score:      single.fasta\n', 
'Model:                         single.fasta.cm\n', 
'GC% of background model:       50%\n', 
'---------------------------------------------------\n', '\n', 
'-32.55 bits : seq1\n', '          seq1 UAGAUGGCUCUCAU\n', 
'          seq1 >.>>......<<.<\n', '\n']

if __name__ == '__main__':
    main()

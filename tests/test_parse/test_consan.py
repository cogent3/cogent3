#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.consan import consan_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ConsanParserTest(TestCase):
    """Provides tests for CONSAN RNA secondary structure parsers"""

    def setUp(self):
        """Setup function"""
        
        #output
        self.consan_out = CONSAN
        #expected
        self.consan_exp = [{'seq1':'GGACACGUCGCUCA','seq2':'G.ACAGAUCGCUCA'},
                            [(0,13),(2,11),(5,8)]]

        
    def test_consan_output(self):
        """Test for consan format parser"""
        
        obs = consan_parser(self.consan_out)
        self.assertEqual(obs,self.consan_exp)

CONSAN = ['Using standard STA scoring\n', 
'Using QRADIUS constraints (Quality > 0.9500 SPACED 20) !\n', 
'# STOCKHOLM 1.0\n', '\n', '#=GF SC\t 1.000000\n', '#=GF PI\t 0.720000\n', 
'#=GF ME\t QRadius  Num: 3  Win: 20  Cutoff: 0.95\n', 
'\n', 'seq1                \tGGACACGUCGCUCA\n', 
'seq2                \tG.ACAGAUCGCUCA\n', 
'#=GC SS_cons\t\t    >.>..>..<..<.<\n', 
'#=GC PN \t\t        .......*......\n', '\n', '\n', '\n', 
'#=GF TIME  43.0\n', '\n', '//\n']

if __name__ == '__main__':
    main()

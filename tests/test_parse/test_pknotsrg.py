#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.parse.pknotsrg import pknotsrg_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class PknotsrgParserTest(TestCase):
    """Provides tests for pknotsRG RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function"""
        
        #output
        self.pknotsrg_out = PKNOTSRG
        #expected
        self.pknotsrg_exp = [['UGCAUAAUAGCUCC',[(0,8),(3,11),(5,13)],-22.40]]

        
    def test_pknotsrg_output(self):
        """Test for pknotsrg format parser"""
        
        obs = pknotsrg_parser(self.pknotsrg_out)
        self.assertEqual(obs,self.pknotsrg_exp)

PKNOTSRG = ['UGCAUAAUAGCUCC\n', '(..{.[..)..}.] (-22.40)\n'] 

if __name__ == '__main__':
    main()

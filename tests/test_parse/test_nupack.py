#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.parse.nupack   import nupack_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class NupackParserTest(TestCase):
    """Provides tests for NUPACK RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function"""
        
        #output
        self.nupack_out = NUPACK
        #expected
        self.nupack_exp = [['GGCUAGUCCCUUCU',[(0,9),(1,8),(4,13),(5,12)],-23.30]]

        
    def test_nupack_output(self):
        """Test for nupack format"""
        
        obs = nupack_parser(self.nupack_out)
        self.assertEqual(obs,self.nupack_exp)

NUPACK = ['****************************************************************\n', 'NUPACK 1.2\n', 
'Copyright 2007-2009 2003, 2004 by Robert M. Dirks & Niles A. Pierce\n', 
'California Institute of Technology\n', 'Pasadena, CA 91125 USA\n', '\n', 
'Last Modified: 03/18/2004\n', 
'****************************************************************\n', '\n', 
'\n', 'Fold.out Version 1.2: Complexity O(N^5) (pseudoknots enabled)\n', 
'Reading Input File...\n', 'Sequence Read.\n', 'Energy Parameters Loaded\n', 
'SeqLength = 14\n', 'Sequence and a Minimum Energy Structure:\n', 
'GGCUAGUCCCUUCU\n', '((..{{..))..}}\n', 'mfe = -23.30 kcal/mol\n', 
'pseudoknotted!\n']

if __name__ == '__main__':
    main()

#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.parse.ct import ct_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class CtParserTest(TestCase):
    """Provides tests for RNA secondary structure parsers"""

    def setUp(self):
        """Setup function"""
        
        #output
        self.carnac_out = CARNAC
        self.dynalign_out = DYNALIGN
        self.mfold_out = MFOLD
        self.sfold_out = SFOLD
        self.unafold_out = UNAFOLD
        self.knetfold_out = KNETFOLD
        #expected
        self.carnac_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)]]]
        self.dynalign_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)],-46.3]]
        self.mfold_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)],-23.47]]
        self.sfold_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)],-22.40]]
        self.unafold_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)],-20.5]]
        self.knetfold_exp = [['GCAGAUGGCUUC',[(0,11),(2,9),(3,8)]]]
        
    def test_carnac_output(self):
        """Test for ct_parser for carnac format"""
 
        obs = ct_parser(self.carnac_out)
        self.assertEqual(obs,self.carnac_exp)

    def test_dynalign_output(self):
        """Test for ct_parser for dynalign format"""
        
        obs = ct_parser(self.dynalign_out)
        self.assertEqual(obs,self.dynalign_exp)

    def test_mfold_output(self):
        """Test for ct_parser for mfold format"""
        
        obs = ct_parser(self.mfold_out)
        self.assertEqual(obs,self.mfold_exp)
        
    def test_sfold_output(self):
        """Test for ct_parser for sfold format"""
    
        obs = ct_parser(self.sfold_out)
        self.assertEqual(obs,self.sfold_exp)
        
    def test_unafold_output(self):
        """Test for ct_parser for unafold format"""
        
        obs = ct_parser(self.unafold_out)
        self.assertEqual(obs,self.unafold_exp)

    def test_knetfold_output(self):
        """Test for ct_parser for knetfold format"""
        
        obs = ct_parser(self.knetfold_out)
        self.assertEqual(obs,self.knetfold_exp)

CARNAC   = ['   12 seq1\n', '    1 G       0    2   12    1\n', 
'    2 C       1    3    0    2\n', '    3 A       2    4   10    3\n', 
'    4 G       3    5    9    4\n', '    5 A       4    6    0    5\n', 
'    6 U       5    7    0    6\n', '    7 G       6    8    0    7\n', 
'    8 G       7    9    0    8\n', '    9 C       8   10    4    9\n', 
'   10 U       6   11    3   10   \n', '   11 U       7   12    0   11\n', 
'   12 C       8   13    1   12\n']
DYNALIGN = ['   72   ENERGY = -46.3  seq 3\n', 
'    1 G       0    2   12    1\n', '    2 C       1    3    0    2\n', 
'    3 A       2    4   10    3\n', '    4 G       3    5    9    4\n', 
'    5 A       4    6    0    5\n', '    6 U       5    7    0    6\n', 
'    7 G       6    8    0    7\n', '    8 G       7    9    0    8\n', 
'    9 C       8   10    4    9\n', '   10 U       6   11    3   10   \n', 
'   11 U       7   12    0   11\n', '   12 C       8   13    1   12\n']
MFOLD    = ['   12   dG = -23.47    [initially -22.40] seq1 \n', 
'    1 G       0    2   12    1\n', '    2 C       1    3    0    2\n', 
'    3 A       2    4   10    3\n', '    4 G       3    5    9    4\n', 
'    5 A       4    6    0    5\n', '    6 U       5    7    0    6\n', 
'    7 G       6    8    0    7\n', '    8 G       7    9    0    8\n', 
'    9 C       8   10    4    9\n', '   10 U       6   11    3   10   \n', 
'   11 U       7   12    0   11\n', '   12 C       8   13    1   12\n']
SFOLD    = ['Structure    1     -22.40       0.63786E-01\n', 
'    1 G       0    2   12    1\n', '    2 C       1    3    0    2\n', 
'    3 A       2    4   10    3\n', '    4 G       3    5    9    4\n', 
'    5 A       4    6    0    5\n', '    6 U       5    7    0    6\n', 
'    7 G       6    8    0    7\n', '    8 G       7    9    0    8\n', 
'    9 C       8   10    4    9\n', '   10 U       6   11    3   10   \n', 
'   11 U       7   12    0   11\n', '   12 C       8   13    1   12\n']
UNAFOLD  = ['12\tdG = -20.5\tseq1\n', '    1 G       0    2   12    1\n', 
'    2 C       1    3    0    2\n', '    3 A       2    4   10    3\n', 
'    4 G       3    5    9    4\n', '    5 A       4    6    0    5\n', 
'    6 U       5    7    0    6\n', '    7 G       6    8    0    7\n', 
'    8 G       7    9    0    8\n', '    9 C       8   10    4    9\n', 
'   10 U       6   11    3   10   \n', '   11 U       7   12    0   11\n', 
'   12 C       8   13    1   12\n']
KNETFOLD = ['   12   \n', '    1 G       0    2   12    1\n', 
'    2 C       1    3    0    2\n', '    3 A       2    4   10    3\n', 
'    4 G       3    5    9    4\n', '    5 A       4    6    0    5\n', 
'    6 U       5    7    0    6\n', '    7 G       6    8    0    7\n', 
'    8 G       7    9    0    8\n', '    9 C       8   10    4    9\n', 
'   10 U       6   11    3   10   \n', '   11 U       7   12    0   11\n', 
'   12 C       8   13    1   12\n']


if __name__ == '__main__':
    main()



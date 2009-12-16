#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.rnaalifold import rnaalifold_parser, MinimalRnaalifoldParser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
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
        
        #output 2
        self.rnaalifold_out_2 = RNAALIFOLD_2_STRUCTS
        #expected
        self.rnaalifold_exp_2 = \
            [['GGCUAGGUAAAUCC',[(0,13),(1,12),(3,10)],\
                                float(-26.50)],\
            ['-GAUCCUAAGCGACGAAGUUYAWSCU------YGKRYARYRWWKKR-',\
                [(6,21),(7,20),(8,19),(9,18),(10,17)],\
                float(-0.80)]]
        
    def test_rnaalifold_output(self):
        """Test for rnaalifold format parser"""
        #Test empty lines
        self.assertEqual(rnaalifold_parser(''),[])
        
        #Test one structure
        obs = rnaalifold_parser(self.rnaalifold_out)
        self.assertEqual(obs,self.rnaalifold_exp)
        
        #Test two structures
        obs_2 = rnaalifold_parser(self.rnaalifold_out_2)
        self.assertEqual(obs_2,self.rnaalifold_exp_2)

class MinimalRnaalifoldParserTest(TestCase):
    """Provides tests for MinimalRnaalifoldParser structure format parser.
    """

    def setUp(self):
        """Setup function """
        
        #output
        self.rnaalifold_out = RNAALIFOLD
        #expected
        self.rnaalifold_exp = [['GGCUAGGUAAAUCC','((.(......).))',\
            float(-26.50)]]
        
        #output 2
        self.rnaalifold_out_2 = RNAALIFOLD_2_STRUCTS
        #expected
        self.rnaalifold_exp_2 = \
            [['GGCUAGGUAAAUCC','((.(......).))',float(-26.50)],\
            ['-GAUCCUAAGCGACGAAGUUYAWSCU------YGKRYARYRWWKKR-',\
            '......(((((......))))).........................',\
            float(-0.80)]]

        
    def test_rnaalifold_output(self):
        """Test for rnaalifold format parser"""
        #Test empty lines
        self.assertEqual(MinimalRnaalifoldParser(''),[])
        
        #Test one structure
        obs = MinimalRnaalifoldParser(self.rnaalifold_out)
        self.assertEqual(obs,self.rnaalifold_exp)
        
        #Test two structures
        obs_2 = MinimalRnaalifoldParser(self.rnaalifold_out_2)
        self.assertEqual(obs_2,self.rnaalifold_exp_2)

RNAALIFOLD = ['GGCUAGGUAAAUCC\n', 
'((.(......).)) (-26.50 = -26.50 +   0.00) \n']

RNAALIFOLD_2_STRUCTS = ['GGCUAGGUAAAUCC\n', 
'((.(......).)) (-26.50 = -26.50 +   0.00) \n',\
'-GAUCCUAAGCGACGAAGUUYAWSCU------YGKRYARYRWWKKR-\n',
'......(((((......)))))......................... ( -0.80 =  -1.30 +   0.50)\n',]


if __name__ == '__main__':
    main()

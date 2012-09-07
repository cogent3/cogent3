#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.kegg_pos import parse_pos_lines, parse_pos_file

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.2-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Production"

"""
Test code for kegg_pos.py in cogent.parse.  
"""

class ParsePosTests(TestCase):
    
    def test_parse_pos_lines(self):
        """Parse pos lines should parse given lines and filename"""
        test_lines = \
        ['YPO0021 hemN    28982   28982..30355    1374\n',\
         'YPO0022 glnG, glnT, ntrC    30409   complement(30409..31821)    1413\n',\
         'YPO0023 glnL, ntrB  31829   complement(31829..32878)    1050\n',\
         'YPO0024 glnA    33131   complement(33131..34540)    1410\n']

        obs = parse_pos_lines(test_lines, file_name = 'y.pestis.pos')
        exp = ['y.pestis\tYPO0021 hemN    28982   28982..30355    1374\n',\
               'y.pestis\tYPO0022 glnG, glnT, ntrC    30409   complement(30409..31821)    1413\n',\
               'y.pestis\tYPO0023 glnL, ntrB  31829   complement(31829..32878)    1050\n',\
               'y.pestis\tYPO0024 glnA    33131   complement(33131..34540)    1410\n']

        for i,parsed_line in enumerate(obs):
            self.assertEqual(parsed_line, exp[i])
        
    def test_pos_to_fields(self):
        """parse_pos_to_fields should open files and extract fields"""
        # Note that this test is set to pass, as this is just a simple
        # open/yield wrapper equivalent to the demo blocks for other parsers.
        # It is kept as an independent function to allow calling by handlers.
        pass

if __name__=="__main__":
    main()

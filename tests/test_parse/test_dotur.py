#!/usr/bin/env python
# test_dotur.py

from cogent.util.unit_test import TestCase, main
from cogent.parse.dotur import get_otu_lists, OtuListParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class DoturParserTests(TestCase):
    """Tests for DoturParser.
    """
    
    def setUp(self):
        """setup for DoturParserTests.
        """
        
        self.otu_list_string = \
"""unique	3	a	b	c
0.00	3	a	b	c
0.59	2	a,c	b
0.78	1	a,c,b
"""
        self.otu_res_list = [
            [0.0,3,[['a'],['b'],['c']]],\
            [0.0,3,[['a'],['b'],['c']]],\
            [float(0.59),2,[['a','c'],['b']]],\
            [float(0.78),1,[['a','c','b']]],\
            ]
        
        self.otu_lists_unparsed=[\
            ['a','b','c'],
            ['a','b','c'],
            ['a,c','b'],
            ['a,c,b']
            ]

        self.otu_lists_parsed=[\
            [['a'],['b'],['c']],
            [['a'],['b'],['c']],
            [['a','c'],['b']],
            [['a','c','b']]
            ]
    
    def test_get_otu_lists_no_data(self):
        """get_otu_lists should function as expected.
        """
        self.assertEqual(get_otu_lists([]),[])
    
    def test_get_otu_lists(self):
        """get_otu_lists should function as expected.
        """
        for obs, exp in zip(self.otu_lists_unparsed,self.otu_lists_parsed):
            self.assertEqual(get_otu_lists(obs),exp)
        
    def test_otulistparser_no_data(self):
        """OtuListParser should return correct result given no data.
        """
        res = OtuListParser([])
        self.assertEqual(list(res),[])
        
    def test_otulistparser_parser(self):
        """OtuListParser should return correct result given basic output.
        """
        res = OtuListParser(self.otu_list_string.split('\n'))
        self.assertEqual(res,self.otu_res_list)


if __name__ == '__main__':
    main()
    
#!/usr/bin/env python
"""Unit tests for the PSL parser.
   Compatible with blat v.34
"""

from cogent.parse.psl import make_header, MinimalPslParser, PslToTable
from cogent.util.unit_test import TestCase, main
from cogent import LoadTable

__author__ = "Gavin Huttley, Anuj Pahwa"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight","Peter Maxwell", "Gavin Huttley", "Anuj Pahwa"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Development"

fname = 'data/test.psl'

class Test(TestCase):
    
    def test_header(self):
        """should return correct header"""
        expect = ['match', 'mis-match', 'rep. match', "N's", 'Q gap count',
            'Q gap bases', 'T gap count', 'T gap bases', 'strand', 'Q name',
            'Q size', 'Q start', 'Q end', 'T name', 'T size', 'T start',
            'T end', 'block count', 'blockSizes', 'qStarts', 'tStarts']
        
        parser = MinimalPslParser(fname)
        version = parser.next()
        header = parser.next()
        self.assertEqual(header, expect)
    
    def test_psl_to_table(self):
        table = PslToTable(fname)
    
    def test_getting_seq_coords(self):
        """get correct sequence coordinates to produce a trimmed sequence"""
        table = PslToTable(fname)
        for row in table:
            query_name = row["Q name"]
            query_strand = row["strand"]
            q_start = row["Q start"]
        
    

if __name__ == "__main__":
    main()

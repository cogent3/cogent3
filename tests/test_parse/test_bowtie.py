#!/usr/bin/env python
"""Unit tests for the bowtie default output parser.
   Compatible with bowtie version 0.12.5
"""

from cogent.parse.bowtie import BowtieOutputParser, BowtieToTable
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

fname = 'data/bowtie_output.map'
expected = [['GAPC_0015:6:1:1283:11957#0/1', '-', 'Mus', 66047927, 'TGTATATATAAACATATATGGAAACTGAATATATATACATTATGTATGTATATATGTATATGTTATATATACATA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['55:A>G', '64:C>A']],
        ['GAPC_0015:6:1:1394:18813#0/1', '+', 'Mus', 77785518, 'ATGAAATTCCTAGCCAAATGGATGGACCTGGAGGGCATCATC', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 447, []],
        ['GAPC_0015:6:1:1560:18056#0/1', '+', 'Mus', 178806665, 'TAGATAAAGGCTCTGTTTTTCATCATTGAGAAATTGTTATTTTTCTGATGTTATA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['9:T>G']],
        ['GAPC_0015:6:1:1565:19849#0/1', '+', 'Mus', 116516430, 'ACCATTTGCTTGGAAAATTGTTTTCCAGCCTTTCACTCTGAG', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 141, []],
        ['GAPC_0015:6:1:1591:17397#0/1', '-', 'Mus', 120440696, 'TCTAAATCTGTTCATTAATTAAGCCTGTTTCCATGTCCTTGGTCTTAAGACCAATCTGTTATGCGGGTGTGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 0, ['70:A>C', '71:G>T']]]

class BowtieOutputTest(TestCase):
    
    def test_parsing(self):
        """make sure that the bowtie output file is parsed properly"""
        parser = BowtieOutputParser(fname)
        header = parser.next()
        index = 0
        for row in parser:
            self.assertEqual(row, expected[index])
            index += 1
    
    def test_psl_to_table(self):
        """make sure that the table is built without any errors"""
        table = BowtieToTable(fname)
    
    def test_getting_seq_coords(self):
        """get correct information from the table"""
        table = BowtieToTable(fname)
        index = 0
        for row in table:
            query_name = row['Query Name']
            strand_direction = row['Strand Direction']
            query_offset = row['Offset']
            self.assertEqual(query_name, expected[index][0])
            self.assertEqual(strand_direction, expected[index][1])
            self.assertEqual(query_offset, expected[index][3])
            index += 1
        
    def test_no_row_converter(self):
        """setting row_converter=None returns strings"""
        # straight parser
        parser = BowtieOutputParser(fname, row_converter=None)
        header = parser.next()
        for index, row in enumerate(parser):
            query_offset = row[3]
            other_matches = row[6]
            self.assertEqual(query_offset, str(expected[index][3]))
            self.assertEqual(other_matches, str(expected[index][6]))
        
        # table
        table = BowtieToTable(fname, row_converter=None)
        for index, row in enumerate(table):
            query_offset = row['Offset']
            other_matches = row['Other Matches']
            self.assertEqual(query_offset, str(expected[index][3]))
            self.assertEqual(other_matches, str(expected[index][6]))
        

if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Tests of Illumina sequence file parser.
"""
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename
from cogent.parse.illumina_sequence import (MinimalIlluminaSequenceParser)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Production"

class ParseIlluminaSequenceTests(TestCase):
    """ Test of top-level Illumina parsing functions """
    
    def setUp(self):
        """ """
        self.illumina_read1 = illumina_read1
        self.illumina_read2 = illumina_read2
        self.expected_read1 = expected_read1
        self.expected_read2 = expected_read2
        
        self.illumina_read1_fp = get_tmp_filename(
         prefix='ParseIlluminaTest',suffix='.txt')
        open(self.illumina_read1_fp,'w').write('\n'.join(self.illumina_read1))
        self.files_to_remove = [self.illumina_read1_fp]
    
    def tearDown(self):
        """ """
        remove_files(self.files_to_remove)
        
    def test_MinimalIlluminaSequenceParser(self):
        """ MinimalIlluminaSequenceParser functions as expected """
        actual_read1 = list(MinimalIlluminaSequenceParser(self.illumina_read1))
        self.assertEqual(actual_read1,self.expected_read1)
        
        actual_read2 = list(MinimalIlluminaSequenceParser(self.illumina_read2))
        self.assertEqual(actual_read2,self.expected_read2)
        
    def test_MinimalIlluminaSequenceParser_handles_filepath_as_input(self):
        """ MinimalIlluminaSequenceParser functions with filepath as input 
        """
        actual_read1 = list(MinimalIlluminaSequenceParser(
                            self.illumina_read1_fp))
        self.assertEqual(actual_read1,self.expected_read1)
        
    def test_MinimalIlluminaSequenceParser_handles_file_as_input(self):
        """ MinimalIlluminaSequenceParser functions with file handle as input
        """
        actual_read1 = list(MinimalIlluminaSequenceParser(
                            open(self.illumina_read1_fp)))
        self.assertEqual(actual_read1,self.expected_read1)

illumina_read1 = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/1:GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')

expected_read1 = [(["HWI-6X_9267","1","1","4","1699#ACCACCC/1"],
 "TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA",
 "abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB"),
 (["HWI-6X_9267","1","1","4","390#ACCTCCC/1"],
 "GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA",
 "aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa")]

illumina_read2 = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')

expected_read2 = [(["HWI-6X_9267","1","1","4","1699#ACCACCC/2"],
 "TTTTAAAAAAAAGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC",
 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB"),
 (["HWI-6X_9267","1","1","4","390#ACCTCCC/2"],
 "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
 "aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb")]
 
if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Unit tests for GFF and related parsers.
"""
from cogent.parse.gff import *
from cogent.util.unit_test import TestCase, main
from StringIO import StringIO

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2009 2008, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

headers = [
"""##gff-version 2 
##source-version <source> <version text> 
##date <date> 
##Type <type> [<seqname>] 
##DNA <seqname>
##acggctcggattggcgctggatgatagatcagacgac
##...
##end-DNA
""",
"""##gff-version 2
""",
"",
]

#    '<seqname>\t<source>\t<feature>\t<start>\t<end>\t<score>\t<strand>\t<frame>\t[attribute]\n'

data_lines = [
('seq1\tBLASTX\tsimilarity\t101\t235\t87.1\t+\t0\tTarget "HBA_HUMAN" 11 55 ; E_value 0.0003\n',
('seq1', 'BLASTX', 'similarity', 100, 235, '87.1', '+', '0', 'Target "HBA_HUMAN" 11 55 ; E_value 0.0003', None)),
('dJ102G20\tGD_mRNA\tcoding_exon\t7105\t7201\t.\t-\t2\tSequence "dJ102G20.C1.1"\n',
('dJ102G20', 'GD_mRNA', 'coding_exon', 7201, 7104, '.', '-', '2', 'Sequence "dJ102G20.C1.1"', None)),
('dJ102G20\tGD_mRNA\tcoding_exon\t7105\t7201\t.\t-\t2\t\n',
('dJ102G20', 'GD_mRNA', 'coding_exon', 7201, 7104, '.', '-', '2', '', None)),
('12345\tSource with spaces\tfeature with spaces\t-100\t3600000000\t1e-5\t-\t.\tSequence "BROADO5" ; Note "This is a \\t tab containing \\n multi line comment"\n',
('12345', 'Source with spaces', 'feature with spaces', 3600000000L, 101, '1e-5', '-', '.', 'Sequence "BROADO5" ; Note "This is a \\t tab containing \\n multi line comment"', None)),
]

class GffTest(TestCase):
    """Setup data for all the GFF parsers."""
    def testGffParserData(self):
        """Test GffParser with valid data lines"""
        for (line,canned_result) in data_lines:
            result = GffParser(StringIO(line)).next()
            self.assertEqual(result,canned_result)
            
    def testGffParserHeaders(self):
        """Test GffParser with valid data headers"""
        data = "".join([x[0] for x in data_lines])
        for header in headers:
            result = list(GffParser(StringIO(header+data)))
            self.assertEqual(result,[x[1] for x in data_lines])
            
    def test_parse_attributes(self):
        """Test parse_attributes"""
        self.assertEqual([parse_attributes(x[1][8]) for x in data_lines],
                    ['HBA_HUMAN', 'dJ102G20.C1.1', '', 'BROADO5'])
                               
if __name__ == '__main__':
    main()

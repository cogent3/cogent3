#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.cut import cut_parser

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class CutParserTest(TestCase):
    """Provides tests for codon usage table parser"""
    def test_cut_parser(self):
        """cut_parser should work on first few lines of supplied file"""
        lines = """#Species: Saccharomyces cerevisiae
#Division: gbpln
#Release: TranstermMay1994
#CdsCount: 24

#Coding GC 44.99%
#1st letter GC 47.28%
#2nd letter GC 40.83%
#3rd letter GC 46.86%

#Codon AA Fraction Frequency Number
GCA    A     0.010     1.040      6
GCC    A     0.240    22.420    130
GCG    A     0.000     0.000      0
GCT    A     0.750    71.610    411
TGC    C     0.070     0.520      3
"""
        result = cut_parser(lines.splitlines())
        self.assertEqual(result, \
            {'GCA':6,'GCC':130,'GCG':0,'GCT':411,'TGC':3})
if __name__ == '__main__':
    main()



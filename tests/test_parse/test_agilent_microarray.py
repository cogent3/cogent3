#!/usr/bin/env python
"""Tests for the Microarray output parser
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.parse.agilent_microarray import *

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class MicroarrayParserTests(TestCase):
    """Tests for MicroarrayParser.
    """

    def setUp(self):
        """Setup function for MicroarrayParser tests.
        """
        self.sample_file = ['first line in file',
                            'second line, useless data',
                            'FEATURES\tFirst\tL\tProbeName\tGeneName\tLogRatio',
                            'DATA\tFirst\tData\tProbe1\tGene1\t0.02',
                            'DATA\tSecond\tData\tProbe2\tGene2\t-0.34']
    def test_MicroarrayParser_empty_list(self):
        #Empty list should return tuple of empty lists
        self.assertEqual(MicroarrayParser([]),([],[],[]))
    def test_MicroarrayParser(self):
        #Given correct file format, return correct results
        self.assertEqual(MicroarrayParser(self.sample_file),
                            (['PROBE1','PROBE2'],
                            ['GENE1','GENE2'],[float(0.02),float(-0.34)]))
        

#run if called from command-line
if __name__ == "__main__":
    main()

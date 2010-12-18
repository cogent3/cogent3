#!/usr/bin/env python
# test_infernal.py

from cogent.util.unit_test import TestCase, main
from cogent.parse.infernal import CmsearchParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class CmsearchParserTests(TestCase):
    """Tests for CmsearchParser.
    """
    
    def setUp(self):
        """setup for CmsearchParserTests.
        """
        self.basic_results_empty = """# command data
# date
# CM summary
# Post search summary
"""
        self.basic_results_hits = """# command data
# date
# CM summary
  Model_1   Target_1  1   10  1   10  25.25   -   50
  Model_1   Target_2  3   13  1   10  14.2    -   49
# Post search summary
        """
        self.basic_res = [['Model_1','Target_1', 1, 10, 1, 10, 25.25, '-', 50],\
                          ['Model_1','Target_2', 3, 13, 1, 10, 14.2, '-',  49]]
        
        self.search_res = [['model1.cm','seq_0', 5, 23, 1, 19, 12.85, '-', 37],\
                           ['model1.cm','seq_1', 1, 19, 1, 19, 14.36, '-', 47]]
    
    def test_cmsearch_parser_no_data(self):
        """CmsearchParser should return correct result given no data.
        """
        parser = CmsearchParser([])
        self.assertEqual(list(parser),[])
        
    def test_cmsearch_parser_no_res(self):
        """CmsearchParser should return correct result given no hits in result.
        """
        parser = CmsearchParser(self.basic_results_empty.split('\n'))
        self.assertEqual(list(parser),[])
        
    def test_cmsearch_parser_basic(self):
        """CmsearchParser should return correct result given basic output.
        """
        parser = CmsearchParser(self.basic_results_hits.split('\n'))
        self.assertEqual(list(parser),self.basic_res)
        
    def test_cmsearch_parser_full(self):
        """CmsearchParser should return correct result given cmsearch output.
        """
        parser = CmsearchParser(SEARCH_DATA.split('\n'))
        self.assertEqual(list(parser),self.search_res)


SEARCH_DATA = """# command:    cmsearch -T 0.0 --tabfile /tmp/tmpQGr0PGVeaEvGUkw2TM3e.txt --informat FASTA /tmp/tmp40hq0MqFPLn2lAymQeAD.txt /tmp/tmplTEQNgv0UA7sFSV0Z2RL.txt
# date:       Mon Nov  8 13:51:12 2010
# num seqs:   3
# dbsize(Mb): 0.000124
#
# Pre-search info for CM 1: model1.cm
#
# rnd  mod  alg  cfg   beta  bit sc cut
# ---  ---  ---  ---  -----  ----------
#   1  hmm  fwd  loc      -        3.00
#   2   cm  cyk  loc  1e-10        0.00
#   3   cm  ins  loc  1e-15        0.00
#
# CM: model1.cm
#                                                         target coord   query coord                         
#                                               ----------------------  ------------                         
# model name                       target name       start        stop  start   stop    bit sc   E-value  GC%
# -------------------------------  -----------  ----------  ----------  -----  -----  --------  --------  ---
  model1.cm  seq_0                 5          23      1     19     12.85         -   37
  model1.cm  seq_1                 1          19      1     19     14.36         -   47
#
# Post-search info for CM 1: /tmp/tmpWmLUo5hsKH6nyib4nGMq.cm
#
# rnd  mod  alg  cfg   beta  bit sc cut  num hits  surv fract
# ---  ---  ---  ---  -----  ----------  --------  ----------
#   1  hmm  fwd  loc      -        3.00         2      0.4113
#   2   cm  cyk  loc  1e-10        0.00         2      0.4113
#   3   cm  ins  loc  1e-15        0.00         2      0.3065
#
#    run time
# -----------
#    00:00:00"""


if __name__ == '__main__':
    main()
    

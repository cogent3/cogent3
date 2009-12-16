#!/usr/bin/env python
# test_infernal.py

from cogent.util.unit_test import TestCase, main
from cogent.parse.infernal import CmsearchParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
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
  Target_1  1   10  1   10  25.25   -   50
  Target_2  3   13  1   10  14.2    -   49
# Post search summary
        """
        self.basic_res = [['Target_1', 1, 10, 1, 10, 25.25, '-', 50],\
                          ['Target_2', 3, 13, 1, 10, 14.2, '-',  49]]
        
        self.search_res = [['seq_a', 5, 23, 1, 19, 12.85, '-', 37],\
                           ['seq_b', 1, 19, 1, 19, 14.36, '-', 47]]
    
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


SEARCH_DATA = """# command:    cmsearch --tabfile aln1_seqs2_searchres.txt aln1.cm seqs2.fasta
# date:       Tue Oct 7 17:11:59 2008
# num seqs:   6
# dbsize(Mb): 0.000230
#
# Pre-search info for CM 1: aln1-1
#
# rnd  mod  alg  cfg   beta  bit sc cut
# ---  ---  ---  ---  -----  ----------
#   1  hmm  fwd  loc      -        3.00
#   2   cm  cyk  loc  1e-07        0.00
#   3   cm  ins  loc  1e-15        0.00
#
# CM: aln1-1
#                        target coord   query coord                         
#              ----------------------  ------------                         
# target name       start        stop  start   stop    bit sc   E-value  GC%
# -----------  ----------  ----------  -----  -----  --------  --------  ---
  seq_a                 5          23      1     19     12.85         -   37
  seq_b                 1          19      1     19     14.36         -   47
#
# Post-search info for CM 1: aln1-1
#
# rnd  mod  alg  cfg   beta  bit sc cut  num hits  surv fract
# ---  ---  ---  ---  -----  ----------  --------  ----------
#   1  hmm  fwd  loc      -        3.00         2      0.2217
#   2   cm  cyk  loc  1e-07        0.00         2      0.2217
#   3   cm  ins  loc  1e-15        0.00         2      0.1652
#
#    run time
# -----------
#    00:00:00"""


if __name__ == '__main__':
    main()
    
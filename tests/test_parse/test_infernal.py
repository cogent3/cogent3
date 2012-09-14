#!/usr/bin/env python
# test_infernal.py

from cogent.util.unit_test import TestCase, main
from cogent.parse.infernal import CmsearchParser,CmalignScoreParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
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

class CmalignScoreParserTests(TestCase):
    """Tests for CmalignScoreParser.
    """
    
    def setUp(self):
        """setup for CmalignScoreParserTests.
        """

        self.basic_results_hits = """# command: data
# date:
#
# cm summary
        1  Target_1     83     55.02      2.94     0.956  00:00:00.01
        2  Target_2     84     53.31      4.42     0.960  00:00:00.01
# post alignment summary
        """
        self.basic_res = [[1,'Target_1',83,55.02,2.94,0.956,'00:00:00.01'],\
                          [2,'Target_2',84,53.31,4.42,0.960,'00:00:00.01']]
        
        self.search_res = \
            [[1,'AABL01002928.1/2363-2445',83,55.02,2.94,0.956,'00:00:00.01'],\
             [2,'AACV01025780.1/26051-26134',84,53.31,4.42,0.960,'00:00:00.01']]
    
    def test_cmalign_score_parser_no_data(self):
        """CmalignScoreParser should return correct result given no data.
        """
        parser = CmalignScoreParser([])
        self.assertEqual(list(parser),[])
        
    def test_cmalign_score_parser_basic(self):
        """CmalignScoreParser should return correct result given basic output.
        """
        parser = CmalignScoreParser(self.basic_results_hits.split('\n'))
        self.assertEqual(list(parser),self.basic_res)
        
    def test_cmalign_score_parser_full(self):
        """CmalignScoreParser should return correct result given cmalign output.
        """
        parser = CmalignScoreParser(ALIGN_DATA.split('\n'))
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

ALIGN_DATA = """# cmalign :: align sequences to an RNA CM
# INFERNAL 1.0.2 (October 2009)
# Copyright (C) 2009 HHMI Janelia Farm Research Campus
# Freely distributed under the GNU General Public License (GPLv3)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# command: cmalign orig_alignments/RF01057_3NPQ_chain_A_results.cm RF01057_two_seqs.fasta
# date:    Thu May 1 13:00:32 2011
#
# cm name                    algorithm  config  sub  bands     tau
# -------------------------  ---------  ------  ---  -----  ------
# RF01057_3NPQ_chain_A_resu    opt acc  global   no    hmm   1e-07
#
#                                                 bit scores                           
#                                             ------------------                       
# seq idx  seq name                      len     total    struct  avg prob      elapsed
# -------  --------------------------  -----  --------  --------  --------  -----------
        1  AABL01002928.1/2363-2445       83     55.02      2.94     0.956  00:00:00.01
        2  AACV01025780.1/26051-26134     84     53.31      4.42     0.960  00:00:00.01

# STOCKHOLM 1.0
#=GF AU Infernal 1.0.2

AABL01002928.1/2363-2445   CGCGCCGAGGAGCGCUGCGACGGCCCG...UCGAGGGCCGCCAGGCUCGG
AACV01025780.1/26051-26134 CCUGCCGAGGGGCGCUGCGACCGGAUCcaaUGAGGCCCGGCCAGGCUCGG
#=GC SS_cons               :::::::::<--<<<--<--<<_____...________>>->--------
#=GC RF                    CuuuCCGAGGAGCGCUGcAACgGgcuc...uuacggcccGCcAGGCUCGG

AABL01002928.1/2363-2445   CGGGG...ACAAucgguUUUCCAACGGCGSUCUGUUUAU
AACV01025780.1/26051-26134 UAAGGuggCUUU.....GUAACAACGGCGCCCGGCUAGA
#=GC SS_cons               -----...----.....--------->>>-->:::::::
#=GC RF                    aaagG...uaaa.....ccuaCAACGGCGCUCAcuCaca
//
#
# CPU time: 0.02u 0.00s 00:00:00.02 Elapsed: 00:00:00"""

if __name__ == '__main__':
    main()
    

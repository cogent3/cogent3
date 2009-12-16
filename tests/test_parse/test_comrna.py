#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.parse.comrna import comRNA_parser,common

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ComrnaParserTest(TestCase):
    """Provides tests for COMRNA RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function """
        
        #output
        self.comrna_out = COMRNA
        #expected
        self.comrna_exp = [['GGCUAGAUAGCUCA',[(0,13),(1,12),(4,9),(5,8)]],
                           ['GGCUAGAUAGCUCA',[(0,13),(1,12),(4,9),(5,8)]],
                           ['GGCUAGAUAGCUCA',[(0,13),(1,12),(4,9),(5,8)]],
                           ['GGCUAGAUAGCUCA',[(0,13),(1,12)]]]

        
    def test_comrna_output(self):
        """Test for comrna format parser"""
        
        obs = comRNA_parser(self.comrna_out)
        self.assertEqual(obs,self.comrna_exp)

    def test_common_func(self):
        """Test common function in comrna parser """
        obs = common(self.comrna_exp)
        exp = [['GGCUAGAUAGCUCA',[(0,13),(1,12),(4,9),(5,8)]],
               ['GGCUAGAUAGCUCA',[(0,13),(1,12)]]]
        self.assertEqual(obs,exp)

COMRNA = ['comRNA input.fasta \n', '\n', 'PARAMETERS: \n', 'L  =     4,   Minimum length of a straight stem;\n', 'E  = -5.00,   Maximum stem energy allowed for a stem to be analyzed, in kcal/mol;\n', 'S  =  0.00,   Minimum stem similarity score b/w two stems compared;\n', 'Sh =  0.60,   Maximum stem similarity score threshold that will be tested;\n', 'Sl =  0.20,   Minimum stem similarity score threshold that will be tested;\n', 'P  =  0.50,   Minimum percentage of sequences in which a common structure should occur;\n', 'n  =    10,   Number of common structures to be reported;\n', 'x  =   999,   Maximum number of pseudoknot crossover pattern allowed between one stem and other stems in a structure;\n', 'a  =     1,   Use anchor region during stem comparison;\n', 'o  =     4,   Maximum number of overlapping nucleotides allowed between two stems;\n', 'c  =  0.30,   Maximum percentage of stem length that is allowed overlapping between two stems;\n', 'j  =  0.70,   Maximum percentage of stems allowed overlapping between two different cliques;\n', 'r  =  0.40,   Minimum percentage of stems required to be same for two cliques to be considered same when reporting structures;\n', 'f  =    10,   Number of flanking nucleotides of a stem to be refolded together during structure refinement;\n', 'v  =     5,   Maximum length of nucleotides allowed for a new loop to deviate from its length in the original structure pattern;\n', 'g  =     0,   Use topological sort to assemble stem blocks;\n', '\n', 'Sequence file name: "input.fasta"\n', '\n', 'Sequences loaded ...\n', '1    seq1                   72 nt\n', '2    seq2                   72 nt\n', '3    seq3                   72 nt\n', '4    seq4                   72 nt\n', '\n', '\n', 'Number of stems in each energy bin for each sequence:\n', '\n', 'energy[kc/m]            -10   -9   -8   -7   -6   -5   -4   -3   -2   -1    0\n', 'seq1                    2    1    1    6    1    5    9    2    1    3    0    1\n', 'seq2                    2    1    1    6    1    5    9    2    1    3    0    1\n', 'seq3                    2    1    1    6    1    5    9    2    1    3    0    1\n', 'seq4                    2    1    1    6    1    5    9    2    1    3    0    1\n', '\n', '\n', 'Pairwise Sequence Identity (%): \n', '\n', '      1   2   3   4\n', '\n', ' 1   -  100 100 100\n', ' 2  100  -  100 100\n', ' 3  100 100  -  100\n', ' 4  100 100 100  - \n', '\n', 'Average Pairwise Sequence Identity (%): 100\n', '\n', 'Comparing stems pairwise ... \n', '\n', 'Number of edges that has stem-similarity-score higher than a certain threshold in the stem graph:\n', '\n', 'Score:           0.8  0.78  0.76  0.74  0.72   0.7  0.68  0.66  0.64  0.62   0.6  0.58  0.56  0.54  0.52   0.5  0.48  0.46  0.44  0.42   0.4  0.38  0.36  0.34  0.32   0.3  0.28  0.26  0.24  0.22   0.2\n', 'Num of edges:     12    12    18    18    18    24    24    54    72    78   102   114   114   120   132   132   144   156   168   168   174   180   180   180   180   180   180   180   180   180   180   180\n', '\n', 'Time spent on comparing stems: 0.03 seconds user CPU time; 0.04 seconds real time.\n', '\n', 'Maximum structure finding time: 1 min\n', '\n', '===========================  S = 0.6  ===========================\n', '\n', 'Find conserved stems (cliques) ... ==== 17 cliques ==== 17 unique ====\n', 'Time spent on finding conserved stems: 0 sec CPU time; 0 sec clock time.\n', '\n', 'Construct clique topological graph ... ==== 53 edges ====\n', 'Assemble conserved stems (cliques) ... ==== 44 structures ====\n', 'Time spent on topologically assembling conserved stems: 0 sec CPU time; 0 sec clock time.\n', '\n', 'Report top 10 structures.\n', '-------------------------------------------\n', 'Structure #1: Score = 10.1, pattern: 41, path: 0 1 3 , comseq: 1 2 3 4 , incompatible_seq: 0() 1() 3() \n', '(a)  Clique 0: OriginalScore = 3.82, ModifiedScore = 3.82\n', '  1, seq1                    1 GGCUAGA   7 ...  66 UCUGGCC  72  [-13 kc/m]\n', '  2, seq2                    1 GGCUAGA   7 ...  66 UCUGGCC  72  [-13 kc/m]\n', '  3, seq3                    1 GGCUAGA   7 ...  66 UCUGGCC  72  [-13 kc/m]\n', '  4, seq4                    1 GGCUAGA   7 ...  66 UCUGGCC  72  [-13 kc/m]\n', '(b)  Clique 1: OriginalScore = 3.45, ModifiedScore = 3.45\n', '  1, seq1                   29 GGAUUGAA  36 ...  54 UUCGAUCC  61  [-11.6 kc/m]\n', '  2, seq2                   29 GGAUUGAA  36 ...  54 UUCGAUCC  61  [-11.6 kc/m]\n', '  3, seq3                   29 GGAUUGAA  36 ...  54 UUCGAUCC  61  [-11.6 kc/m]\n', '  4, seq4                   29 GGAUUGAA  36 ...  54 UUCGAUCC  61  [-11.6 kc/m]\n', '(c)  Clique 3: OriginalScore = 2.82, ModifiedScore = 2.82\n', '  1, seq1                   49 GUCGG  53 UUCGAUC  61 CCGGC  65  [-8.4 kc/m]\n', '  2, seq2                   49 GUCGG  53 UUCGAUC  61 CCGGC  65  [-8.4 kc/m]\n', '  3, seq3                   49 GUCGG  53 UUCGAUC  61 CCGGC  65  [-8.4 kc/m]\n', '  4, seq4                   49 GUCGG  53 UUCGAUC  61 CCGGC  65  [-8.4 kc/m]\n', '\n', '\n', 'seq1                    1 GGCUAGAUAGCUCA 14   \n', '                          aa  bb  bb  aa\n', 'seq2                    1 GGCUAGAUAGCUCA 14  \n', '                          aa  bb  bb  aa\n', 'seq3                    1 GGCUAGAUAGCUCA 14   \n', '                          aa  bb  bb  aa\n', 'seq4                    1 GGCUAGAUAGCUCA 14   \n', '                          aa          aa\n', '\n', '\n', '-------------------------------------------']

if __name__ == '__main__':
    main()

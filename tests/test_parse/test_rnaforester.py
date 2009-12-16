#!/usr/bin/env python

from cogent.util.unit_test    import TestCase, main
from cogent.core.info         import Info
from cogent.parse.rnaforester import rnaforester_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RnaforesterParserTest(TestCase):
    """Provides tests for RNAforester RNA secondary structure format parsers"""

    def setUp(self):
        """Setup function

        application output is not always actual(real) output from 
        the application but output in the format of the application 
        output in question. this to save space and time(mine)"""
        
        #output
        self.rnaforester_out = RNAFORESTER
        #expected
        self.rnaforester_exp = [[{'seq1':'ggccacguagcucagucgguagagcaaaggacugaaaauccuugugucguugguucaauuccaaccguggccacca',
'seq2':'-gccagauagcucagucgguagagcguucgccugaaaagugaaaggucgccgguucgaucccggcucuggccacca'},
'ggccacauagcucagucgguagagcaaacgacugaaaagccaaaggucgccgguucaaucccaacccuggccacca',
[(0, 71), (1, 70), (2, 69), (3, 68), (4, 67), (5, 66), (6, 65), (9, 24), 
(10, 23), (11, 22), (12, 21), (26, 42), (27, 41), (28, 40), (29, 39), 
(31, 38), (32, 37), (48, 64), (49, 63), (50, 62), (51, 61), (52, 60)]]]

        
    def test_rnaforester_output(self):
        """Test for rnaforester format"""
        
        obs = rnaforester_parser(self.rnaforester_out)
        self.assertEqual(obs,self.rnaforester_exp)

RNAFORESTER = ['*** Scoring parameters ***\n', '\n', 
'Scoring type: similarity\n', 'Scoring parameters:\n', 'pm:   10\n', 
'pd:   -5\n', 'bm:   1\n', 'br:   0\n', 'bd:   -10\n', '\n', '\n', 
'Input string (upper or lower case); & to end for multiple alignments, @ to quit\n', 
'....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n', 
'\n', '*** Calculation ***\n', '\n', 'clustering threshold is: 0.7\n', 
'join clusters cutoff is: 0\n', '\n', 'Computing all pairwise similarities\n', 
'2,1: 0.74606\n', '\n', 'joining alignments:\n', '1,2: 0.74606 -> 1\n', 
'Calculate similarities to other clusters\n', '\n', '\n', '*** Results ***\n', 
'\n', 'Minimum basepair probability for consensus structure (-cmin): 0.5\n', 
'\n', 'RNA Structure Cluster Nr: 1\n', 'Score: 264.25\n', 'Members: 2\n', '\n', 'seq1                     ggccacguagcucagucgguagagcaaaggacugaaaauccuugugucguugguu\n', 
'seq2                     -gccagauagcucagucgguagagcguucgccugaaaagugaaaggucgccgguu\n', 
'                          ****  ******************    * *******       ****  ****\n', '\n', 
'seq1                     caauuccaaccguggccacca\n', 
'seq2                     cgaucccggcucuggccacca\n', 
'                         * ** **  *  *********\n', '\n', 
'seq1                     (((((((..((((........)))).(((((.......))))).....(((((..\n', 
'seq2                     -((((((..((((........)))).((((.((....)))))).....(((((..\n', 
'                          *****************************   **** *****************\n', '\n', 'seq1                     .....))))))))))))....\n', 
'seq2                     .....))))))))))).....\n', 
'                         **************** ****\n', '\n', '\n', 
'Consensus sequence/structure:\n', 
'                    100%  ****  ******************    * *******       ****  ****\n', '                     90%  ****  ******************    * *******       ****  ****\n', '                     80%  ****  ******************    * *******       ****  ****\n', '                     70%  ****  ******************    * *******       ****  ****\n', '                     60%  ****  ******************    * *******       ****  ****\n', '                     50% *******************************************************\n', '                     40% *******************************************************\n', '                     30% *******************************************************\n', '                     20% *******************************************************\n', '                     10% *******************************************************\n', '                         ggccacauagcucagucgguagagcaaacgacugaaaagccaaaggucgccgguu\n', '                         (((((((..((((........)))).((((.((....)))))).....(((((..\n', '                     10% *******************************************************\n', '                     20% *******************************************************\n', '                     30% *******************************************************\n', '                     40% *******************************************************\n', '                     50% *******************************************************\n', '                     60% *******************************************************\n', '                     70%  ******************************  ****  ****************\n', '                     80%  ******************************  ****  ****************\n', '                     90%  ******************************  ****  ****************\n', '                    100%  ******************************  ****  ****************\n', '\n', '                    100% * ** **  *  *********\n', '                     90% * ** **  *  *********\n', '                     80% * ** **  *  *********\n', '                     70% * ** **  *  *********\n', '                     60% * ** **  *  *********\n', '                     50% *********************\n', '                     40% *********************\n', '                     30% *********************\n', '                     20% *********************\n', '                     10% *********************\n', '                         caaucccaacccuggccacca\n', '                         .....))))))))))))....\n', '                     10% *********************\n', '                     20% *********************\n', '                     30% *********************\n', '                     40% *********************\n', '                     50% *********************\n', '                     60% *********************\n', '                     70% **************** ****\n', '                     80% **************** ****\n', '                     90% **************** ****\n', '                    100% **************** ****\n', '\n', '\n']

if __name__ == '__main__':
    main()

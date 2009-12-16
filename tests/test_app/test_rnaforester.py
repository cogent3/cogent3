#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.rnaforester import RNAforester

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RnaforesterTest(TestCase):
    """Tests for Rnaforester application controller"""

    def setUp(self):
        self.input = rnaforester_input
        
        
    def test_stdout_input_as_lines(self):
        """Test rnaforester stdout input as lines"""

        r = RNAforester(InputHandler='_input_as_lines')
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in rnaforester_stdout])

        res = r(self.input)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test rnaforester stdout input as string"""

        r = RNAforester()
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in rnaforester_stdout])

        f = open('/tmp/input.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = r('/tmp/input.fasta')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/input.fasta')

    def test_get_result_path(self):
        """Tests rnaforester result path"""

        r = RNAforester(InputHandler='_input_as_lines')
        res = r(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

rnaforester_input = ['>seq1\n', 
'GGCCACGTAGCTCAGTCGGTAGAGCAAAGGACTGAAAATCCTTGTGTCGTTGGTTCAATTCCAACCGTGGCCACCA\n', '>seq2\n', 
'GCCAGATAGCTCAGTCGGTAGAGCGTTCGCCTGAAAAGTGAAAGGTCGCCGGTTCGATCCCGGCTCTGGCCACCA\n']
rnaforester_stdout = ['*** Scoring parameters ***\n', '\n', 'Scoring type: similarity\n', 'Scoring parameters:\n', 'pm:   10\n', 'pd:   -5\n', 'bm:   1\n', 'br:   0\n', 'bd:   -10\n', '\n', '\n', 'Input string (upper or lower case); & to end for multiple alignments, @ to quit\n', '....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n', '\n', '*** Calculation ***\n', '\n', 'clustering threshold is: 0.7\n', 'join clusters cutoff is: 0\n', '\n', 'Computing all pairwise similarities\n', '2,1: 0.74606\n', '\n', 'joining alignments:\n', '1,2: 0.74606 -> 1\n', 'Calculate similarities to other clusters\n', '\n', '\n', '*** Results ***\n', '\n', 'Minimum basepair probability for consensus structure (-cmin): 0.5\n', '\n', 'RNA Structure Cluster Nr: 1\n', 'Score: 264.25\n', 'Members: 2\n', '\n', 'seq1                     ggccacguagcucagucgguagagcaaaggacugaaaauccuugugucguugguu\n', 'seq2                     -gccagauagcucagucgguagagcguucgccugaaaagugaaaggucgccgguu\n', '                          ****  ******************    * *******       ****  ****\n', '\n', 'seq1                     caauuccaaccguggccacca\n', 'seq2                     cgaucccggcucuggccacca\n', '                         * ** **  *  *********\n', '\n', 'seq1                     (((((((..((((........)))).(((((.......))))).....(((((..\n', 'seq2                     -((((((..((((........)))).((((.((....)))))).....(((((..\n', '                          *****************************   **** *****************\n', '\n', 'seq1                     .....))))))))))))....\n', 'seq2                     .....))))))))))).....\n', '                         **************** ****\n', '\n', '\n', 'Consensus sequence/structure:\n', '                    100%  ****  ******************    * *******       ****  ****\n', '                     90%  ****  ******************    * *******       ****  ****\n', '                     80%  ****  ******************    * *******       ****  ****\n', '                     70%  ****  ******************    * *******       ****  ****\n', '                     60%  ****  ******************    * *******       ****  ****\n', '                     50% *******************************************************\n', '                     40% *******************************************************\n', '                     30% *******************************************************\n', '                     20% *******************************************************\n', '                     10% *******************************************************\n', '                         ggccacauagcucagucgguagagcaaacgacugaaaagccaaaggucgccgguu\n', '                         (((((((..((((........)))).((((.((....)))))).....(((((..\n', '                     10% *******************************************************\n', '                     20% *******************************************************\n', '                     30% *******************************************************\n', '                     40% *******************************************************\n', '                     50% *******************************************************\n', '                     60% *******************************************************\n', '                     70%  ******************************  ****  ****************\n', '                     80%  ******************************  ****  ****************\n', '                     90%  ******************************  ****  ****************\n', '                    100%  ******************************  ****  ****************\n', '\n', '                    100% * ** **  *  *********\n', '                     90% * ** **  *  *********\n', '                     80% * ** **  *  *********\n', '                     70% * ** **  *  *********\n', '                     60% * ** **  *  *********\n', '                     50% *********************\n', '                     40% *********************\n', '                     30% *********************\n', '                     20% *********************\n', '                     10% *********************\n', '                         caaucccaacccuggccacca\n', '                         .....))))))))))))....\n', '                     10% *********************\n', '                     20% *********************\n', '                     30% *********************\n', '                     40% *********************\n', '                     50% *********************\n', '                     60% *********************\n', '                     70% **************** ****\n', '                     80% **************** ****\n', '                     90% **************** ****\n', '                    100% **************** ****\n', '\n', '\n']


if __name__ == '__main__':
    main()

#!/usr/bin/env python

from os import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.carnac import Carnac

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class CarnacTest(TestCase):
    """Tests for Carnac application controller"""

    def setUp(self):
        self.input = carnac_input
        
        
    def test_stdout_input_as_lines(self):
        """Test carnac stdout input as lines

        
        If error check computation time in carnac_stdout!! 
        Usually 00:00:00 but on slower systems may be different"""

        c = Carnac(InputHandler='_input_as_lines')
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in carnac_stdout])
        res = c(self.input)

        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test carnac stdout input as string

        If error check computation time in carnac_stdout!! 
        Usually 00:00:00 but on slower systems may be different"""

        c = Carnac()
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in carnac_stdout])
        f = open('/tmp/input.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = c('/tmp/input.fasta')

        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/input.fasta')

    def test_get_result_path(self):
        """Tests carnac result path"""

        c = Carnac(InputHandler='_input_as_lines')
        res = c(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',\
        'ct1','eq1','ct2','eq2','out_seq1','out_seq2','graph','align'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['ct1'] is not None
        assert res['eq1'] is not None
        assert res['out_seq1'] is not None
        res.cleanUp()

carnac_input = ['>seq1\n', 
'GGCCACGTAGCTCAGTCGGTAGAGCAAAGGACTGAAAATCCTTGTGTCGTTGGTTCAATTCCAACCGTGGCCACCA\n', '>seq2\n', 
'GCCAGATAGCTCAGTCGGTAGAGCGTTCGCCTGAAAAGTGAAAGGTCGCCGGTTCGATCCCGGCTCTGGCCACCA\n']
carnac_stdout = [' Sequences \n', '\n', 
'  sequence  1 (length   76,  gc 52): seq1\n', 
'  sequence  2 (length   75,  gc 60): seq2\n', '\n', 
' Finding all potential stems \n', '\n', 
'  sequence  1 :  12 potential stems\n', 
'  sequence  2 :  18 potential stems\n', '\n', ' Pairwise foldings \n', '\n', 
'  seq  1 / seq  2:   3 vs   3 stems\n', '\n', 
' Combination of pairwise foldings: 6 classes of stems in 3 connex components\n', '\n', 
'  sequence  1:  1 cofoldings with   12  stems ->   3 selected stems +   2 remaining stems \n', 
'  sequence  2:  1 cofoldings with   18  stems ->   3 selected stems +   1 remaining stems \n', '\n', ' Overall computation time : 00:00:00\n', '\n', 
' Parameter values :\n', '\n', '  AP_THD         8\n', 
'  SUB            -5\n', '  SIZE_HAIRPIN   3\n', '  CORRECT_THD    1\n', 
'  INI_THD        -500\n', '  DIST_1         50\n', '  DIST_2         300\n', 
'  Energy tresholds :          -800 -  -1300\n', 
'  Allowing single hairpins : no\n', '  SIZE_MAX_HP    8\n', 
'  THD_HP         -1500\n', '  FLT            1\n']


if __name__ == '__main__':
    main()

#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.consan import Consan

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ConsanTest(TestCase):
    """Tests for Consan application controller"""

    def setUp(self):
        self.input1 = consan_input1
        self.input2 = consan_input2
        
        
    def test_stdout_input_as_lines(self):
        """Test Consan stdout input as lines"""

        c = Consan(InputHandler='_input_as_lines')
        input = self.input1
        input.extend(self.input2)
        res = c(input)

        #Impossible to compare stdout since copmutation time is in the output
        #which may differ between runs
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test Consan stdout input as string"""

        c = Consan()
        f = open('/tmp/seq1.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input1])
        f.write(txt)
        f.close()
        s = open('/tmp/seq2.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input2])
        s.write(txt)
        s.close()
        res = c(['/tmp/seq1.fasta','/tmp/seq2.fasta'])

        #Impossible to compare stdout since copmutation time is in the output
        #which may differ between runs
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
                         
        res.cleanUp()

    def test_get_result_path(self):
        """Tests Consan result path"""

        c = Consan(InputHandler='_input_as_lines')
        input = self.input1
        input.extend(self.input2)
        res = c(input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

consan_input1 = ['>seq1\n', 
'GGCCACGTAGCTCAGTCGGTAGAGCAAAGGACTGAAAATCCTTGTGTCGTTGGTTCAATTCCAACCGTGGCCACCA']
consan_input2 = ['>seq2\n', 
'GCCAGATAGCTCAGTCGGTAGAGCGTTCGCCTGAAAAGTGAAAGGTCGCCGGTTCGATCCCGGCTCTGGCCACCA']



if __name__ == '__main__':
    main()

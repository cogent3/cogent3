#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.contrafold import Contrafold

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ContrafoldTest(TestCase):
    """Tests for Contrafold application controller"""

    def setUp(self):
        self.input = contrafold_input
        
        
    def test_stdout_input_as_lines(self):
        """Test contrafold stdout input as lines"""

        c = Contrafold(InputHandler='_input_as_lines')
        exp= '%s\n' % '\n'.join([str(i).strip('\n') for i in contrafold_stdout])

        res = c(self.input)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test contrafold stdout input as string"""

        c = Contrafold()
        exp= '%s\n' % '\n'.join([str(i).strip('\n') for i in contrafold_stdout])
        f = open('/tmp/single.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = c('/tmp/single.fasta')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests contrafold result path"""

        c = Contrafold(InputHandler='_input_as_lines')
        res = c(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

contrafold_input = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']
contrafold_stdout = ['1 G 0\n', '2 G 71\n', '3 C 70\n', '4 U 9\n', '5 A 0\n', 
'6 G 0\n', '7 A 0\n', '8 U 0\n', '9 A 4\n', '10 G 0\n', '11 C 0\n', 
'12 U 0\n', '13 C 0\n', '14 A 0\n', '15 G 0\n', '16 A 0\n', '17 U 0\n', 
'18 G 0\n', '19 G 0\n', '20 U 69\n', '21 A 68\n', '22 G 67\n', '23 A 66\n', 
'24 G 65\n', '25 C 64\n', '26 A 0\n', '27 G 62\n', '28 A 0\n', '29 G 61\n', 
'30 G 60\n', '31 A 59\n', '32 U 58\n', '33 U 57\n', '34 G 56\n', '35 A 55\n', 
'36 A 54\n', '37 G 51\n', '38 A 50\n', '39 U 49\n', '40 C 0\n', '41 C 0\n', 
'42 U 0\n', '43 U 0\n', '44 G 0\n', '45 U 0\n', '46 G 0\n', '47 U 0\n', 
'48 C 0\n', '49 G 39\n', '50 U 38\n', '51 C 37\n', '52 G 0\n', '53 G 0\n', 
'54 U 36\n', '55 U 35\n', '56 C 34\n', '57 G 33\n', '58 A 32\n', '59 U 31\n', 
'60 C 30\n', '61 C 29\n', '62 C 27\n', '63 G 0\n', '64 G 25\n', '65 C 24\n', 
'66 U 23\n', '67 C 22\n', '68 U 21\n', '69 G 20\n', '70 G 3\n', '71 C 2\n']


if __name__ == '__main__':
    main()

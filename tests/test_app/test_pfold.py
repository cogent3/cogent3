#!/usr/bin/env python

"""Provides Tests for pfold application controller.

IMPORTANT!!!!!
'dir' variable must be set in pfold app controller file for application to work!
"""

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.pfold      import fasta2col,findphyl,mltree,scfg

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

class Fasta2colTest(TestCase):
    """Tests for fasta2col application controller"""

    def setUp(self):
        self.input = f2c_input
        
        
    def test_stdout_input_as_lines(self):
        """Test fasta2col stdout input as lines"""

        c = fasta2col(InputHandler='_input_as_lines')
        res = c(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test fasta2col stdout input as string"""

        c = fasta2col()
        f = open('/tmp/single.col','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = c('/tmp/single.col')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.col')

    def test_get_result_path(self):
        """Tests fasta2col result path"""

        c = fasta2col(InputHandler='_input_as_lines')
        res = c(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)

        res.cleanUp()

class FindphylTest(TestCase):
    """Tests for findphyl application controller"""

    def setUp(self):
        self.input = fp_input
               
    def test_input_as_lines(self):
        """Test findphyl input as lines"""

        fp = findphyl(InputHandler='_input_as_lines')
        res = fp(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test findphyl input as string"""

        fp = findphyl()
        f = open('/tmp/single.col','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = fp('/tmp/single.col')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.col')

    def test_get_result_path(self):
        """Tests findphyl result path"""

        fp = findphyl(InputHandler='_input_as_lines')
        res = fp(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)

        res.cleanUp()

class MltreeTest(TestCase):
    """Tests for Pfold application controller"""

    def setUp(self):
        self.input = ml_input
        
        
    def test_input_as_lines(self):
        """Test mltree input as lines"""

        m = mltree(InputHandler='_input_as_lines')
        res = m(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test mltree input as string"""

        m = mltree()
        f = open('/tmp/single.col','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = m('/tmp/single.col')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.col')

    def test_get_result_path(self):
        """Tests mltree result path"""

        m = mltree(InputHandler='_input_as_lines')
        res = m(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)

        res.cleanUp()

class ScfgTest(TestCase):
    """Tests for scfg application controller"""

    def setUp(self):
        self.input = scfg_input
        
        
    def test_input_as_lines(self):
        """Test scfg input as lines"""

        s = scfg(InputHandler='_input_as_lines')
        res = s(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test scfg input as string"""

        s = scfg()
        f = open('/tmp/single.col','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = s('/tmp/single.col')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.col')

    def test_get_result_path(self):
        """Tests scfg result path"""

        s = scfg(InputHandler='_input_as_lines')
        res = s(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)

        res.cleanUp()

f2c_input   = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']
f2c_stdout  = ['; generated by fasta2col\n', 
'; ============================================================\n', 
'; TYPE              RNA\n', '; COL 1             label\n', 
'; COL 2             residue\n', '; COL 3             seqpos\n', 
'; COL 4             alignpos\n', '; ENTRY             seq1\n', 
'; ----------\n', 'N  G      1      1\n', 'N  G      2      2\n', 
'N  C      3      3\n', 'N  U      4      4\n', 'N  A      5      5\n', 
'N  G      6      6\n', 'N  A      7      7\n', 'N  U      8      8\n', 
'N  A      9      9\n', 'N  G     10     10\n', 'N  C     11     11\n', 
'N  U     12     12\n', 'N  C     13     13\n', 'N  A     14     14\n', 
'N  G     15     15\n', 'N  A     16     16\n', 'N  U     17     17\n', 
'N  G     18     18\n', 'N  G     19     19\n', 'N  U     20     20\n', 
'N  A     21     21\n', 'N  G     22     22\n', 'N  A     23     23\n', 
'N  G     24     24\n', 'N  C     25     25\n', 'N  A     26     26\n', 
'N  G     27     27\n', 'N  A     28     28\n', 'N  G     29     29\n', 
'N  G     30     30\n', 'N  A     31     31\n', 'N  U     32     32\n', 
'N  U     33     33\n', 'N  G     34     34\n', 'N  A     35     35\n', 
'N  A     36     36\n', 'N  G     37     37\n', 'N  A     38     38\n', 
'N  U     39     39\n', 'N  C     40     40\n', 'N  C     41     41\n', 
'N  U     42     42\n', 'N  U     43     43\n', 'N  G     44     44\n', 
'N  U     45     45\n', 'N  G     46     46\n', 'N  U     47     47\n', 
'N  C     48     48\n', 'N  G     49     49\n', 'N  U     50     50\n', 
'N  C     51     51\n', 'N  G     52     52\n', 'N  G     53     53\n', 
'N  U     54     54\n', 'N  U     55     55\n', 'N  C     56     56\n', 
'N  G     57     57\n', 'N  A     58     58\n', 'N  U     59     59\n', 
'N  C     60     60\n', 'N  C     61     61\n', 'N  C     62     62\n', 
'N  G     63     63\n', 'N  G     64     64\n', 'N  C     65     65\n', 
'N  U     66     66\n', 'N  C     67     67\n', 'N  U     68     68\n', 
'N  G     69     69\n', 'N  G     70     70\n', 'N  C     71     71\n', 
'; **********\n']
fp_input    = f2c_stdout
fp_stdout   = ['; generated by fasta2col\n', 
'; ============================================================\n', 
'; TYPE              TREE\n', '; COL 1             label\n', '; COL 2             number\n', '; COL 3             name\n', '; COL 4             uplen\n', 
'; COL 5             child\n', '; COL 6             brother\n', 
'; ENTRY             tree\n', '; root              1\n', '; ----------\n', 
' N     1         seq1 0.000000     .     .\n', '; **********\n', 
'; TYPE              RNA\n', '; COL 1             label\n', 
'; COL 2             residue\n', '; COL 3             seqpos\n', 
'; COL 4             alignpos\n', '; ENTRY             seq1\n', 
'; ----------\n', 'N  G      1      1\n', 'N  G      2      2\n', 
'N  C      3      3\n', 'N  U      4      4\n', 'N  A      5      5\n', 
'N  G      6      6\n', 'N  A      7      7\n', 'N  U      8      8\n', 
'N  A      9      9\n', 'N  G     10     10\n', 'N  C     11     11\n', 
'N  U     12     12\n', 'N  C     13     13\n', 'N  A     14     14\n', 
'N  G     15     15\n', 'N  A     16     16\n', 'N  U     17     17\n', 
'N  G     18     18\n', 'N  G     19     19\n', 'N  U     20     20\n', 
'N  A     21     21\n', 'N  G     22     22\n', 'N  A     23     23\n', 
'N  G     24     24\n', 'N  C     25     25\n', 'N  A     26     26\n', 
'N  G     27     27\n', 'N  A     28     28\n', 'N  G     29     29\n', 
'N  G     30     30\n', 'N  A     31     31\n', 'N  U     32     32\n', 
'N  U     33     33\n', 'N  G     34     34\n', 'N  A     35     35\n', 
'N  A     36     36\n', 'N  G     37     37\n', 'N  A     38     38\n', 
'N  U     39     39\n', 'N  C     40     40\n', 'N  C     41     41\n', 
'N  U     42     42\n', 'N  U     43     43\n', 'N  G     44     44\n', 
'N  U     45     45\n', 'N  G     46     46\n', 'N  U     47     47\n', 
'N  C     48     48\n', 'N  G     49     49\n', 'N  U     50     50\n', 
'N  C     51     51\n', 'N  G     52     52\n', 'N  G     53     53\n', 
'N  U     54     54\n', 'N  U     55     55\n', 'N  C     56     56\n', 
'N  G     57     57\n', 'N  A     58     58\n', 'N  U     59     59\n', 
'N  C     60     60\n', 'N  C     61     61\n', 'N  C     62     62\n', 
'N  G     63     63\n', 'N  G     64     64\n', 'N  C     65     65\n', 
'N  U     66     66\n', 'N  C     67     67\n', 'N  U     68     68\n', 
'N  G     69     69\n', 'N  G     70     70\n', 'N  C     71     71\n', 
'; **********\n']
ml_input    = fp_stdout
ml_stdout   = ['; generated by fasta2col\n', 
'; ============================================================\n', 
'; TYPE              TREE\n', '; COL 1             label\n', 
'; COL 2             number\n', '; COL 3             name\n', 
'; COL 4             uplen\n', '; COL 5             child\n', 
'; COL 6             brother\n', '; ENTRY             tree\n', 
'; root              1\n', '; ----------\n', 
' N     1         seq1 0.001000     .     .\n', '; **********\n', 
'; TYPE              RNA\n', '; COL 1             label\n', 
'; COL 2             residue\n', '; COL 3             seqpos\n', 
'; COL 4             alignpos\n', '; ENTRY             seq1\n', 
'; ----------\n', 'N  G      1      1\n', 'N  G      2      2\n', 
'N  C      3      3\n', 'N  U      4      4\n', 'N  A      5      5\n', 
'N  G      6      6\n', 'N  A      7      7\n', 'N  U      8      8\n', 
'N  A      9      9\n', 'N  G     10     10\n', 'N  C     11     11\n', 
'N  U     12     12\n', 'N  C     13     13\n', 'N  A     14     14\n', 
'N  G     15     15\n', 'N  A     16     16\n', 'N  U     17     17\n', 
'N  G     18     18\n', 'N  G     19     19\n', 'N  U     20     20\n', 
'N  A     21     21\n', 'N  G     22     22\n', 'N  A     23     23\n', 
'N  G     24     24\n', 'N  C     25     25\n', 'N  A     26     26\n', 
'N  G     27     27\n', 'N  A     28     28\n', 'N  G     29     29\n', 
'N  G     30     30\n', 'N  A     31     31\n', 'N  U     32     32\n', 
'N  U     33     33\n', 'N  G     34     34\n', 'N  A     35     35\n', 
'N  A     36     36\n', 'N  G     37     37\n', 'N  A     38     38\n', 
'N  U     39     39\n', 'N  C     40     40\n', 'N  C     41     41\n', 
'N  U     42     42\n', 'N  U     43     43\n', 'N  G     44     44\n', 
'N  U     45     45\n', 'N  G     46     46\n', 'N  U     47     47\n', 
'N  C     48     48\n', 'N  G     49     49\n', 'N  U     50     50\n', 
'N  C     51     51\n', 'N  G     52     52\n', 'N  G     53     53\n', 
'N  U     54     54\n', 'N  U     55     55\n', 'N  C     56     56\n', 
'N  G     57     57\n', 'N  A     58     58\n', 'N  U     59     59\n', 
'N  C     60     60\n', 'N  C     61     61\n', 'N  C     62     62\n', 
'N  G     63     63\n', 'N  G     64     64\n', 'N  C     65     65\n', 
'N  U     66     66\n', 'N  C     67     67\n', 'N  U     68     68\n', 
'N  G     69     69\n', 'N  G     70     70\n', 'N  C     71     71\n', 
'; **********\n']
scfg_input  = ml_stdout
scfg_stdout = ['; generated by fasta2col\n', 
'; ============================================================\n', 
'; TYPE              TREE\n', '; COL 1             label\n', 
'; COL 2             number\n', '; COL 3             name\n', 
'; COL 4             uplen\n', '; COL 5             child\n', 
'; COL 6             brother\n', '; ENTRY             tree\n', 
'; root              1\n', '; ----------\n', 
' N     1         seq1 0.001000     .     .\n', '; **********\n', 
'; TYPE              RNA\n', '; COL 1             label\n', 
'; COL 2             residue\n', '; COL 3             seqpos\n', 
'; COL 4             alignpos\n', '; COL 5             align_bp\n', 
'; COL 6             certainty\n', '; ENTRY             seq1\n', 
'; ----------\n', 'N  G      1      1     . 0.8723\n',
 'N  G      2      2    71 0.5212\n', 'N  C      3      3    70 0.5697\n', 
'N  U      4      4    69 0.5377\n', 'N  A      5      5    68 0.5193\n', 
'N  G      6      6    67 0.4899\n', 'N  A      7      7     . 0.6499\n', 
'N  U      8      8     . 0.9159\n', 'N  A      9      9     . 0.7860\n', 
'N  G     10     10     . 0.4070\n', 'N  C     11     11     . 0.3208\n', 
'N  U     12     12     . 0.4499\n', 'N  C     13     13     . 0.5170\n', 
'N  A     14     14     . 0.8507\n', 'N  G     15     15     . 0.8156\n', 
'N  A     16     16     . 0.8715\n', 'N  U     17     17     . 0.8722\n', 
'N  G     18     18     . 0.7489\n', 'N  G     19     19     . 0.7477\n', 
'N  U     20     20     . 0.7500\n', 'N  A     21     21     . 0.7109\n', 
'N  G     22     22     . 0.3999\n', 'N  A     23     23     . 0.3800\n', 
'N  G     24     24     . 0.3241\n', 'N  C     25     25     . 0.3072\n', 
'N  A     26     26     . 0.7266\n', 'N  G     27     27     . 0.5672\n', 
'N  A     28     28     . 0.4695\n', 'N  G     29     29    41 0.5402\n', 
'N  G     30     30    40 0.5552\n', 'N  A     31     31    39 0.5167\n', 
'N  U     32     32    38 0.4277\n', 'N  U     33     33     . 0.4967\n', 
'N  G     34     34     . 0.7119\n', 'N  A     35     35     . 0.7504\n', 
'N  A     36     36     . 0.8149\n', 'N  G     37     37     . 0.6416\n', 
'N  A     38     38    32 0.4277\n', 'N  U     39     39    31 0.5167\n', 
'N  C     40     40    30 0.5552\n', 'N  C     41     41    29 0.5402\n', 
'N  U     42     42     . 0.4754\n', 'N  U     43     43     . 0.6645\n', 
'N  G     44     44     . 0.7777\n', 'N  U     45     45     . 0.8122\n', 
'N  G     46     46     . 0.6969\n', 'N  U     47     47     . 0.6836\n', 
'N  C     48     48     . 0.5963\n', 'N  G     49     49     . 0.4748\n', 
'N  U     50     50     . 0.5209\n', 'N  C     51     51     . 0.4149\n', 
'N  G     52     52     . 0.3754\n', 'N  G     53     53     . 0.4969\n', 
'N  U     54     54     . 0.7206\n', 'N  U     55     55     . 0.7112\n', 
'N  C     56     56     . 0.5039\n', 'N  G     57     57     . 0.5171\n', 
'N  A     58     58     . 0.5683\n', 'N  U     59     59     . 0.6093\n', 
'N  C     60     60     . 0.5462\n', 'N  C     61     61     . 0.3332\n', 
'N  C     62     62     . 0.3746\n', 'N  G     63     63     . 0.3898\n', 
'N  G     64     64     . 0.3047\n', 'N  C     65     65     . 0.2899\n', 
'N  U     66     66     . 0.3037\n', 'N  C     67     67     6 0.4899\n', 
'N  U     68     68     5 0.5193\n', 'N  G     69     69     4 0.5377\n', 
'N  G     70     70     3 0.5697\n', 'N  C     71     71     2 0.5212\n', 
'; **********\n']


if __name__ == '__main__':
    main()

#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.dynalign   import Dynalign

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class DynalignTest(TestCase):
    """Tests for Dynalign application controller"""

    def setUp(self):
        self.input1 = dynalign_input1
        self.input2 = dynalign_input2
        
        
    def test_stdout_input_as_lines(self):
        """Test Dynalign stdout input as lines"""

        d = Dynalign(InputHandler='_input_as_lines')
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in dynalign_stdout])

        res = d([self.input1,self.input2])
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test Dynalign stdout input as string"""

        d = Dynalign()
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in dynalign_stdout])

        f = open('/tmp/dyn1','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input1])
        f.write(txt)
        f.close()
        s = open('/tmp/dyn2','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input2])
        s.write(txt)
        s.close()
        res = d(['/tmp/dyn1','/tmp/dyn2'])
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_get_result_path(self):
        """Tests Dynalign result path"""

        d = Dynalign(InputHandler='_input_as_lines')
        res = d([self.input1,self.input2])
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',\
        'seq_1_ct','seq_2_ct','alignment'])
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

dynalign_input1 = [';\n', 'testSeq1\n', '\n', 
'GGCTAGATAGCTCAGGTCGGTTCGATCCCGGCTCTGGCC1']
dynalign_input2 = [';\n', 'testSeq2\n', '\n', 
'GGCTAGATAGCTGTGTCGTCGGTTCGATCCCGGCTCTGGCC1']
dynalign_stdout = ['12%\n', '25%\n', '38%\n', '51%\n', '64%\n', '76%\n', 
'89%\n']


if __name__ == '__main__':
    main()

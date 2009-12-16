#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.mfold      import Mfold

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class MfoldTest(TestCase):
    """Tests for Mfold application controller"""

    def setUp(self):
        self.input = mfold_input
        
        
    def test_stdout_input_as_lines(self):
        """Test Mfold stdout input as lines"""

        m = Mfold(InputHandler='_input_as_lines')
        res = m(self.input)

        #Impossible to compare stdout since tmp filenames in app controller
        #can't be predicted and are in stdout
        #Test exitstatus = 0 and stdout is not none
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test Mfold stdout input as string"""

        m = Mfold()
        f = open('/tmp/single.fasta','w')
        f.write('\n'.join(self.input))
        f.close()
        res = m('/tmp/single.fasta')
 
        #Impossible to compare stdout since tmp filenames in app controller
        #can't be predicted and are in stdout
        #Test exitstatus = 0 and stdout is not none
        assert res['StdOut'] is not None
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests mfold result path"""

        m = Mfold(InputHandler='_input_as_lines')
        res = m(self.input)


        self.assertEqualItems(res.keys(),['ps','_1.ps','_1.ss','StdOut', 'StdErr', 'ExitStatus', 'ct_all', 'ct1', 'log', 'ann', 'h-num', 'det', 'pnt', 'sav', 'ss-count', '-local.seq', 'rnaml', 'out', 'plot','pdf1'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['ct_all'] is not None
        assert res['log'] is not None
        assert res['ann'] is not None
        assert res['h-num'] is not None
        assert res['det'] is not None
        assert res['pnt'] is not None
        assert res['sav'] is not None
        assert res['ss-count'] is not None
        assert res['-local.seq'] is not None
        assert res['rnaml'] is not None
        assert res['out'] is not None
        assert res['plot'] is not None

        res.cleanUp()

mfold_input = ['>seq1\n', 'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', '\n']


if __name__ == '__main__':
    main()

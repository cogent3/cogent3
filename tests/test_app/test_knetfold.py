#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.knetfold   import Knetfold

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class KnetfoldTest(TestCase):
    """Tests for Knetfold application controller"""

    def setUp(self):
        self.input = knetfold_input
        
        
    def test_input_as_lines(self):
        """Test Knetfold stdout input as lines"""

        k = Knetfold(InputHandler='_input_as_lines')
        res = k(self.input)
       
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None

        res.cleanUp()

    def test_input_as_string(self):
        """Test Knetfold stdout input as string"""

        k = Knetfold()

        f = open('/tmp/single.fasta','w')
        f.write('\n'.join(self.input))
        f.close()
        res = k('/tmp/single.fasta')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None

        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests knetfold result path"""

        k = Knetfold(InputHandler='_input_as_lines')
        res = k(self.input)

        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',
        'ct','coll','sec','fasta','pdf','mx0','mx1','mx2','mx3'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['ct'] is not None
        assert res['coll'] is not None
        assert res['sec'] is not None
        assert res['fasta'] is not None

        res.cleanUp()

knetfold_input = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']


if __name__ == '__main__':
    main()

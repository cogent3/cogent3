#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.rnaalifold import RNAalifold

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RnaalifoldTest(TestCase):
    """Tests for Rnaalifold application controller"""

    def setUp(self):
        self.input = RNAALIFOLD_INPUT
        
        
    def test_input_as_lines(self):
        """Test rnaalifold stdout input as lines"""

        r = RNAalifold(InputHandler='_input_as_lines')
        res = r(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test rnaalifold stdout input as string"""

        r = RNAalifold()
        f = open('/tmp/clustal','w')
        f.write('\n'.join(self.input))
        f.close()
        res = r('/tmp/clustal')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/clustal')

    def test_get_result_path(self):
        """Tests rnaalifold result path"""

        r = RNAalifold(InputHandler='_input_as_lines')
        res = r(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus','SS'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None

        res.cleanUp()

RNAALIFOLD_INPUT = ['CLUSTAL\n', '\n', 'seq1 GGCTAGATAGCTCAGATGGT-AGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC----\n']

if __name__ == '__main__':
    main()

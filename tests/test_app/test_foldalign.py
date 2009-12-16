#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.foldalign  import foldalign

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class FoldalignTest(TestCase):
    """Tests for Foldalign application controller"""

    def setUp(self):
        self.input = FOLDALIGN_INPUT
        
        
    def test_input_as_lines(self):
        """Test foldalign stdout input as lines"""

        f = foldalign(InputHandler='_input_as_lines')
        res = f(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test foldalign stdout input as string"""

        f = foldalign()
        t = open('/tmp/single.col','w')
        t.write('\n'.join(self.input))
        t.close()
        res = f('/tmp/single.col')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.col')

    def test_get_result_path(self):
        """Tests foldalign result path"""

        f = foldalign(InputHandler='_input_as_lines')
        res = f(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None

        res.cleanUp()

FOLDALIGN_INPUT = ['>seq1\n', 'GGCCACGTAGCTCAGTCGGTAGAGCAAAGGACTGAAAATCCTTGTGTCGTTGGTTCAATTCCAACCGTGGCCACCA','>seq2\n', 'GCCAGATAGCTCAGTCGGTAGAGCGTTCGCCTGAAAAGTGAAAGGTCGCCGGTTCGATCCCGGCTCTGGCCACCA']


if __name__ == '__main__':
    main()

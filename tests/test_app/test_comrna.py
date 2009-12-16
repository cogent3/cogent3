#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.comrna     import comRNA

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ComrnaTest(TestCase):
    """Tests for comRNA application controller"""

    def setUp(self):
        self.input = comrna_input
             
    def test_input_as_lines(self):
        """Test comrna input as lines"""

        c = comRNA(InputHandler='_input_as_lines')
        res = c(self.input)

        #Can't compare stdout since comRNA app controller uses tmp filenames
        #that are impossible to predict. 
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_input_as_string(self):
        """Test comrna input as string"""

        c = comRNA()
        f = open('/tmp/single.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = c('/tmp/single.fasta')
        #Can't compare stdout since comRNA app controller uses tmp filenames
        #that are impossible to predict. 
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests comrna result path"""

        c = comRNA(InputHandler='_input_as_lines')
        res = c(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

comrna_input = ['>seq1\n', 'GGCTAGATAGCTCAGATGGTAGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC\n', '>seq2\n', 'GGCTAGATAGCTCAGATGGTAGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC\n', '>seq3\n', 'GGCTAGATAGCTCAGATGGTAGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC\n', '>seq4\n', 'GGCTAGATAGCTCAGATGGTAGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC']



if __name__ == '__main__':
    main()

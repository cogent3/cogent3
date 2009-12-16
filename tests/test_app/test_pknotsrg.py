#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.pknotsrg   import PknotsRG

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

class PknotsrgTest(TestCase):
    """Tests for Pknotsrg application controller"""

    def setUp(self):
        self.input = pknotsrg_input
        
        
    def test_stdout_input_as_lines(self):
        """Test pknotsrg stdout input as lines"""

        p = PknotsRG(InputHandler='_input_as_lines')
        exp= '%s\n' % '\n'.join([str(i).strip('\n') for i in pknotsrg_stdout])

        res = p(self.input)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test pknotsrg stdout input as string"""

        p = PknotsRG()
        exp= '%s\n' % '\n'.join([str(i).strip('\n') for i in pknotsrg_stdout])
        f = open('/tmp/single.plain','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = p('/tmp/single.plain')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/single.plain')

    def test_get_result_path(self):
        """Tests pknotsrg result path"""

        p = PknotsRG(InputHandler='_input_as_lines')
        res = p(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

pknotsrg_input = ['GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n'] 
pknotsrg_stdout = ['GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', '.((((((..((((........)))).(((((((...))))))).....(((((.......)))))))))))  (-22.40)\n']


if __name__ == '__main__':
    main()

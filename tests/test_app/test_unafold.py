#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.unafold    import hybrid_ss_min

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class UnafoldTest(TestCase):
    """Tests for Unafold application controller"""

    def setUp(self):
        self.input = unafold_input
        
    def test_stdout_input_as_lines(self):
        """Test Unafold stdout input as lines"""

        u = hybrid_ss_min(InputHandler='_input_as_lines')
        exp = '\n'.join(unafold_stdout)

        res = u(self.input)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test Unafold stdout input as string"""

        u = hybrid_ss_min()
        exp = '\n'.join(unafold_stdout)

        f = open('/tmp/single.fasta','w')
        f.write('\n'.join(self.input))
        f.close()
        res = u('/tmp/single.fasta')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests unafold result path"""

        u = hybrid_ss_min(InputHandler='_input_as_lines')
        res = u(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',\
        'ct','dG','run','plot_37','ext_37'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['ct'] is not None
        assert res['dG'] is not None
        assert res['run'] is not None
        assert res['plot_37'] is not None
        assert res['ext_37'] is not None
        res.cleanUp()

unafold_input = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']
unafold_stdout = ['Calculating for seq1, t = 37\n']


if __name__ == '__main__':
    main()

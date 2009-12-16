#!/usr/bin/env python

"""Provides Tests for Sfold application controller.

IMPORTANT!!! don't forget to set param_dir variable in sfold
application controller IMPORTANT!!!
"""

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.core.info      import Info
from cogent.app.sfold      import Sfold

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class SfoldTest(TestCase):
    """Tests for Sfold application controller"""

    def setUp(self):
        self.input = sfold_input
        
    def test_stdout_input_as_lines(self):
        """Test Sfold stdout input as lines"""

        s = Sfold(InputHandler='_input_as_lines')
        res = s(self.input)

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """Test Sfold stdout input as string"""

        s = Sfold()
        f = open('/tmp/single.fasta','w')
        f.write('\n'.join(self.input))
        f.close()
        res = s('/tmp/single.fasta')

        self.assertEqual(res['ExitStatus'],0)
        assert res['StdOut'] is not None
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests sfold result path"""

        s = Sfold(InputHandler='_input_as_lines')
        res = s(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',
        '10structure','10structure_2','Dharmacon_thermo','bp',
        'bprob','cdf','fe','loopr','oligo','oligo_f','pdf',
        'sample','sample_1000','sirna','sirna_f','sirna_s',
        'smfe','sstrand','stability'])
        self.assertEqual(res['ExitStatus'],0)
        assert res['10structure'] is not None
        assert res['10structure_2'] is not None
        assert res['Dharmacon_thermo'] is not None
        assert res['bp'] is not None
        assert res['cdf'] is not None
        assert res['fe'] is not None
        assert res['loopr'] is not None
        assert res['oligo'] is not None
        assert res['oligo_f'] is not None
        assert res['pdf'] is not None
        assert res['sample'] is not None
        assert res['sample_1000'] is not None
        assert res['sirna'] is not None
        assert res['sirna_f'] is not None
        assert res['sirna_s'] is not None
        assert res['smfe'] is not None
        assert res['sstrand'] is not None
        assert res['stability'] is not None

        res.cleanUp()

sfold_input = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']

if __name__ == '__main__':
    main()

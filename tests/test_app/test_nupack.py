#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.nupack     import Nupack

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

class NupackTest(TestCase):
    """Tests for Nupack application controller"""

    def setUp(self):
        self.input = nupack_input
                
    def test_input_as_lines(self):
        """Test nupack input as lines"""

        n = Nupack(InputHandler='_input_as_lines')
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in nupack_stdout])

        res = n(self.input)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_input_as_string(self):
        """Test nupack input as string"""

        n = Nupack()
        exp = '%s\n' % '\n'.join([str(i).strip('\n') for i in nupack_stdout])

        f = open('/tmp/single.fasta','w')
        txt = '\n'.join([str(i).strip('\n') for i in self.input])
        f.write(txt)
        f.close()
        res = n('/tmp/single.fasta')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/single.fasta')

    def test_get_result_path(self):
        """Tests nupack result path"""

        n = Nupack(InputHandler='_input_as_lines')
        res = n(self.input)
        self.assertEqualItems(res.keys(),['StdOut','StdErr','ExitStatus',
        'pair','ene'])
        assert res['pair'] is not None
        assert res['ene'] is not None
        res.cleanUp()

nupack_input = ['>seq1\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'\n']

nupack_stdout = ['****************************************************************\n', 
'NUPACK 1.2\n', 'Copyright 2003, 2004 by Robert M. Dirks & Niles A. Pierce\n', 
'California Institute of Technology\n', 'Pasadena, CA 91125 USA\n', '\n', 
'Last Modified: 03/18/2004\n', 
'****************************************************************\n', 
'\n', '\n', 'Fold.out Version 1.2: Complexity O(N^5) (pseudoknots enabled)\n', 
'Reading Input File...\n', 'Sequence Read.\n', 'Energy Parameters Loaded\n', 
'SeqLength = 71\n', 'Sequence and a Minimum Energy Structure:\n', 
'GGCUAGAUAGCUCAGAUGGUAGAGCAGAGGAUUGAAGAUCCUUGUGUCGUCGGUUCGAUCCCGGCUCUGGC\n', 
'.....((..((.((...{.{{{{{{.{.{{{{{{{{......)).))..))..}}}}}}}}}.}}}}}}.}\n', 
'mfe = -23.30 kcal/mol\n', 'pseudoknotted!\n']

if __name__ == '__main__':
    main()

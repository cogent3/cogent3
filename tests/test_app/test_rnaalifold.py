#!/usr/bin/env python

from os                    import remove
from cogent.util.unit_test import TestCase, main
from cogent.app.rnaalifold import RNAalifold, rnaalifold_from_alignment
from cogent.core.alignment import DataError

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Shandy Wikman","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RnaalifoldTest(TestCase):
    """Tests for Rnaalifold application controller"""

    def setUp(self):
        self.input = RNAALIFOLD_INPUT
        self.unaligned = \
            {'seq_0': 'GGUAGGUCGCUGGACUUGUCUCCUUGACUGUCCGGAAGGAGCGGU',
             'seq_1': 'GGUAGGUCGCUGGAUUGAUAUGAGUAUUGUCCGGAAGGAGCGGA',
             'seq_2': 'GGUAGGACGCGGGACUUCUGUUCAGGACUGUCCCGAAGGUGCGGU',
             'seq_3': 'GGUAGGUCGCCGCACGUCGCUUCAGGACUGUGCGGAAGGAGCGGU',
             'seq_4': 'GGUAGGUCGCUGUACUUCUAUCAGGACUGUACGGAAGGAGCGGU',
             }
        self.alignment = \
            {'seq_0': 'GGUAGGUCGCUGGAC-UUGUCUCCUUGACU-GUCCGGAAGGAGCGGU',
             'seq_1': 'GGUAGGUCGCUGGAU-UGAUAUGAGUAUU--GUCCGGAAGGAGCGGA',
             'seq_2': 'GGUAGGACGCGGGAC-UUCUGUUCAGGACU-GUCCCGAAGGUGCGGU',
             'seq_3': 'GGUAGGUCGCCGCACGUCGCUUCAGGAC--UGUGCGGAAGGAGCGGU',
             'seq_4': 'GGUAGGUCGCUGUAC-UUCUAUCAGGACU--GUACGGAAGGAGCGGU',
             }
        self.old_struct = '....(.((.((((((..(((....)))....))))))...)).)...'  
        #from version 1.4
        self.new_struct = '....(.((.((((((................))))))...)).)...'  
        #from version 1.8 - no easy way to check version so putting both in
        
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
    
    def test_rnaalifold_from_alignment_unaligned(self):
        """rnaalifold_from_alignment should handle unaligned seqs.
        """
        self.assertRaises(DataError,rnaalifold_from_alignment,self.unaligned)
    
    def test_rnaalifold_from_alignment(self):
        """rnaalifold_from_alignment should give correct result.
        """
        [[seq, struct,energy]] = rnaalifold_from_alignment(aln=self.alignment)
        try:
            self.assertEqual(struct,self.old_struct)
        except AssertionError:
            self.assertEqual(struct,self.new_struct)
        

RNAALIFOLD_INPUT = ['CLUSTAL\n', '\n', 'seq1 GGCTAGATAGCTCAGATGGT-AGAGCAGAGGATTGAAGATCCTTGTGTCGTCGGTTCGATCCCGGCTCTGGCC----\n']

if __name__ == '__main__':
    main()

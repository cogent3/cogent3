#!/usr/bin/env python

import os
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.struct.selection import einput
    from cogent.parse.pdb import PDBParser
    from cogent.parse.stride import stride_parser
    from cogent.app.stride import Stride, stride_xtra
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.struct.selection import einput
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.parse.stride import stride_parser
    from zenpdb.cogent.app.stride import Stride, stride_xtra

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class StrideParseTest(TestCase):
    """Tests for Stride application controller."""

    def setUp(self):
        input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(input_file))
        stride_app = Stride()
        res = stride_app(self.input_structure)
        self.lines = res['StdOut'].readlines()

    def test_stride_parser(self):
        """tests if output is parsed fully"""
        id_xtra = stride_parser(self.lines)

        assert len(id_xtra) < len(self.input_structure[(0,)][('A',)]) + \
                              len(self.input_structure[(0,)][('B',)])

        self.input_structure[(0,)][('A',)].remove_hetero()
        self.input_structure[(0,)][('B',)].remove_hetero()

        assert len(id_xtra) == len(self.input_structure[(0,)][('A',)]) + \
               len(self.input_structure[(0,)][('B',)])

    def test_stride_xtra(self):
        """tests if residues get annotated with parsed data."""
        stride_xtra(self.input_structure)
        self.assertEquals(\
            self.input_structure[(0,)][('A',)][(('H_HOH', 138, ' '),)].xtra, {})
        self.assertAlmostEquals(\
            self.input_structure[(0,)][('A',)][(('ILE', 86, ' '),)].xtra['STRIDE_ASA'], 13.9)
        self.input_structure[(0,)][('A',)].remove_hetero()
        self.input_structure[(0,)][('B',)].remove_hetero()
        all_residues = einput(self.input_structure, 'R')
        a = all_residues.data_children('STRIDE_ASA', xtra=True, forgiving=False)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.parse.pdb import PDBParser
    from cogent.struct.annotation import xtradata
    from cogent.struct.selection import einput
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.struct.annotation import xtradata
    from zenpdb.cogent.struct.selection import einput


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class AnnotationTest(TestCase):
    """tests if annotation get into xtra."""

    def setUp(self):

        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))

    def test_xtradata(self):
        """tests if an full_id's in the data dict are correctly parsed."""

        structure = einput(self.input_structure, 'S')[('2E12',)]
        model = einput(self.input_structure, 'M')[('2E12', 0)]
        chain = einput(self.input_structure, 'C')[('2E12', 0, 'B')]
        residue = einput(self.input_structure, 'R')[('2E12', 0, 'B', ('LEU', 24, ' '))]
        atom = einput(self.input_structure, 'A')[('2E12', 0, 'B', ('LEU', 24, ' '), ('CD1', ' '))]

        data_model = {(None, 0):{'model':1}}
        xtradata(data_model, structure)
        self.assertEquals(model.xtra, {'model': 1})

        data_chain = {(None, None, 'B'):{'chain':1}}
        xtradata(data_chain, model)
        self.assertEquals(chain.xtra, {'chain': 1})

        data_chain = {(None, 0, 'B'):{'chain': 2}}
        xtradata(data_chain, structure)
        self.assertEquals(chain.xtra['chain'], 2)

        data_residue = {(None, None, 'B', ('LEU', 24, ' ')):{'residue':1}}
        xtradata(data_residue, model)
        self.assertEquals(residue.xtra, {'residue': 1})

        data_residue = {(None, 0, 'B', ('LEU', 24, ' ')):{'residue':2}}
        xtradata(data_residue, structure)
        self.assertEquals(residue.xtra, {'residue': 2})


if __name__ == '__main__':
    main()


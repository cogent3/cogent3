#!/usr/bin/env python

import os
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.parse.pdb import PDBParser
    from cogent.struct.selection import einput, select
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.struct.selection import einput, select


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class AnnotationTest(TestCase):
    """tests selecting entities"""

    def setUp(self):
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))

    def test_einput(self):
        """tests einput."""
        structures = einput(self.input_structure, 'S')
        models = einput(self.input_structure, 'M')
        chains = einput(self.input_structure, 'C')
        residues = einput(self.input_structure, 'R')
        atoms = einput(self.input_structure, 'A')
        self.assertEquals(structures.level, 'H')
        self.assertEquals(models.level, 'S')
        self.assertEquals(chains.level, 'M')
        self.assertEquals(residues.level, 'C')
        self.assertEquals(atoms.level, 'R')
        atoms2 = einput(models, 'A')
        self.assertEquals(atoms, atoms2)
        atoms3 = einput(chains, 'A')
        self.assertEquals(atoms, atoms3)
        structures2 = einput(atoms, 'S')
        self.assertEquals(self.input_structure, structures2.values()[0])
        residues2 = einput(atoms, 'R')
        self.assertEquals(residues, residues2)

    def test_select(self):
        """tests select."""
        water = select(self.input_structure, 'R', 'H_HOH', 'eq', 'name')
        for residue in  water:
            self.assertTrue(residue.name == 'H_HOH')
        non_water = select(self.input_structure, 'R', 'H_HOH', 'ne', 'name')
        for residue in  non_water:
            self.assertTrue(residue.name != 'H_HOH')


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os
import numpy as np
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.parse.pdb import PDBParser
    from cogent.struct.selection import einput
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.struct.selection import einput


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class asaTest(TestCase):
    """Tests for surface calculations."""

    def setUp(self):
        self.arr = np.random.random(3000).reshape((1000, 3))
        self.point = np.random.random(3)
        self.center = np.array([0.5, 0.5, 0.5])

    def test_0import(self):
        # sort by name
        """tests if can import _contact cython extension."""
        global _contact
        from cogent.struct import _contact
        assert 'cnt_loop' in dir(_contact)

    def test_1import(self):
        # sort by name
        """tests if can import contact."""
        global contact
        from cogent.struct import contact

    def test_chains(self):
        """compares contacts diff chains"""
        self.input_file = os.path.join('data', '1A1X.pdb') # one chain
        self.input_structure = PDBParser(open(self.input_file))
        res = contact.contacts_xtra(self.input_structure)
        self.assertTrue(res == {})
        self.input_file = os.path.join('data', '2E12.pdb') # one chain
        self.input_structure = PDBParser(open(self.input_file))
        res = contact.contacts_xtra(self.input_structure)
        self.assertTrue(res)
        self.assertFloatEqual(\
        res[('2E12', 0, 'B', ('THR', 17, ' '), ('OG1', ' '))]['CONTACTS']\
        [('2E12', 0, 'A', ('ALA', 16, ' '), ('CB', ' '))][0], 5.7914192561064004)


    def test_symmetry(self):
        """compares contacts diff symmetry mates"""
        self.input_file = os.path.join('data', '2E12.pdb') # one chain
        self.input_structure = PDBParser(open(self.input_file))
        res = contact.contacts_xtra(self.input_structure, \
                                    symmetry_mode='uc',
                                    contact_mode='diff_sym')
        self.assertTrue(res)
        self.assertFloatEqual(\
        res[('2E12', 0, 'B', ('GLU', 77, ' '), ('OE2', ' '))]['CONTACTS']\
           [('2E12', 0, 'B', ('GLU', 57, ' '), ('OE2', ' '))][0], \
           5.2156557833123873)

    def test_crystal(self):
        """"compares contacts diff unit-cell-mates"""
        pass


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os, tempfile
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.parse.pdb import PDBParser
    from cogent.format.pdb import PDBWriter
    from cogent.struct.selection import einput
    from cogent.struct.manipulation import copy, clean_ical, \
    expand_symmetry, expand_crystal

except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from cogent.format.pdb import PDBWriter
    from zenpdb.cogent.struct.selection import einput
    from zenpdb.cogent.struct.manipulation import copy, clean_ical, \
    exapnd_symmetry, expand_crystal


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


class ManipulationTest(TestCase):
    """tests manipulationg entities"""

    def setUp(self):
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))

    def test_clean_ical(self):
        """tests the clean ical function which cleans structures."""
        chainB = self.input_structure.table['C'][('2E12', 0, 'B')]
        leu25 = self.input_structure.table['R'][('2E12', 0, 'B', \
                                                ('LEU', 25, ' '))]
        leu25icA = copy(leu25)
        self.assertTrue(leu25icA.parent is None)
        self.assertTrue(leu25icA is not leu25)
        self.assertTrue(leu25icA[(('N', ' '),)] is not leu25[(('N', ' '),)])
        leu25icA.setIc('A')
        self.assertEquals(leu25icA.getId(), (('LEU', 25, 'A'),))
        chainB.addChild(leu25icA)
        self.assertFalse(chainB[(('LEU', 25, 'A'),)] is \
                         chainB[(('LEU', 25, ' '),)])
        self.assertEquals(clean_ical(self.input_structure), \
                          ([], [('2E12', 0, 'B', ('LEU', 25, 'A'))]))
        clean_ical(self.input_structure, pretend=False)

        self.assertTrue(chainB[(('LEU', 25, 'A'),)] is leu25icA)
        self.assertFalse((('LEU', 25, 'A'),) in chainB.keys())
        self.assertFalse((('LEU', 25, 'A'),) in chainB)
        self.assertTrue((('LEU', 25, 'A'),) in chainB.keys(unmask=True))

        self.input_structure.setUnmasked(force=True)
        self.assertEquals(clean_ical(self.input_structure), \
                          ([], [('2E12', 0, 'B', ('LEU', 25, 'A'))]))
        clean_ical(self.input_structure, pretend=False, mask=False)
        self.assertFalse((('LEU', 25, 'A'),) in chainB.keys())
        self.assertFalse((('LEU', 25, 'A'),) in chainB)
        self.assertFalse((('LEU', 25, 'A'),) in chainB.keys(unmask=True))

    def test_0expand_symmetry(self):
        """tests the expansion of a asu to a unit-cell."""
        global fn
        mh = expand_symmetry(self.input_structure[(0,)])
        fd, fn = tempfile.mkstemp('.pdb')
        os.close(fd)
        fh = open(fn, 'w')
        PDBWriter(fh, mh, self.input_structure.raw_header)
        fh.close()

    def test_1expand_crystal(self):
        """tests the expansion of a unit-cell to a crystal"""
        fh = open(fn, 'r')
        input_structure = PDBParser(fh)
        self.assertTrue(input_structure.values(), 4) # 4 models
        sh = expand_crystal(input_structure)
        self.assertTrue(len(sh) == 27)
        fd, fn2 = tempfile.mkstemp('.pdb')
        os.close(fd)
        fh = open(fn2, 'w')
        a1 = einput(input_structure, 'A')
        a2 = einput(sh.values()[3], 'A')
        k = a1.values()[99].getFull_id()
        name = sh.values()[3].name
        a1c = a1[k].coords
        a2c = a2[(name,) + k[1:]].coords
        self.assertTrue(len(a1), len(a2))
        self.assertRaises(AssertionError, self.assertFloatEqual, a1c, a2c)
        PDBWriter(fh, sh)
        fh.close()
        os.unlink(fn)
        os.unlink(fn2)

if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os, tempfile
from unittest import main
from cogent.util.unit_test import TestCase

from cogent import FromFilenameStructureParser
from cogent.struct.selection import einput
from cogent.format import xyzrn

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class XyzrnTest(TestCase):
    """Tests conversion of PDB files into the informal xyzrn format."""

    def setUp(self):
        self.structure = FromFilenameStructureParser('data/1A1X.pdb')
        self.residues = einput(self.structure, 'R')
        self.atoms = einput(self.structure, 'A')
        self.residue8 = self.residues.values()[8]
        self.atom17 = self.atoms.values()[17]
        self.atom23 = self.atoms.values()[23]
        
    def test_write_atom(self):
        fd, fn = tempfile.mkstemp()
        os.close(fd)
        handle = open(fn, 'wb')
        xyzrn.XYZRNWriter(handle, [self.atom17])
        handle.close()
        handle = open(fn, 'rb')
        coords_radius = [float(n) for n in handle.read().split()[:4]]
        self.atom17.setRadius()
        radius = self.atom17.getRadius()
        self.assertFloatEqualRel(self.atom17.coords, coords_radius[:3])
        self.assertFloatEqualRel(radius, coords_radius[3])
        
if __name__ == '__main__':
    main()

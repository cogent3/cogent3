#!/usr/bin/env python

import os, tempfile
from cogent.util.unit_test import TestCase, main
from cogent.parse.pdb import PDBParser
from cogent.app.msms import Msms, surface_xtra
from cogent.struct.selection import einput

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

class MsmsTest(TestCase):
    """Tests for Msms application controller"""

    def setUp(self):

        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))

    def test_stdout_input_from_entity(self):
        """Test Stride when input is an entity"""

        s = Msms()
        res = s(self.input_structure)
        self.assertEqual(res['ExitStatus'], 0)
        stdout = res['StdOut'].read()
        assert stdout.find('1634 spheres 0 collision only, radii  1.600 to  1.850') != -1
        assert not res['StdOut'].read()
        assert res['ExitStatus'] == 0
        assert list(sorted(res.keys())) == list(sorted(['FaceFile', 'StdOut', 'AreaFile', 'StdErr', 'VertFile', 'ExitStatus']))
        af = res['AreaFile'].readlines()
        assert len(af) == 1635, len(af)
        ff = res['FaceFile'].readlines()
        assert ff[1].strip() == "#faces  #sphere density probe_r"
        assert ff[2].strip() == "25310    1634  1.00  1.50" or \
               ff[2].strip() == "51712    1634  1.00  1.50"
        assert len(ff) == 25313 or len(ff) == 51715
        vf = res['VertFile'].readlines()
        assert vf[1] == '#vertex #sphere density probe_r\n'
        assert vf[2].strip() == '12657    1634  1.00  1.50' or \
               vf[2].strip() == '25858    1634  1.00  1.50'
        
        assert len(vf) == 12660 or len(vf) == 25861
        res.cleanUp()
        
    def test_surface_xtra(self):
        res = surface_xtra(self.input_structure)
        assert res.shape == (12658, 3) or res.shape == (25859, 3)
        assert res is self.input_structure.xtra['MSMS_SURFACE']
        chains = einput(self.input_structure, 'C')
        chainA, chainB = chains.sortedvalues()
        resA = surface_xtra(chainA)
        assert len(resA) == 6223 or len(resA) == 12965
        resB = surface_xtra(chainB)
        assert len(resB) == 6620 or len(resB) == 13390
        assert chainB.xtra['MSMS_SURFACE'] is resB is \
        self.input_structure[(0,)][('B',)].xtra['MSMS_SURFACE']
        
if __name__ == '__main__':
    main()

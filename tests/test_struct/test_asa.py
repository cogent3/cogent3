#!/usr/bin/env python

import os
import numpy as np
from numpy import sum
try:
    from cogent.util.unit_test import TestCase, main
    from cogent.app.util import ApplicationNotFoundError
    from cogent.parse.pdb import PDBParser
    from cogent.struct.selection import einput
    from cogent.maths.stats.test import correlation
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main
    from zenpdb.cogent.parse.pdb import PDBParser
    from zenpdb.cogent.struct.selection import einput
    from zenpdb.maths.stats.test import correlation


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
        """tests if can import _asa cython extension."""
        global _asa
        from cogent.struct import _asa
        assert 'asa_loop' in dir(_asa)

    def test_1import(self):
        # sort by name
        """tests if can import asa."""
        global asa
        from cogent.struct import asa

    def test_asa_loop(self):
        """tests if inner asa_loop (cython) performs correctly"""
        self.lcoords = np.array([[-4., 0, 0], [0, 0, 0], [4, 0, 0], [10, 0, 0]])
        self.qcoords = np.array([[0., 0, 0], [4., 0, 0]])
        self.lradii = np.array([2., 3.])
        self.qradii = np.array([3., 2.])
        #spoints, np.ndarray[DTYPE_t, ndim =1] box,\
        # DTYPE_t probe, unsigned int bucket_size, MAXSYM =200000)
        self.spoints = np.array([[1., 0., 0.], [-1., 0., 0.], [0., 1., 0.], \
                                 [0., -1., 0.], [0., 0., 1.], [0., 0., -1.]])
        output = _asa.asa_loop(self.qcoords, self.lcoords, self.qradii, \
                               self.lradii, self.spoints, \
                    np.array([-100., -100., -100., 100., 100., 100.]), 1., 10)
        self.assertFloatEqual(output, np.array([ 75.39822369, 41.88790205]))

    def test_asa_xtra(self):
        """test internal asa"""
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        self.assertRaises(ValueError, asa.asa_xtra, self.input_structure, mode='a')
        result = asa.asa_xtra(self.input_structure)
        a = einput(self.input_structure, 'A')
        for i in range(len(result)):
            self.assertEquals(result.values()[i]['ASA'], a[result.keys()[i]].xtra['ASA'])
        r = einput(self.input_structure, 'R')
        for water in  r.selectChildren('H_HOH', 'eq', 'name').values():
            self.assertFalse('ASA' in water.xtra)
        for residue in  r.selectChildren('H_HOH', 'ne', 'name').values():
            for a in residue:
                self.assertTrue('ASA' in a.xtra)
        result = asa.asa_xtra(self.input_structure, xtra_key='SASA')
        for residue in  r.selectChildren('H_HOH', 'ne', 'name').values():
            for a in residue:
                a.xtra['ASA'] == a.xtra['SASA']

    def test_asa_xtra_stride(self):
        """test asa via stride"""
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        try:
            result = asa.asa_xtra(self.input_structure, 'stride')
        except ApplicationNotFoundError: 
            return
        self.assertAlmostEqual(self.input_structure[(0,)][('B',)]\
                               [(('LEU', 35, ' '),)].xtra['STRIDE_ASA'], 17.20)

    def test_compare(self):
        """compares internal asa to stride."""
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        try:
            asa.asa_xtra(self.input_structure, mode='stride')
        except ApplicationNotFoundError: 
            return            
        asa.asa_xtra(self.input_structure)
        self.input_structure.propagateData(sum, 'A', 'ASA', xtra=True)
        residues = einput(self.input_structure, 'R')
        asa1 = []
        asa2 = []
        for residue in  residues.selectChildren('H_HOH', 'ne', 'name').values():
            asa1.append(residue.xtra['ASA'])
            asa2.append(residue.xtra['STRIDE_ASA'])
        self.assertAlmostEqual(correlation(asa1, asa2)[1], 0.)

    def test_uc(self):
        """compares asa within unit cell."""
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        asa.asa_xtra(self.input_structure, symmetry_mode='uc', xtra_key='ASA_UC')
        asa.asa_xtra(self.input_structure)
        self.input_structure.propagateData(sum, 'A', 'ASA', xtra=True)
        self.input_structure.propagateData(sum, 'A', 'ASA_UC', xtra=True)
        residues = einput(self.input_structure, 'R')
        x = residues[('2E12', 0, 'B', ('GLU', 77, ' '))].xtra.values()
        self.assertTrue(x[0] != x[1])

    def test_uc2(self):
        self.input_file = os.path.join('data', '1LJO.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        asa.asa_xtra(self.input_structure, symmetry_mode='uc', xtra_key='ASA_XTAL')
        asa.asa_xtra(self.input_structure)
        self.input_structure.propagateData(sum, 'A', 'ASA', xtra=True)
        self.input_structure.propagateData(sum, 'A', 'ASA_XTAL', xtra=True)
        residues = einput(self.input_structure, 'R')
        r1 = residues[('1LJO', 0, 'A', ('ARG', 65, ' '))]
        r2 = residues[('1LJO', 0, 'A', ('ASN', 46, ' '))]
        self.assertFloatEqual(r1.xtra.values(),
                              [128.94081270529105, 22.807700865674093])
        self.assertFloatEqual(r2.xtra.values(),
                              [115.35738419425566, 115.35738419425566])

    def test_crystal(self):
        """compares asa within unit cell."""
        self.input_file = os.path.join('data', '2E12.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        asa.asa_xtra(self.input_structure, symmetry_mode='uc', crystal_mode=2, xtra_key='ASA_XTAL')
        asa.asa_xtra(self.input_structure)
        self.input_structure.propagateData(sum, 'A', 'ASA', xtra=True)
        self.input_structure.propagateData(sum, 'A', 'ASA_XTAL', xtra=True)
        residues = einput(self.input_structure, 'R')
        r1 = residues[('2E12', 0, 'A', ('ALA', 42, ' '))]
        r2 = residues[('2E12', 0, 'A', ('VAL', 8, ' '))]
        r3 = residues[('2E12', 0, 'A', ('LEU', 25, ' '))]
        self.assertFloatEqual(r1.xtra.values(), \
                                [32.041070749038823, 32.041070749038823])
        self.assertFloatEqual(r3.xtra.values(), \
                               [0., 0.])
        self.assertFloatEqual(r2.xtra.values(), \
                                [28.873559956056916, 0.0])

    def test_bio(self):
        """compares asa within a bio unit."""
        self.input_file = os.path.join('data', '1A1X.pdb')
        self.input_structure = PDBParser(open(self.input_file))
        asa.asa_xtra(self.input_structure, symmetry_mode='bio', xtra_key='ASA_BIO')
        asa.asa_xtra(self.input_structure)
        self.input_structure.propagateData(sum, 'A', 'ASA', xtra=True)
        self.input_structure.propagateData(sum, 'A', 'ASA_BIO', xtra=True)
        residues = einput(self.input_structure, 'R')
        r1 = residues[('1A1X', 0, 'A', ('GLU', 37, ' '))]
        r2 = residues[('1A1X', 0, 'A', ('TRP', 15, ' '))]
        self.assertFloatEqual(r1.xtra.values(), \
                                [20.583191467544726, 78.996394472066541])
        self.assertFloatEqual(r2.xtra.values(), \
                                [136.41436710386989, 136.41436710386989])





if __name__ == '__main__':
    main()

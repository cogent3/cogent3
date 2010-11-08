#!/usr/bin/env python
"""Unit tests for the pdb parser.
"""
from cogent.util.unit_test import TestCase, main
from cogent.core.entity import Structure
from cogent.parse.structure import FromFilenameStructureParser, FromFileStructureParser

__author__ = "Marcin Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Production"


class structuresTests(TestCase):
    """Tests of cogent.parse.structure UI functions."""
    
    def test_FromFilenameStructureParser(self):
        structure = FromFilenameStructureParser('data/1LJO.pdb', 'pdb')
        self.assertRaises(TypeError, FromFilenameStructureParser, open('data/1LJO.pdb'), 'pdb')
        assert isinstance(structure, Structure)
    
    def test_FromFileStructureParser(self):
        structure = FromFileStructureParser(open('data/1LJO.pdb'), 'pdb')
        assert isinstance(structure, Structure)
        self.assertRaises(TypeError, FromFileStructureParser, 'data/1LJO.pdb', 'pdb')
        assert isinstance(structure, Structure)         
    
    
if __name__ == '__main__':
    main()

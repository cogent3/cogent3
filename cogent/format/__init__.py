#!/usr/bin/env python
"""Package cogent.format: provides modules for writing specific file formats.

Currently provides:
    mage: writers for the MAGE 3D visualization program
    xml:  xml base class
    file: general functions to read and write files
"""
__all__ = ['alignment', 'clustal', 'fasta', 'mage', 'motif', 'nexus', 
           'pdb_color', 'pdb', 'phylip', 'rna_struct', 'stockholm', 'structure', 
           'table', 'text_tree', 'xyzrn']

__author__ = ""
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Gavin Huttley", "Matthew Wakefield",
               "Rob Knight", "Sandra Smit", "Peter Maxwell", "Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

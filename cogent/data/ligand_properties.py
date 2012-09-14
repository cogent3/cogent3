#!/usr/bin/env python

"""Properties of ligands are data from 'tables' about ligands and their 
atoms."""

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Production"

HOH_NAMES = ['H_HOH', 'H_WAT', 'H_DOH', 'H_HOD', 'H_DOD']
WATER_NAMES = HOH_NAMES

LIGAND_ATOM_PROPERTIES = {
    ('H_HOH', ' O  '): [1.60]
    }

LIGAND_AREAIMOL_VDW_RADII = dict([(k, v[0]) for k, v in LIGAND_ATOM_PROPERTIES.iteritems()])

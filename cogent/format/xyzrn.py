#!/usr/bin/env python

"""Function for XYZRN (coordinates followed by radius and id_hash) format 
output. This is a rather free implementation of the no-standard, but is 
compatible with MSMS."""

from itertools import chain
from cogent.struct.selection import einput
from cogent.data.protein_properties import AREAIMOL_VDW_RADII
from cogent.data.ligand_properties import LIGAND_AREAIMOL_VDW_RADII

XYZRN_COORDS_STRING = "%8.3f %8.3f %8.3f %8.3f %d %s\n"
AREAIMOL_VDW_RADII.update(LIGAND_AREAIMOL_VDW_RADII)

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

def write_header(header):
    """Write the header. Not implemented."""
    return (header or [])

def write_coords(atoms, radius_type):
    """Write coordinate lines from a list of ``Atom`` instances. Takes a string
    which identifies the radius type."""
    lines = []
    radius_type = eval(radius_type)
    for atom in atoms.sortedvalues():
        residue_name = atom.parent.name
        atom_name = atom.name
        radius = radius_type[(residue_name, atom_name)]
        if radius == 0.0:
            continue
        (x, y, z) = atom.coords
        args = (x, y, z, radius, 1, hash((x, y, z)))
        line = XYZRN_COORDS_STRING % args
        lines.append(line)
    return lines

def write_trailer(trailer):
    """Write the trailer. Not implemented."""
    return (trailer or [])

def XYZRNWriter(f, entities, radius_type=None, header=None, trailer=None):
    """Function which writes XYZRN files from ``Entity`` instances."""
    radius_type = (radius_type or 'AREAIMOL_VDW_RADII')

    structure = einput(entities, level='A', name='structure')
    header = write_header(header)
    coords = write_coords(structure, radius_type)
    trailer = write_trailer(trailer)

    for part in chain(header, coords, trailer):
        f.writelines(part)

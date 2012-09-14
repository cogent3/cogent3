"""Classes and functions for computing and manipulating accessible 
surface areas (ASA)."""

from cogent.app.stride import Stride
from cogent.parse.stride import stride_parser
from cogent.struct.selection import einput
from cogent.struct.annotation import xtradata
from cogent.maths.geometry import sphere_points, coords_to_symmetry, \
                                  coords_to_crystal
from _asa import asa_loop
from numpy import array, r_


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


def _run_asa(atoms, lattice_coords, spoints, probe=1.4, bucket_size=5, \
             MAXSYM=200000):
    """Runs an ASA calculation. This function takes a selection of atoms 
    (in the most common case all atoms in a structure) and lattice coordinates 
    (in the most common a 3x3 box of unit-cells).
    
    This function should be considered low-level and not part of the interface.
    
    Arguments:
    
        - atoms: Holder of atom entities.
        - lattice_coords: Numpy array of coordinates.
        - spoints: an array of coordinates on the unit sphere, defines the 
          accuracy of ASA calculation.
        - probe: size of the probe i.e. solvent molecule.
        - bucket_size: see: ``KDTree``.
        - MAXSYM (int): maximum number of symmetry generated atoms.
    """
    # get array of radii inflated by probe size of the selection of atoms.
    atom_radii = array(atoms.getData('radius', forgiving=False)) + probe
    # get array of coordinates
    atom_coords = array(atoms.getData('coords', forgiving=False))
    # calculate bounding box in a form of an array
    search_limit = 2 * (2.0 + probe)    # 2.0 is maximum atom radius
    atom_box = r_[atom_coords.min(axis=0) - search_limit, \
                  atom_coords.max(axis=0) + search_limit]
    # the lattice coordinates are atom coordinates after transformations in a 
    # 4D-array this array gets reshaped into an all_atoms x 3 array.
    shape = lattice_coords.shape
    lattice_coords = \
        lattice_coords.reshape((shape[0] * shape[1] * shape[2], shape[3]))
    # this calls the cython code which loops over all query atoms, surface 
    # points, and lattice atoms
    return asa_loop(atom_coords, lattice_coords, atom_radii, atom_radii, \
                    spoints, atom_box, probe, bucket_size, MAXSYM)

def _prepare_entities(entities):
    """Prepares input entities for ASA calculation, which includes masking water
    molecules and water chains.
    """
    # First we mask all water residues and chains with all residues masked 
    # (water chains).
    lattice_residues = einput(entities, 'R')
    lattice_residues.maskChildren('H_HOH', 'eq', 'name')
    lattice_chains = einput(entities, 'C')
    lattice_chains.maskChildren([], 'eq', 'values', method=True)
    # if no residues or chains are left - no atoms to work with, 
    # abort with warning.
    if not lattice_chains.values():
        # the following makes sure that masking changes by the above
        # tests are reverted.
        lattice_structures = einput(entities, 'S')
        lattice_structures.setUnmasked(force=True)
        raise ValueError('No unmasked atoms to build lattice.')
    # these are all atoms we can work with
    lattice_atoms = einput(entities, 'A')
    lattice_atoms.dispatch('setRadius')

def _postpare_entities(entities):
    """Restores entities after ASA calculation, which includes unmasking."""
    structures = einput(entities, 'S')
    structures.setUnmasked(force=True)


def _prepare_asa(entities, symmetry_mode=None, crystal_mode=None, points=960, \
                **kwargs):
    """Prepares the atomic solvent-accessible surface area (ASA) calculation.
    
    Arguments:
    
        - entities: input entities for ASA calculation (most commondly a 
          structure entity).
        - symmetry_mode (str): One of 'uc', 'bio' or 'table'. This defines the 
          transformations of applied to the coordinates of the input entities. 
          It is one of 'bio', 'uc' or 'table'. Where 'bio' and 'uc' are 
          transformations to create the biological molecule or unit-cell from 
          the PDB header. The 'table' uses transformation matrices derived from 
          space-group information only using crystallographic tables(requires 
          ``cctbx``).
        - crystal_mode (int): Defines the number of unit-cells to expand the 
          initial unit-cell into. The number of unit-cells in each direction 
          i.e. 1 is makes a total of 27 unit cells: (-1, 0, 1) == 3, 3^3 == 27
        - points: number of points on atom spheres higher is slower but more 
          accurate.    

    Additional keyworded arguments are passed to the ``_run_asa`` function.
    """
    # generate uniform points on the unit-sphere
    spoints = sphere_points(points)
    # prepare entities for asa calculation
    # free-floating area mode
    result = {}
    atoms = einput(entities, 'A')
    if not symmetry_mode and not crystal_mode:

        coords = array(atoms.getData('coords', forgiving=False))
        coords = array([[coords]]) # fake 3D and 4D
        idx_to_id = dict(enumerate(atoms.getData('getFull_id', \
                                                forgiving=False, method=True)))
        asas = _run_asa(atoms, coords, spoints, **kwargs)
        for idx in xrange(asas.shape[0]):
            result[idx_to_id[idx]] = asas[idx]
    # crystal-contact area mode    
    elif symmetry_mode in ('table', 'uc'):
        structure = einput(entities, 'S').values()[0]
        sh = structure.header
        coords = array(atoms.getData('coords', forgiving=False))
        idx_to_id = dict(enumerate(atoms.getData('getFull_id', \
                                                forgiving=False, method=True)))
        # expand to unit-cell, real 3D
        coords = coords_to_symmetry(coords, \
                                            sh[symmetry_mode + '_fmx'], \
                                            sh[symmetry_mode + '_omx'], \
                                            sh[symmetry_mode + '_mxs'], \
                                            symmetry_mode)
        # expand to crystal, real 4D
        if crystal_mode:
            coords = coords_to_crystal(coords, \
                                               sh[symmetry_mode + '_fmx'], \
                                               sh[symmetry_mode + '_omx'], \
                                               crystal_mode) # real 4D
        else:
            coords = array([coords]) # fake 4D
        asas = _run_asa(atoms, coords, spoints, **kwargs)
        for idx in xrange(asas.shape[0]):
            result[idx_to_id[idx]] = asas[idx]

     # biological area mode
    elif symmetry_mode == 'bio':
        structure = einput(entities, 'S').values()[0]
        chains = einput(entities, 'C')
        sh = structure.header
        start = 0
        for chain_ids, mx_num in sh['bio_cmx']:
            sel = chains.selectChildren(chain_ids, 'contains', 'id').values()
            atoms = einput(sel, 'A')
            coords = array(atoms.getData('coords', forgiving=False))
            idx_to_id = dict(enumerate(atoms.getData('getFull_id', \
                                              forgiving=False, method=True)))
            stop = start + mx_num
            coords = coords_to_symmetry(coords, \
                                               sh['uc_fmx'], \
                                               sh['uc_omx'], \
                                               sh['bio_mxs'][start:stop], \
                                               symmetry_mode)
            coords = array([coords])
            start = stop
            asas = _run_asa(atoms, coords, spoints, **kwargs)
            for idx in xrange(asas.shape[0]):
                result[idx_to_id[idx]] = asas[idx]
    return result


def asa_xtra(entities, mode='internal', xtra_key=None, **asa_kwargs):
    """Calculates accessible surface areas (ASA) and puts the results into the
    xtra dictionaries of entities.
    
    Arguments:
    
        - entities: an entity or sequence of entities
        - mode(str): 'internal' for calculations using the built-in cython code 
          or 'stride' if the stride binary should be called to do the job.
        - xtra_key(str): Key in the xtra dictionary to hold the result for each
          entity
    
    Additional keyworded arguments are passed to the ``_prepare_asa`` and 
    ``_run_asa`` functions.
    """
    xtra_key = xtra_key or 'ASA'
    structures = einput(entities, 'S')
    if len(structures.values()) > 1:
        raise ValueError('Entities from multiple structures are not supported.')
    if mode == 'internal':
        _prepare_entities(entities) # mask waters
        result = _prepare_asa(entities, **asa_kwargs) # calculate ASA
        _postpare_entities(entities) # unmask waters
        result = dict([(id, {xtra_key:v}) for id, v in result.iteritems()])
        xtradata(result, structures)
    elif mode == 'stride':
        models = einput(entities, 'M')
        stride_app = Stride()
        result = stride_app(entities)['StdOut'].readlines()
        result = stride_parser(result)
        xtradata(result, structures.values()[0][(0,)])
    else:
        raise ValueError('Not a valid mode: "%s"' % mode)
    return result

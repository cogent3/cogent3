"""Classes and functions for computing and manipulating accessible 
surface areas (ASA)."""

from cogent.struct.selection import einput
from cogent.struct.annotation import xtradata
from cogent.maths.geometry import coords_to_symmetry, \
                                  coords_to_crystal
from _contact import cnt_loop
from collections import defaultdict
from numpy import array, r_, sqrt, int64


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


def _prepare_contacts(query, model=None, level='A', search_limit=6.0, \
                      contact_mode='diff_chain', symmetry_mode=None, \
                      crystal_mode=None, **kwargs):
    """Prepares distance contact calculations.
    
    Arguments:
    
        - query(entitie[s]): query entitie[s] for contact calculation 
          (most commonly a structure entity).
        - model(entity): a Model entity which will be transformed according to
          symmetry_mode and crystal_mode. (most commonly it is the same as the
          query)
        - level(str): The level in the hierarchy at which distances will be 
          calculated (most commonly 'A' for atoms) 
        - search_limit(float): maximum distance in Angstrom's
        - contact_mode(str): One of "diff_cell", "diff_sym", "diff_chain".
          Defines the allowed contacts i.e. requires that contacts are by 
          entities, which have: "diff_cell" different unit cells; "diff_sym" 
          different symmetry operators (if in the same unit cell) "diff_chain" 
          with different chain ids (if in the same unit cell and symmetry).
        - symmetry_mode (str): One of 'uc', 'bio' or 'table'. This defines the 
          transformations of applied to the coordinates of the input entities. 
          It is one of 'bio', 'uc' or 'table'. Where 'bio' and 'uc' are 
          transformations to create the biological molecule or unit-cell from 
          the PDB header. The 'table' uses transformation matrices derived from 
          space-group information only using crystallographic tables(requires 
          ``cctbx``).        
        - crystal_mode (int): Defines the number of unit-cells to expand the 
          initial unit-cell into. The number of unit cells in each direction 
          i.e. 1 is makes a total of 27 unit cells: (-1, 0, 1) == 3, 3^3 == 27
    
    Additional arguments are passed to the ``cnt_loop`` Cython function.
    """

    contact_mode = {'diff_asu'  :0,
                    'diff_sym' :1,
                    'diff_chain':2 }[contact_mode]

    # determine unique structure
    structure = einput(query, 'S').values()[0]
    sh = structure.header
    # if not specified otherwise the lattice is the first model
    lattice = model or structure[(0,)]
    lents = einput(lattice, level)
    lents_ids = lents.getData('getFull_id', forgiving=False, method=True)
    lcoords = array(lents.getData('coords', forgiving=False))
    qents = einput(query, level)
    qents_ids = qents.getData('getFull_id', forgiving=False, method=True)
    qcoords = array(qents.getData('coords', forgiving=False))

    if symmetry_mode:
        if symmetry_mode == 'table':
            lcoords = coords_to_symmetry(lcoords, \
                                         sh['table_fmx'], \
                                         sh['table_omx'], \
                                         sh['table_mxs'], \
                                         symmetry_mode)
        elif symmetry_mode == 'uc':
            lcoords = coords_to_symmetry(lcoords, \
                                         sh['uc_fmx'], \
                                         sh['uc_omx'], \
                                         sh['uc_mxs'], \
                                         symmetry_mode)
        elif symmetry_mode == 'bio':
            # TODO see asa
            raise ValueError("Unsupported symmetry_mode: %s" % symmetry_mode)
        else:
            raise ValueError("Unsupported symmetry_mode: %s" % symmetry_mode)
    else:
        lcoords = array([lcoords]) # fake 3D
    if crystal_mode:
        zero_tra = {1:13, 2:62, 3:171}[crystal_mode]
        # 0,0,0 translation is: Thickened cube numbers: 
        # a(n)=n*(n^2+(n-1)^2)+(n-1)*2*n*(n-1).
        # 1, 14, 63, 172, 365, 666, 1099, 1688, 2457, 3430, 4631, 6084, 7813 ...
        if symmetry_mode == 'table':
            lcoords = coords_to_crystal(lcoords, \
                                        sh['table_fmx'], \
                                        sh['table_omx'], \
                                        crystal_mode)
        elif symmetry_mode == 'uc':
            lcoords = coords_to_crystal(lcoords, \
                                        sh['uc_fmx'], \
                                        sh['uc_omx'], \
                                        crystal_mode)
        else:
            raise ValueError('crystal_mode not possible for "bio" symmetry')
    else:
        zero_tra = 0
        lcoords = array([lcoords]) # fake 4D
    shape = lcoords.shape
    lcoords = lcoords.reshape((shape[0] * shape[1] * shape[2], shape[3]))
    box = r_[qcoords.min(axis=0) - search_limit, \
             qcoords.max(axis=0) + search_limit]
    lc = [] # lattice chain
    qc = [] # query chain
    lchains = [i[2] for i in lents_ids]
    qchains = [i[2] for i in qents_ids]
    allchains = set()
    allchains.update(lchains)
    allchains.update(qchains)
    chain2id = dict(zip(allchains, range(len(allchains))))
    for lent_id in lents_ids:
        lc.append(chain2id[lent_id[2]])
    for qent_id in qents_ids:
        qc.append(chain2id[qent_id[2]])
    lc = array(lc, dtype=int64)
    qc = array(qc, dtype=int64)
    # here we leave python
    (idxc, n_src, n_asu, n_sym, n_tra, n_dst) = cnt_loop(\
                            qcoords, lcoords, qc, lc, shape[1], shape[2], \
                            zero_tra, contact_mode, search_limit, box, \
                            **kwargs)

    result = defaultdict(dict)
    for contact in xrange(idxc):
        qent_id = qents_ids[n_src[contact]]
        lent_id = lents_ids[n_asu[contact]]
        result[qent_id][lent_id] = (sqrt(n_dst[contact]), n_tra[contact], n_sym[contact])
    return result


def contacts_xtra(query, xtra_key=None, **cnt_kwargs):
    """Finds distance contacts between entities. This function searches for 
    contacts for query entities (query) either within the asymmetric unit,
    biological molecule, unit-cell or crystal.
    
    Arguments:
        
        - query (entitie[s]): query entity or sequence entities
        - xtra_key (str): name of the key
        
    Additional keyworded arguments are passed to the ``_prepare_contacts`` 
    functon.
    """
    xtra_key = xtra_key or 'CONTACTS'
    structures = einput(query, 'S')
    if len(structures.values()) > 1:
        raise ValueError('Entities from multiple structures are not supported.')
    result = _prepare_contacts(query, **cnt_kwargs) # calculate CONTACTS
    result = dict([(id, {xtra_key:v}) for id, v in result.iteritems()])
    xtradata(result, structures)
    return result

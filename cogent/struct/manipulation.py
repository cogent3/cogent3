"""Contains functions to manipulate i.e. modify macromolecular entities."""


from numpy import array
from itertools import izip
from cogent.core.entity import HIERARCHY, copy, StructureHolder, ModelHolder
from cogent.maths.geometry import coords_to_symmetry, coords_to_crystal
from cogent.struct.selection import einput


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


def clean_ical(entities, pretend=True, mask=True):
    """Removes or masks entities with ambiguous (i)nsertion (c)odes or 
    (a)lternate (l)ocations.
    
    Arguments:
    
        - entities: universal input see: ``cogent.struct.selection.einput``
        - pretend: If ``True`` only reports icals and does not mask or remove 
          anything.
        - mask (boolean): If pretend is ``False`` masks entities instead of 
          removing them.
   
    This function does not check for occupancy. I retains the residue which is 
    first when sorted by id number, insertion code and finally name. Residues 
    without IC come first. Atoms within a retained residue are sorted according 
    to PDB rules and the first one is chosen. If The first entity has an IC or 
    alt_loc different from ' ' it will be changed to ' '.
    """
    conflicts = []
    changes = []
    residues = einput(entities, 'R')
    id_r = [[None, None, None]]
    for r in residues.sortedvalues():  # sort by id, ic, name
        id_a = [[None, None]]
        if r.res_id == id_r[0][1]:      # on collision choose first ...
            conflicts.append(r.getFull_id())
            if not pretend:
                if mask:
                    r.setMasked(True)
                else:
                    r.parent.delChild(r.id)
            continue                    # an entity could be in other holders
        # keep it there as-is
        for a in r.sortedvalues():     # sort by id, alt_loc (' ', 'A' ...)
            if a.at_id == id_a[0][0]:   # on collision choose first
                conflicts.append(a.getFull_id())
                if not pretend:
                    if mask:
                        a.setMasked(True)
                    else:
                        r.delChild(a.id)
            else:
                if a.id[0][1] != ' ':
                    changes.append((a.getFull_id(), ((a.id[0][0], ' '),)))
                    if not pretend:
                        a.setAlt_loc(' ')
                        try:
                            a.parent.updateIds()
                        except AttributeError:
                            pass
                id_a = a.id
        if r.id[0][2] != ' ':
            changes.append((r.getFull_id(), ((r.id[0][0], r.id[0][1], ' '),)))
            if not pretend:
                r.set_res_ic(' ')
                try:
                    r.parent.updateIds()
                except AttributeError:
                    pass
        id_r = r.id
    return (changes, conflicts)

def expand_symmetry(model, mode='uc', name='UC', **kwargs):
    """Applies the symmetry operations defined by the header of the PDB files to
    the given ``Model`` entity instance. Returns a ``ModelHolder`` entity.
    
    Arguments:
    
        - model: model entity to expand
        - mode: 'uc', 'bio' or 'raw'
        - name: optional name of the ``ModelHolder`` instance.
    
    Requires a PDB file with a correct CRYST1 field and space group information.
    """
    structure = model.getParent('S')
    sh = structure.header
    fmx = sh['uc_fmx']
    omx = sh['uc_omx']
    mxs = sh['uc_mxs']
    # get initial coordinates
    atoms = einput(model, 'A')
    coords = array(atoms.getData('coords'))
    # expand the coordinates to symmetry
    all_coords = coords_to_symmetry(coords, fmx, omx, mxs, mode)
    models = ModelHolder(name)

    for i in xrange(0, len(mxs)):
        # copy model
        new_model = copy(model) # with additional models which
        new_atoms = einput(new_model, 'A')
        # patch with coordinates
        new_coords = all_coords[i]
        for (atom_id, new_coord) in izip(atoms.keys(), new_coords):
            new_atoms[atom_id[1:]].coords = new_coord
        # give it an id: the models are numbered by the symmetry operations with
        # identity being the first model
        new_model.setName(i)
        models.addChild(new_model)
    return models

def expand_crystal(structure, n=1, name='XTAL'):
    """Expands the contents of a structure to a crystal of a given size. Returns
    a `` StructureHolder`` entity instance. 

    Arguments:
    
        - structure: ``Structure`` entity instance. 
        - n: number number of unit-cell layers.
        - name: optional name.

    Requires a PDB file with correct CRYST1 field and space group information.
    """
    sh = structure.header
    sn = structure.name
    fmx = sh['uc_fmx']
    omx = sh['uc_omx']
    # get initial coorinates
    atoms = einput(structure, 'A')
    coords = array([atoms.getData('coords')]) # fake 3D
    # expand the coordinates to crystal
    all_coords = coords_to_crystal(coords, fmx, omx, n)
    structures = StructureHolder(name)
    rng = range(-n, n + 1) # a range like -2, -1, 0, 1, 2
    vectors = [(x, y, z) for x in rng for y in rng for z in rng]
    for i, (u, v, w) in enumerate(vectors):
        new_structure = copy(structure)
        new_atoms = einput(new_structure, 'A')
        new_coords = all_coords[i, 0]
        for (atom_id, new_coord) in izip(atoms.keys(), new_coords):
            new_atoms[atom_id].coords = new_coord
        new_structure.setName("%s_%s%s%s" % (sn, u, v, w))
        structures.addChild(new_structure)
    return structures


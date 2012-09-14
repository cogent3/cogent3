"""Functions to create PDB files from entities optionally with data in the
B-factor and Q (occupancy) columns.
"""

from collections import defaultdict
from itertools import chain

from cogent.struct.selection import einput
from cogent.data.protein_properties import AA_NAMES
from cogent.parse.pdb import dict2pdb, dict2ter

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

def write_header(header):
    if not isinstance(header, dict):
        return (header or [])
    xt = ('REMARK 200 EXPERIMENT TYPE : %s\n', 'experiment_type')
    sg = ('REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: %s\n', 'space_group')
    tests = [xt, sg]
    results = []
    for test in tests:
        string, key = test
        try:
            value = header[key]
            results.append(string % value)
        except (KeyError, TypeError):
            pass
    return results

def write_coords(atoms):
    old_fields = defaultdict(str)
    coords = ['MODEL        1\n']
    for atom in atoms.sortedvalues():
        fields = atom.getDict()
        if (old_fields['chain_id'] != fields['chain_id']) and old_fields['chain_id'] and old_fields['res_name'] in AA_NAMES:
                coords.append(dict2ter(old_fields))
        if (old_fields['model'] != fields['model']) and old_fields['model'] != '': # model can be 0 :)
                if old_fields['chain_id'] and old_fields['res_name'] in AA_NAMES:
                    coords.append(dict2ter(old_fields))
                coords.append('ENDMDL\n')
                coords.append('MODEL     %4i\n' % (fields['model'] + 1))
        coords.append(dict2pdb(fields))
        old_fields = fields
    if fields['res_name'] in AA_NAMES:
        coords.append(dict2ter(fields))
    coords.append('ENDMDL\n')
    coords.append('END   \n')
    return coords

def write_trailer(trailer):
    if not isinstance(trailer, dict):
        return (trailer or [])


def number(data, rest_val):
    try:
        return float(data)
    except TypeError:
        return rest_val

def iterable(data, rest_val):
    try:
        return float(len(data))
    except TypeError:
        return rest_val

def PDBWriter(f, entities, header_=None, trailer_=None):
        structure = einput(entities, level='A', name='atoms')
        # hierarchy: args, dicts
        try:
            header = (header_ or entities.raw_header or entities.header)
        except AttributeError:
            header = header_
        try:
            trailer = (trailer_ or entities.raw_trailer or entities.trailer)
        except AttributeError:
            trailer = trailer_
        coords = write_coords(structure)
        header = write_header(header)
        trailer = write_trailer(trailer)
        for part in chain(header, coords, trailer):
            f.writelines(part)
        # did not open do not close
        # f.close()

def PDBXWriter(f, entities, level, b_key, b_mode=None, b_val=0.0, q_key=None, \
               q_mode=None, q_val=0.0):
    """Writes data from the ``xtra`` dictionary into B and Q columns of a 
    PDB file. The level from which the dictionary is taken can be specified.
    The b_key and q_key specifies should be a key in the dictionaries, b_val and
    q_val are the default values, b_mode and q_mode can be "number" - 
    ``float`` will be called to transform the data or "iterable" which will 
    return the length (``len``) of the sequence.
    The B and Q columns can contain only numeric values, thus any data which we 
    wish to store in those columns needs to be converted. The following 
    functions convert data to numeric form. boolean type can also be treated as
    number.
    """

    b_mode = (b_mode or 'number')   # B
    q_mode = (q_mode or 'number')   # Q

    entities = einput(entities, level)

    for entity in entities:
        q_data = eval(q_mode)(entity.xtra.get(q_key), b_val) # occupancy
        b_data = eval(b_mode)(entity.xtra.get(b_key), q_val) # b-factor

        if level != 'A':
            atoms = einput(entity, 'A')
            for atom in atoms:
                if b_key:
                    atom.setBfactor(b_data)
                if q_key:
                    atom.setOccupancy(q_data)
        else:
            if b_key:
                entity.setBfactor(b_data)
            if q_key:
                entity.setOccupancy(q_data)
    PDBWriter(f, entities)  # we try to preserve the headers


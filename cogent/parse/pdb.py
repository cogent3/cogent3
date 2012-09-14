#!/usr/bin/env python
"""PDB parser class and parsing utility functions."""

from re import compile
from numpy import array, linalg

from cogent.data.protein_properties import AA_NAMES
from cogent.core.entity import StructureBuilder, ConstructionWarning, ConstructionError

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

match_coords = compile('^HETATM|ATOM|MODEL')    # start of coordinates
match_trailer = compile('^CONECT')      # end of coordinates

# default PDB format string -> \n neccesery for writelines
# try not to parse a second \n eg. by trying to parse charge
PDB_COORDS_STRING = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n"
PDB_TER_STRING = "%s%5i %-4s%c%3s %c%4i%c\n"

def dict2pdb(d):
    """Transform an atom dictionary into a valid PDB line."""
    (x, y, z) = d['coords']
    args = (d['at_type'], d['ser_num'], d['at_name'], d['alt_loc'],
            d['res_name'][-3:], d['chain_id'], d['res_id'], d['res_ic'],
            x , y , z , d['occupancy'], d['bfactor'], d['seg_id'], d['element'])
    return PDB_COORDS_STRING % args

def dict2ter(d):
    """Transforms an atom dictionary into a valid TER line AFTER it."""
    args = ('TER   ', d['ser_num'] + 1, '   ', ' ',
            d['res_name'][-3:], d['chain_id'], d['res_id'], d['res_ic'])
    return PDB_TER_STRING % args

def pdb2dict(line):
    """Parses a valid PDB line into an atomic dictionary."""
    at_type = line[0:6]
    ser_num = int(line[6:11])   # numbers are ints's?
    at_name = line[12:16]       # " N B "
    at_id = at_name.strip()     # "N B"
    alt_loc = line[16:17]       # always keep \s
    res_name = line[17:20]      # non standard is 4chars long

    chain_id = line[21:22]
    res_id = int(line[22:26].strip())   # pdb requirement int
    res_ic = line[26:27]

    x = line[30:38]             # no float conversion necessery
    y = line[38:46]             # it gets loaded straight into an
    z = line[46:54]             # numpy array.
    occupancy = float(line[54:60])
    bfactor = float(line[60:66])
    seg_id = line[72:76]
    element = line[76:78]

    h_flag = ' '                # this is the default
    if at_type == 'HETATM':     # hetatms get the het flag (it's not used for writing)
        if res_name not in ('MSE', 'SEL'): # SeMet are not ligands
            h_flag = 'H'
            res_name = '%s_%s' % (h_flag, res_name)

    at_long_id = (at_id, alt_loc)
    res_long_id = (res_name, res_id, res_ic)
    coords = array((x, y, z)).astype("double")

    result = {
    'at_type': at_type, 'ser_num': ser_num, 'at_name': at_name,
    'at_id': at_id, 'alt_loc': alt_loc, 'res_name': res_name,
    'chain_id': chain_id, 'res_id': res_id, 'res_ic': res_ic,
    'h_flag': h_flag, 'coords': coords, 'occupancy': occupancy,
    'bfactor': bfactor, 'seg_id': seg_id, 'res_long_id': res_long_id,
    'at_long_id': at_long_id, 'element':element}
    return result

def get_symmetry(header):
    """Extracts symmetry operations from header, either by parsing of
    conversion matrices (SMTRY or BIOMT) or using CCTBX based on the
    space group group data name and a,b,c,alpha,beta,gamma.
    """
    header_result = {}
    for (mode, remark) in (('uc', 'REMARK 290   SMTRY'), ('bio', 'REMARK 350   BIOMT')):
        # parsing the raw-header matrices
        # uc: parsing the symmetry matrices for a given space group
        # bio: parsing the symmetry matrices to construct the biological molecule
        #'REMARK 290   SMTRY3  96 -1.000000  0.000000  0.000000        0.00000'
        # how should we deal with the ORIGx matrix?
        # needed function to check if the SCALEn matrix is correct with
        # respect to the CRYST1 card.
        mx_num = 1
        new_chain = False
        fmx, mxs, mx, cmx = [], [], [], []
        for line in header:
            if line.startswith('SCALE') and mode == 'uc':
                data = map(float, line[6:].split()[:-1])
                fmx.append(data)
            elif line.startswith(remark):
                m_line = map(float, line[20:].split())
                if mx_num != m_line[0] or new_chain:
                    mx.append([ 0., 0., 0., 1.])
                    mxs.append(mx)
                    mx = []
                    if mode == 'bio':
                        if new_chain:
                            new_chain = False
                        else:
                            cmx[-1][1] += 1
                mx_num = m_line[0]
                mx.append(m_line[1:])
            elif line.startswith('REMARK 350 APPLY') and mode == 'bio':
                if cmx: new_chain = True; cmx[-1][1] += 1
                chains = [(c.strip(),) for c in line[42:].split(',')]
                cmx.append([chains, 0])
        mx.append([ 0., 0., 0., 1.]) # finish the last matrix
        mxs.append(mx)
        mxs = array(mxs)        # symmetry matrices
        header_result[mode + "_mxs"] = mxs
        if mode == 'uc':
            fmx = array(fmx)        # fractionalization_matrix
            omx = linalg.inv(fmx)   # orthogonalization_matrix
            header_result["uc_fmx"] = fmx
            header_result["uc_omx"] = omx
        elif mode == 'bio':
            cmx[-1][1] += 1
            header_result[mode + "_cmx"] = cmx
    return header_result

def get_coords_offset(line_list):
    """Determine the line number where coordinates begin."""
    i = 0
    for i, line in enumerate(line_list):
        if match_coords.match(line):
            break
    return i

def get_trailer_offset(line_list):
    """Determine the line number where coordinates end."""
    i = 0
    for i, line in enumerate(line_list):
        if match_trailer.match(line):
            break
    return i

def parse_header(header):
    """Parse parts of the PDB header."""
    id = (compile('HEADER\s{4}\S+.*\S+\s*\S{9}\s+(\S{4})\s*$'), 'id')
    dt = (compile('HEADER\s{4}\S+.*\S+\s+(\S{9})\s+\S{4}\s*$'), 'date')
    nm = (compile('HEADER\s{4}(\S+.*\S+)\s+\S{9}\s+\S{4}\s*$'), 'name')
    mc = (compile('REMARK 280\s+MATTHEWS COEFFICIENT,\s+VM\s+\(ANGSTROMS\*\*3/DA\):\s+(\d+\.\d+)'), 'matthews')
    sc = (compile('REMARK 280\s+SOLVENT CONTENT,\s+VS\s+\(%\):\s+(\d+\.\d+)'), 'solvent_content')
    sg = (compile('REMARK 290\s+SYMMETRY OPERATORS FOR SPACE GROUP:\s+(.*\S)'), 'space_group')
    xt = (compile('REMARK 200\s+EXPERIMENT TYPE\s+:\s+(.*\S)'), 'experiment_type')
    rs = (compile('REMARK   2\s+RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROMS\.'), 'resolution')
    rf = (compile('REMARK   3\s+FREE R VALUE\s+:\s+(\d+\.\d+)'), 'r_free')
    xd = (compile('EXPDTA\s+([\w\-]+).*'), 'expdta')
    c1 = (compile('CRYST1\s+((\d+\.\d+\s+){6})'), 'cryst1')
    ra = (compile('DBREF\s+\S{4}\s+\S\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\S+\s+\d+\s+\d+\s+$'), 'dbref_acc')
    rx = (compile('DBREF\s+\S{4}\s+\S\s+\d+\s+\d+\s+\S+\s+\S+\s+(\S+)\s+\d+\s+\d+\s+$'), 'dbref_acc_full')
    #CRYST1   60.456   60.456   82.526  90.00  90.00  90.00 P 41          4
    #DBREF  1UI9 A    1   122  UNP    Q84FH6   Q84FH6_THETH     1    122             \n'
    tests = [id, dt, nm, mc, sc, sg, xt, rs, rf, xd, c1, ra, rx]
    results = {}
    for line in header:
        for (regexp, name) in tests:
            try:
                results.update({name:regexp.search(line).group(1).strip()})
                continue # taking only the first grep.
            except AttributeError:
                pass
    return results

def parse_coords(builder, coords, forgive=1):
    """Parse coordinate lines."""
    current_model_id = 0
    model_open = False
    current_chain_id = None
    current_seg_id = None
    current_res_long_id = None
    current_res_name = None

    for line in coords:
        record_type = line[0:6]

        if record_type == 'MODEL ':
            builder.initModel(current_model_id)
            current_model_id += 1
            model_open = True
            current_chain_id = None
            current_res_id = None

        elif record_type == 'ENDMDL':
            model_open = False
            current_chain_id = None
            current_res_id = None

        elif record_type == 'ATOM  ' or record_type == 'HETATM':
            if not model_open:
                builder.initModel(current_model_id)
                current_model_id += 1
                model_open = 1
            new_chain = False
            fields = pdb2dict(line)

            if current_seg_id != fields['seg_id']:
                current_seg_id = fields['seg_id']
                builder.initSeg(current_seg_id)

            if current_chain_id != fields['chain_id']:
                current_chain_id = fields['chain_id']
                new_chain = True
                try:
                    builder.initChain(current_chain_id)
                except ConstructionWarning:
                    if not forgive:
                        raise ConstructionError

            if current_res_name != fields['res_name'] or current_res_long_id != fields['res_long_id'] or new_chain:
                current_res_long_id = fields['res_long_id']
                current_res_name = fields['res_name']
                try:
                    builder.initResidue(fields['res_long_id'], fields['h_flag'])
                except ConstructionWarning:
                    if not forgive:
                        raise ConstructionError

                new_chain = False
            try:
                builder.initAtom(fields['at_long_id'], fields['at_name'], fields['ser_num'], \
                                  fields['coords'], fields['occupancy'], fields['bfactor'], \
                                  fields['element'])
            except ConstructionError:
                if not forgive > 1:
                    raise ConstructionError


    return builder.getStructure()

def parse_trailer(trailer):
    return {}

def PDBParser(open_file, structure_id=None, forgive=2):
    """Parse a PDB file and return a Structure object."""
    file_ = open_file.readlines()
    builder = StructureBuilder()

    c_offset = get_coords_offset(file_)
    t_offset = get_trailer_offset(file_)

    raw_coords = file_[c_offset:t_offset]
    raw_trailer = file_[t_offset:]
    raw_header = file_[:c_offset]

    parsed_header = parse_header(raw_header)
    parsed_trailer = parse_trailer(raw_trailer)
    structure_id = (structure_id or parsed_header.get('id'))
    builder.initStructure(structure_id)
    structure = parse_coords(builder, raw_coords, forgive)

    # only X-ray structures will contain crystallographic data
    if parsed_header.get('expdta') == 'X-RAY':
        symetry_info = get_symmetry(raw_header)
        parsed_header.update(symetry_info)

    structure.header = parsed_header
    structure.trailer = parsed_trailer
    structure.raw_header = raw_header
    structure.raw_trailer = raw_trailer

    return structure

#!/usr/bin/env python
"""Parser for the stride commandline tool.
"""


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Production"


def stride_parser(lines):
    """Parses stride output. Returns a ``full_id:{data}`` dictionary.
    
    The parser does not currently support the following (useful) options:
        
        - hydrogen bond assignments (-h)
        - contact order calculation (-k)
        
    The returned data is a dictionary with the following keys:
        'STRIDE_SS': secondary structure,
        'STRIDE_PHI': phi angle,
        'STRIDE_PSI': psi angle,
        'STRIDE_ASA': accessible surface area
    """
    data = {}
    for line in lines:
        if line.startswith('ASG'):
            res_name = line[5:8] # we use 3
            chain_id = line[9]
            if chain_id == '-': # stride ' ' -> '-' rename
                chain_id = ' '
            try:
                res_id = int(float(line[10:15]))
                res_ic = ' '
            except ValueError:
                res_id = int(float(line[10:14]))
                res_ic = line[14]
            ss_code = line[24]
            phi = float(line[43:49])
            psi = float(line[53:59])
            asa = float(line[62:69])
            data[(None, None, chain_id, (res_name, res_id, res_ic), None)] = \
                {
                 'STRIDE_SS': ss_code,
                 'STRIDE_PHI': phi,
                 'STRIDE_PSI': psi,
                 'STRIDE_ASA': asa
                }
    return data


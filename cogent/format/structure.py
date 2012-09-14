#!/usr/bin/env python

import os
from cogent.format.pdb import PDBWriter, PDBXWriter
from cogent.format.xyzrn import XYZRNWriter
from cogent.parse.record import FileFormatError

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

def save_to_filename(entities, filename, format, **kw):
    """Saves a  structure in a specified format into a given file name.
    Arguments:
        - entities: structure or entities to be written
        - filename: name of the structure file
        - format: structure file format
    """
    f = open(filename, 'w')
    try:
        write_to_file(f, entities, format, **kw)
    except Exception:
        try:
            os.unlink(filename)
        except Exception:
            pass
        raise
    f.close()

def write_to_file(f, entities, format, **kw):
    """Saves a structure in a specified format into a given file handle.
    Arguments:
        - entities: structure or entities to be written
        - filename: name of the structure file
        - format: structure file format
    """
    format = format.lower()
    if format not in WRITERS:
        raise FileFormatError("Unsupported file format %s" % format)
    writer = WRITERS[format]
    writer(f, entities, **kw)

# to add a new file format add it's suffix and class name here
WRITERS = {
        'pdb': PDBWriter,
        'pdbx': PDBXWriter,
        'xyzrn': XYZRNWriter
        }


#!/usr/bin/env python
"""Classes for reading macromolecular structure files in different formats."""

import xml.dom.minidom
from cogent.parse import pdb
from cogent.parse.sequence import format_from_filename
from cogent.parse.record import FileFormatError

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


def FromFilenameStructureParser(filename, format=None, **kw):
    """
    Returns a structure parser for a specified format for given filename.
    Arguments:
        - filename: name of the structure file
        - format: the structure file format
    """
    format = format_from_filename(filename, format)
    f = open(filename, 'U')
    return FromFileStructureParser(f, format, **kw)

def FromFileStructureParser(f, format, dialign_recode=False, **kw):
    """
    Returns a structure parser for a specified format for given filename.
    Arguments:
        - filename: name of the structure file
        - format: the structure file format
    """
    if not type(f) is file:
        raise TypeError('%s is not a file' % f)
    format = format.lower()
    if format in XML_PARSERS:
        doctype = format
        format = 'xml'
    else:
        doctype = None
    if format == 'xml':
        source = dom = xml.dom.minidom.parse(f)
        if doctype is None:
            doctype = str(dom.doctype.name).lower()
        if doctype not in XML_PARSERS:
            raise FileFormatError("Unsupported XML doctype %s" % doctype)
        parser = XML_PARSERS[doctype]
    else:
        if format not in PARSERS:
            raise FileFormatError("Unsupported file format %s" % format)
        parser = PARSERS[format]
        source = f
    return parser(source, **kw)

PARSERS = {
        'pdb': pdb.PDBParser,
        }

XML_PARSERS = {
        }


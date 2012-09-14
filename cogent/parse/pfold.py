#!/usr/bin/env python

from string import split,strip
from cogent.struct.rna2d  import Pairs,ViennaStructure
from cogent.parse.column import column_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def pfold_parser(lines):
    """Parser for Pfold output
    """
    tree,lines = tree_struct_sep(lines)
    result = column_parser(lines)

    return result

def tree_struct_sep(lines):
    """Separates tree structure from rest of the data. 

    This is done to get an excepted format for the column_parser
    """
    indx = None
    for line in lines:
        if line.startswith('; ********'):
            indx = lines.index(line)+1
            break
    return lines[:indx],lines[indx:]

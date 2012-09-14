#!/usr/bin/env python

from cogent.parse.bpseq import _parse_residues

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def contrafold_parser(lines):
    """Parser Contarfold output
    
    Returns a list containing sequence and structure(in pair format)
    Ex: [[sequence,[structure]]]

    Tested in tests for bpseq)
    """
    result = []
    seq,struct = _parse_residues(lines,True)
    result.append([seq,struct])
    return result
        

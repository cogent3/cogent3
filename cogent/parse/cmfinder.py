#!/usr/bin/env python

from cogent.util.transform import make_trans
from cogent.struct.rna2d   import wuss_to_vienna, Pairs
from cogent.parse.rfam     import RfamParser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def CMfinderParser(lines):
    """Parser for CMfinder output format

    Parser tested through RfamParser test
    """
    for info, alignment, struct in RfamParser(lines,strict=False):
        struct = wuss_to_vienna(struct)
        pairs = struct.toPairs()
    return [alignment, pairs]
    

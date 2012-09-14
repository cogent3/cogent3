#!/usr/bin/env python

from cogent.parse.ct import ct_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def unafold_parser(lines=None):
    """Parser for unafold output"""
    result = ct_parser(lines)
    return result 

def order_structs(result):
    """Order structures according to energy value

    Order the structures so that the structure with lowest energy is ranked
    first and so on...
    
    Unafold returns results in the same order as the input files
    """
    for i in result:
        i.reverse()
    result.sort()
    #result.reverse()  #to test with the lowest negetiv value as the best struct
    for i in result:
        i.reverse()
    return result

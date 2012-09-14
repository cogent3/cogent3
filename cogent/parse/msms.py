#!/usr/bin/env python
"""Parsers for the MSMS commandline applications.
"""
import numpy as np

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Production"


def parse_VertFile(VertFile):
    """Read the vertex file (with vert extension) into a numpy array.
    
    Arguments:
    
        * VertFile - open vertex file as returned by the ``Msms`` application
          controller.
    
    Returns a numpy array of vertices.
    """
    vertex_list = []
    for line in VertFile.readlines():
        elements = line.split()
        try:
            vertex = map(float, elements[0:3])
        except ValueError:
            continue
        vertex_list.append(vertex)
    return np.array(vertex_list)


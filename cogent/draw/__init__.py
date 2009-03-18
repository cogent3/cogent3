#!/usr/bin/env python

try:
    from linear import colors, \
            TrackDefn, Display, DisplayPolicy, \
            Area, Arrow, BluntArrow, Box, Diamond, Line, Scale

    from legend import Legend
except ImportError: #fails if reportlab not installed
    pass

__all__ = ['dendrogram', 'dotplot', 'legend', 'linear', 'colors', 'TrackDefn',
           'Display', 'DisplayPolicy', 'Area', 'Arrow', 'BluntArrow', 'Box',
           'Diamond', 'Line', 'Scale', 'Legend']

__author__ = ""
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", 
                    "Zongzhi Liu", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

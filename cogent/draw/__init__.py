#!/usr/bin/env python

__all__ = ['dendrogram', 'dotplot', 'legend', 'linear', 'colors', 'TrackDefn',
           'Display', 'DisplayPolicy', 'Area', 'Arrow', 'BluntArrow', 'Box',
           'Diamond', 'Line', 'Scale', 'Legend'] + [
            'arrow_rates', 'dinuc', 'fancy_arrow', 'codon_usage','util']

__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", 
                    "Zongzhi Liu", "Matthew Wakefield", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__status__ = "Production"

try:
    from linear import colors, \
            TrackDefn, Display, DisplayPolicy, \
            Area, Arrow, BluntArrow, Box, Diamond, Line, Scale

    from legend import Legend
except ImportError: #fails if reportlab not installed
    pass

try:
    import cogent.draw.matplotlib
except ImportError:
    pass
else:
    import warnings
    warnings.warn("You still have a cogent/draw/matplotlib subpackage" +
        " present, that will cause problems for matplotlib imports")

#!/usr/bin/env python

sub_modules = ['alltests',
               'test_draw',
               'test_phylo']

for sub_module in sub_modules:
    exec("from %s import %s" % (__name__, sub_module))

__author__ = ""
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
               "Matthew Wakefield", "Andrew Butterfield", "Edward Lang"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

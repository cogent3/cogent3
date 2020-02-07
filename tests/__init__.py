#!/usr/bin/env python

sub_modules = ["test_draw", "test_phylo"]

for sub_module in sub_modules:
    exec("from %s import %s" % (__name__, sub_module))

__author__ = ""
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Rob Knight",
    "Matthew Wakefield",
    "Andrew Butterfield",
    "Edward Lang",
]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

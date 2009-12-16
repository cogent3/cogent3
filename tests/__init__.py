#!/usr/bin/env python

sub_modules = ['alltests',
               'benchmark',
               'benchmark_aligning',
               'test_draw',
               'test_phylo',
               'timetrial']

for sub_module in sub_modules:
    exec ("from %s import %s" % (__name__, sub_module))

__author__ = ""
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Matthew Wakefield", "Andrew Butterfield", "Edward Lang"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

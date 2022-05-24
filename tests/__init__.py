#!/usr/bin/env python
import os
import pathlib


os.chdir(pathlib.Path(__file__).parent)

sub_modules = ["test_draw", "test_phylo"]

for sub_module in sub_modules:
    exec(f"from {__name__} import {sub_module}")

__author__ = ""
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Rob Knight",
    "Matthew Wakefield",
    "Andrew Butterfield",
    "Edward Lang",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

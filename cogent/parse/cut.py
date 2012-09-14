#!/usr/bin/env python
"""Parser for the EMBOSS .cut codon usage table.

Currently just reads the codons and their counts into a dict.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def cut_parser(lines):
    """cut format parser

    Takes lines from a cut file as input.
    returns dict of {codon:count}.
    """
    result = {}
    for line in lines:
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        fields = line.split()
        result[fields[0]] = float(fields[-1])
    return result

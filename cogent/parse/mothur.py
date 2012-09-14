#!/usr/bin/env python
#file cogent.parse.mothur.py
"""Parses Mothur otu list"""

from record_finder import is_empty

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"

    
def parse_otu_list(lines, precision=0.0049):
    """Parser for mothur *.list file

    To ensure all distances are of type float, the parser returns a
    distance of 0.0 for the unique groups.  However, if some sequences
    are very similar, mothur may return a grouping at zero distance.
    What Mothur really means by this, however, is that the clustering
    is at the level of Mothur's precision.  In this case, the parser
    returns the distance explicitly.

    If you are parsing otu's with a non-default precision, you must
    specify the precision here to ensure that the parsed distances are
    in order.

    Returns an iterator over (distance, otu_list)
    """
    for line in lines:
        if is_empty(line):
            continue
        tokens = line.strip().split('\t')

        distance_str = tokens.pop(0)
        if distance_str.lstrip().lower().startswith('u'):
            distance = 0.0
        elif distance_str == '0.0':
            distance = float(precision)
        else:
            distance = float(distance_str)

        num_otus = int(tokens.pop(0))
        otu_list = [t.split(',') for t in tokens]

        yield (distance, otu_list)


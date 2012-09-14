#!/usr/bin/env python
#file cogent.parse.dotur.py
"""Parses various Dotur output formats."""

from record_finder import is_empty

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"


def get_otu_lists(data):
    """Returns list of lists of OTUs given data.
        - data: list of OTUs in following format:
            ['seq_1,seq_2,seq_3','seq_4,seq_5','seq_6','seq_7,seq_8']
    """
    return [i.split(',') for i in data]
    
    
def OtuListParser(lines,ignore=is_empty):
    """Parser for *.list file format dotur result.
        
        - Result will be list of lists with following order:
            [[OTU distance, number of OTUs, [list of OTUs]],
             [OTU distance, number of OTUs, [list of OTUs]],
             [etc...]]
        
    """
    result = []
    if not lines:
        return result
    for line in lines:
        if ignore(line):
            continue
        curr_data = line.strip().split()
        #Get distance.  Replace 'unique' string with 0
        distance = float(curr_data[0].upper().replace('UNIQUE','0'))
        #number of OTUs is second column
        num_otus = int(curr_data[1])
        #remaining columns contain lists of OTUs
        otu_list = get_otu_lists(curr_data[2:])
        result.append([distance,num_otus,otu_list])
    
    return result

    

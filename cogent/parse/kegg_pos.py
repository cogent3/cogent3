#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Release"

"""
Parser for kegg .pos files 

Currently this is quite bare-bones, and primarily useful 
for associating the species name with the results, which is
essential if combining multiple .pos files into a single
database.
"""
# Pos file parsers
def parse_pos_file(fname):
    """Opens fname, extracts pos fields and prepends filename"""
    curr_file = open(fname,"U")
    for line in parse_pos_lines(curr_file,fname):
        yield line

def parse_pos_lines(lines, file_name):
    """Parse lines from a KEGG .pos file, yielding tab-
    delimited strings
    
    file name -- the file name, for deriving the 
    species for the pos file (this is not available within
    the pos file, but important for mapping to other KEGG
    data)
    """
    species_name = file_name.split('/')[-1].rsplit('.',1)[0]
    for line in lines:
        yield species_name + '\t' + line[:-1] + "\n"

if __name__ == '__main__':
    from sys import argv
    filename = argv[1]
    for result_line in parse_pos_file(filename):
        print result_line.strip()

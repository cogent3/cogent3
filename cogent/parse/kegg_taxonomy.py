#!/usr/bin/env python
__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from sys import argv
from string import strip
from os import listdir,path 
from optparse import OptionParser
from datetime import datetime

def parse_kegg_taxonomy(lines):
    """Returns successive taxonomy entries from lines. 

    Format of return value is four levels of taxonomy (sometimes empty),
    unique id, three-letter kegg code, abbreviated name, full name,
    genus, species, and common name if present.
    
    Taxonomic level information is implicit in the number of hashes read
    at the beginning of the last line with hashes. Need to keep track of
    the last level read.

    Each hash line has a number of hashes indicating the level, and a name
    for that taxonomic level. Note that this is not as detailed as the real
    taxonomy in the genome file!

    Maximum taxonomic level as of this writing is 4: exclude any levels more
    detailed than this.

    Each non-taxon line is tab-delimited: has a unique id of some kind, then
    the three-letter KEGG code, then the short name (should be the same as
    the names of the individual species files for genes, etc.), then the
    genus and species names which may have a common name in parentheses
    afterwards.
"""


    max_taxonomy_length = 4
    taxonomy_stack = []
    for line in lines:
        #bail out if it's a blank line
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('#'):    #line defining taxonomic level
            hashes, name = line.split(None, 1)
            name = name.strip()
            level = len(hashes)
            if level == len(taxonomy_stack):    #new entry at same level
                taxonomy_stack[-1] = name
            elif level > len(taxonomy_stack):   #add level: assume sequential
                taxonomy_stack.append(name)
            else:   #level must be less than stack length: truncate
                del taxonomy_stack[level:]
                taxonomy_stack[level-1] = name
        else:   #line defining an individual taxonomy entry
            fields = map(strip, line.split('\t'))
            #add genus, species, and common name as three additional fields
            raw_species_name = fields[-1]
            species_fields = raw_species_name.split()
            if not species_fields:
                print "ERROR"
                print line
            genus_name = species_fields[0]
            if len(species_fields) > 1:
                species_name = species_fields[1]
            else:
                species_name = ''
            #check for common name
            if '(' in raw_species_name:
                prefix, common_name = raw_species_name.split('(', 1)
                common_name, ignored = common_name.split(')', 1)
            else:
                common_name = ''
            output_taxon = taxonomy_stack + \
              ['']*(max_taxonomy_length-len(taxonomy_stack)) \
              + fields + [genus_name, species_name, common_name]
            yield "\t".join(output_taxon) + "\n"

if __name__ == '__main__':
    from sys import argv
    filename = argv[1]
    for result_line in parse_kegg_taxonomy(open(filename,"U")):
        print result_line.strip() 

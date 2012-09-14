#!/usr/bin/env python

from string import strip
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Production"

"""
Parser for KEGG fasta files 

This code is useful for parsing the KEGG .nuc or .pep files
"""

def parse_fasta(lines):
    """lightweight parser for KEGG FASTA format sequences"""
    for label, seq in MinimalFastaParser(lines):
        yield '\t'.join(list(kegg_label_fields(label)) \
          + [seq] + ["\n"])

def kegg_label_fields(line):
    """Splits line into KEGG label fields.

    Format is species:gene_id [optional gene_name]; description.
    """
    fields = map(strip, line.split(None, 1))
    id_ = fields[0] 
    species, gene_id = map(strip, id_.split(':',1))
    #check if we got a description
    gene_name = description = ''
    if len(fields) > 1:
        description = fields[1]
        if ';' in description:
            gene_name, description = map(strip, description.split(';',1))
    return id_, species, gene_id, gene_name, description



if __name__ == '__main__':
    from sys import argv
    filename = argv[1]
    for result_line in parse_fasta(open(filename)):
        print result_line.strip()

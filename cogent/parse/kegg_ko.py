#!/usr/bin/env python

from sys import argv
from string import strip
from os import listdir,path 
from optparse import OptionParser
from datetime import datetime
import tarfile

_author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


"""A parser for the KEGG 'ko' file 
containing information on KEGG orthology groups and their associated pathways.
"""

def parse_ko_file(filepath,dir_prefix=None,debug = True): 
    """Parse the NCBI KO file lines, and output several tab-delimited files
    
    filepath - the full filepath to the input KO file from KEGG
    dir_prefix - the directory to which tab-delimited output files will be saved.
    debug - if set to True, pring debugging output to the screen
    """
    
    lines = open(filepath,"U")
    ko_gene_fname = 'ko_to_gene.tab'
    ko_fname = 'ko.tab'
    ko_pathway_fname = 'ko_to_pathway.tab'
    pathway_fname = 'pathway.tab'
    ko_cog_fname = 'ko_to_cog.tab'
    ko_cazy_fname = 'ko_to_cazy.tab'
    ko_go_fname = 'ko_to_go.tab'
    fnames = [ko_gene_fname, ko_fname, ko_pathway_fname,\
              pathway_fname, ko_cog_fname, ko_cazy_fname,\
              ko_go_fname]
    
    if dir_prefix:
        fnames = [dir_prefix + '/' + f for f in fnames]
    
    if debug:
        for res_fp in fnames:
            print "Outputting parsed info to: %s" %(res_fp) 
    
    ko_gene, ko, ko_pathway, pathway, ko_cog, ko_cazy, ko_go = \
        [open(i, 'w') for i in fnames]

    #figure out what fields we want (and get them), and get pathway data
    fields = ['ENTRY', 'NAME', 'DEFINITION']
    ko_to_pathway = {}
    for rec in parse_ko(lines):
        ko.write('\t'.join([rec.get(f,'') for f in fields]))
        ko.write('\n')
        entry = rec['ENTRY']
        if 'GENES' not in rec:
            continue    #apparently, some records don't have genes...
        genes = rec['GENES']

        for species, gene_list in genes.items():
            for g in gene_list:
                ko_gene.write('%s\t%s:%s\n' % (entry, species.lower(), g))
        
        if 'CLASS' not in rec:
            continue    #apparently they also lack classes...
        
        ko_to_pathway[entry] = rec['CLASS']

        dblinks = rec.get('DBLINKS', None)
        
        if dblinks:
            cogs = dblinks.get('COG', None)
            cazy = dblinks.get('CAZy', None)
            go   = dblinks.get('GO', None)
	    
        if cogs:
            for c in cogs:
                ko_cog.write("%s\t%s\n" % (entry, c))
        if go:
            for g in go:
                ko_go.write("%s\t%s\n" % (entry, g))
        if cazy:
            for c in cazy:
                ko_cazy.write("%s\t%s\n" % (entry,c))
    #postprocess the ko_to_pathway data to find out what the pathway terms
    #are and to write them out into a join file
    max_terms = 10
    unique_recs = {}    #will hold tuple(fields) -> unique_id
    curr_uid = 0
    for ko, classes in ko_to_pathway.items():
        for (id_, fields) in classes:
            if fields not in unique_recs:
                unique_recs[fields] = curr_uid
                fields_for_output = fields[:]
                if len(fields_for_output) > max_terms:
                    fields_for_output = fields_for_output[:max_terms]
                elif len(fields_for_output) < max_terms:
                    fields_for_output += \
                        ('',)*(max_terms - len(fields_for_output))
                pathway.write('\t'.join((str(curr_uid),str(id_)) +\
                    fields_for_output)+'\n')
                curr_uid += 1
            uid = unique_recs[fields]
            ko_pathway.write(str(ko)+ '\t'+ str(uid) + '\n')


def make_tab_delimited_line_parser(columns_to_convert):
    """Generates a function that parses a tab-delimited line
    
    columns_to_convert:  a list of column indexes to convert into integers
    by splitting on ':' and taking the second entry (e.g. to convert listings
    like GO:0008150 to 0008150 or ncbi-gi:14589889 to 14589889)"""

    def parse_tab_delimited_line(line):
        """Parse a tab-delimited line taking only the second item of cols %s""" %\
          str(columns_to_convert)
        fields = line.split("\t")
        for index in columns_to_convert:
            fields[index] = fields[index].split(":")[1]
        return "\t".join(fields)

    return parse_tab_delimited_line

def ko_default_parser(lines):
    """Handle default KEGG KO entry lines
     lines -- default format of space separated lines.
             Examples include the NAME and DEFINITION 
             entries
     
     Strips out newlines and joins lines together."""

    return ' '.join(map(strip, lines)).split(None, 1)[1]

def ko_first_field_parser(lines):
    """Handles KEGG KO entries where only the first field is of interest

    For example, ENTRY fields like:
    'ENTRY       K01559             KO\n'
    
    Strips out newlines and joins lines together for the first field only."""
    return ' '.join(map(strip, lines)).split()[1]

def delete_comments(line):
    """Deletes comments in parentheses from a line."""
    fields = line.split(')')
    result = []
    for f in fields:
        if '(' in f:
            result.append(f.split('(',1)[0])
        else:
            result.append(f)
    return ''.join(result)

def ko_colon_fields(lines, without_comments=True):
    """Converts line to (key, [list of values])
    
    lines -- colon fields such as DBLINKS or GENES 
    in the KEGG KO file.

    Example:
    '            BXE: Bxe_B0037 Bxe_C0683 Bxe_C1002 Bxe_C1023\n'
    """
    merged = ' '.join(map(strip, lines))
    if without_comments:
        merged = delete_comments(merged)
    key, remainder = merged.split(':',1)
    vals = remainder.split()
    return key, vals

def ko_colon_delimited_parser(lines, without_comments=True):
    """For lines of the form LABEL: id: values.

    Returns dict of id:values.
    """
    first_line = lines[0]
    without_first_field = first_line.split(None, 1)[1]
    data_start = len(first_line) - len(without_first_field)
    result = {}
    curr = []
    for line in lines:
        line = line[data_start:]
        if line[0] != ' ':  #start of new block
            if curr:
                key, vals = ko_colon_fields(curr, without_comments)
                result[key] = vals
                curr = []
        curr.append(line)
    if curr:
        key, vals = ko_colon_fields(curr, without_comments)
        result[key] = vals
    return result

def _is_new_kegg_rec_group(prev, curr):
    """Check for irregular record group terminators"""
    
    return curr[0].isupper() and not prev.endswith(';') and \
    not curr.startswith('CoA biosynthesis') and not prev.endswith(' and') and \
        not prev.endswith('-') and not prev.endswith(' in') and not \
        prev.endswith(' type') and not prev.endswith('Bindng') and not \
        prev.endswith('Binding')

def group_by_end_char(lines, end_char = ']', \
    is_new_rec=_is_new_kegg_rec_group):
    """Yields successive groups of lines that end with the specified char.
    
    Note: also returns the last group of lines whether or not the end char
    is present.
    """
    curr_lines = []
    prev_line = ''
    for line in lines:
        stripped = line.strip()
        #unfortunately, not all records in kegg actually end with the
        #terminator, so need to check for termination condition
        if is_new_rec(prev_line, stripped):
            if curr_lines:
                yield curr_lines
            curr_lines = []
        #if the line ends with the character we're looking for, assume we've
        #found a new record
        if stripped.endswith(end_char):
            yield curr_lines + [line]
            curr_lines = []
        else:
            curr_lines.append(line)
        prev_line = stripped
    if curr_lines:
        yield curr_lines

def class_lines_to_fields(lines):
    """Converts a list of lines in a single pathway within one KO class definition.
    """
    rec = ' '.join(map(strip, lines))
    #need to split off the class declaration if it is present
    if rec.startswith('CLASS'):
        rec = rec.split(None,1)[1]
    #figure out if it has an id and process accordingly
    if rec.endswith(']'):
        rec, class_id = rec.rsplit('[', 1)
        class_id = class_id[:-1]
    else:
        class_id = None
    rec_fields = map(strip, rec.split(';'))
    return class_id, tuple(rec_fields)

def ko_class_parser(lines, without_comments='ignored'):
    """For the CLASS declaration lines.

    These take the form of multi-line semicolon-delimited fields (where
    each field is a successive entry in the KEGG pathway hierarchy), ending
    in a field of the form [PATH:ko00071].

    Strategy:
    - iterate over groups of lines that end in ] (each represents one pathway)
    - for each line:
        - split off and extract the pathway id
        - split the rest of the terms on semicolon
        - return a tuple of (pathway_id, [terms_in_order])
    
    Don't consolidate the terms in this parser because each KO group has
    its own class declaration so we would have to merge them for each class:
    instead, merge at higher level.
    """
    for group in group_by_end_char(lines):
        yield class_lines_to_fields(group)

def parse_ko(lines):
    """Parses a KO record into fields."""
    
    # Here we define records by their category 
    # to allow parsers to be reused on 
    # similar entries.
    
    default_fields = ['NAME', 'DEFINITION']
    colon_fields = ['DBLINKS', 'GENES']
    first_field_only = ['ENTRY']
    class_fields = ['CLASS']
    
    for rec in ko_record_iterator(lines):
        split_fields = ko_record_splitter(rec)
        result = {}
        for k, v in split_fields.items():
            if k in default_fields:
                result[k] = ko_default_parser(v)
            elif k in colon_fields:
                result[k] = ko_colon_delimited_parser(v)
            elif k in first_field_only:
                result[k] = ko_first_field_parser(v)
            elif k in class_fields:
                result[k] = list(ko_class_parser(v))
        yield result

#parse_ko: lightweight standalone ko parser
def ko_record_iterator(lines):
    """Iterates over KO records, delimited by '///'"""
    curr = []
    for line in lines:
        if line.startswith('///') and curr:
            yield curr
            curr = []
        else:
            curr.append(line)
    if curr:
        yield curr

def ko_record_splitter(lines):
    """Splits KO lines into dict of groups keyed by type."""
    result = {}
    curr_label = None
    curr = []
    i = 0
    for line in lines:
        i+= 1
        if line[0] != ' ':
            if curr_label is not None:
                result[curr_label] = curr
            fields = line.split(None, 1)
            
            # Annoyingly, can have blank REFERENCE lines
            # Lacking PMID, these still have auth/title info, however...
            if len(fields) == 1:
                curr_label = fields[0]
                curr_line = ''
            else:    
                curr_label, curr_line = fields
            curr = [line]
        else:
            curr.append(line)
    if curr:
        result[curr_label] = curr
    return result

if __name__ == '__main__':
    from sys import argv
    filename = argv[1]
    out_dir = argv[2]
    parse_ko_file(filename, \
      dir_prefix = out_dir, \
      debug = True)

#!/usr/bin/env python
#bpseq.py
"""Provides parser for bpseq files downloaded from Gutell's comparative 
RNA website (CRW): http://www.rna.icmb.utexas.edu/

The file format is:

Filename: d.16.b.E.coli.bpseq
Organism: Escherichia coli
Accession Number: J01695
Citation and related information available at http://www.rna.icmb.utexas.edu
1 A 0
2 A 0
3 A 0
4 U 0
5 U 0
6 G 0
7 A 0
8 A 0
9 G 25
10 A 24
11 G 23
12 U 22
13 U 21
14 U 0

So, header of four lines (Filename, Organism, Accession Number, Citation)
Sequence, structure information in tuples of residue position, residue name,
residue partner. The residue partner is 0 if the base is unpaired.

Numbering is 1-based!
"""
from __future__ import division
from string import strip
from cogent.struct.rna2d import Vienna, Pairs
from cogent.struct.knots import opt_single_random
from cogent.core.info import Info
from cogent.core.sequence import RnaSequence

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class BpseqParseError(Exception):
    """Exception raised when an error occurs during parsing a bpseq file"""
    pass

def parse_header(header_lines):
    """Return Info object from header information.
   
    header_lines -- list of lines or anything that behaves like it.

    Parses only the first three header lines with Filename, Organism, and
    Accession number. In general lines that contain a colon will be parsed.
    There's no error checking in here. If it fails to split on ':', the 
    information is simply not added to the dictionary. The expected format
    for header lines is "key: value". The citation lane is parsed differently. 
    """
    info = {}
    for line in header_lines:
        if line.startswith('Citation'):
            info['Citation'] = line.split()[-1].strip()
        elif ':' in line:
            try:
                field, value = map(strip,line.split(':',1))
                info[field] = value
            except ValueError:
                #no interesting header line
                continue
        else:
            continue
    return Info(info)

def construct_sequence(seq_dict):
    """Construct RnaSequence from dict of {pos:residue}.

    seq_dict -- dictionary of {position: residue}
    
    Checks whether the first residue is 0. Checks whether all residues
    between min and max index are present. No checking on validity of
    residue symbols. 
    """
    all_pos = seq_dict.keys()
    min_pos, max_pos = min(all_pos), max(all_pos)
    if min_pos != 0:
        raise BpseqParseError(\
            "Something went wrong with adjusting the numbering")
    # make sure all positions are in the dictionary
    for idx in range(min_pos, max_pos+1):
        if idx not in seq_dict:
            raise BpseqParseError(\
                "Description of residue with index %s is missing"%(idx))
    seq = [] 
    for idx in range(min_pos, max_pos+1):
        seq.append(seq_dict[idx])
    # Return as a simple string
    return ''.join(seq)

def parse_residues(residue_lines, num_base, unpaired_symbol):
    """Return RnaSequence and Pairs object from residue lines.

    residue_lines -- list of lines or anything that behaves like it. 
        Lines should contain:
        residue_position, residue_identiy, residue_partner.
    num_base -- int, basis of the residue numbering. In bpseq files from
        the CRW website, the numbering starts at 1.
    unpaired_symbol -- string, symbol in the 'partner' column that indicates
        that a base is unpaired. In bpseq files from the CRW website, the
        unpaired_symbol is '0'. This parameter should be a string to allow
        other symbols that can't be casted to an integer to indicate
        unpaired bases.
    
    Checks for double entries both in the sequence and the structure, and
    checks that the structre is valid in the sense that if (up,down) in there,
    that (down,up) is the same.
    """
    #create dictionary/list for sequence and structure
    seq_dict = {}
    pairs = Pairs()
    
    for line in residue_lines:
        try:
            pos, res, partner = line.strip().split()
            if partner == unpaired_symbol:
                # adjust pos, not partner
                pos = int(pos) - num_base
                partner = None
            else:
                # adjust pos and partner
                pos = int(pos) - num_base
                partner = int(partner) - num_base
            pairs.append((pos,partner))
            
            #fill seq_dict
            if pos in seq_dict:
                raise BpseqParseError(\
                    "Double entry for residue %s (%s in bpseq file)"\
                    %(str(pos), str(pos+1)))
            else:
                seq_dict[pos] = res
        
        except ValueError:
            raise BpseqParseError("Failed to parse line: %s"%(line))
    
    #check for conflicts, remove unpaired bases 
    if pairs.hasConflicts():
        raise BpseqParseError("Conflicts in the list of basepairs")
    pairs = pairs.directed()
    pairs.sort()
    
    # construct sequence from seq_dict
    seq = RnaSequence(construct_sequence(seq_dict))
    
    return seq, pairs

def MinimalBpseqParser(lines):
    """Separate header and content (residue lines).
    
    lines -- a list of lines or anything that behaves like that.

    The standard bpseq header (from the CRW website) is recognized. Also,
    lines that contain a colon are accepted as header lines.
    Header lines that aren't accepted as header, but that can be split into 
    three parts are residue lines (sequence and structure description). 
    Lines that don't fall into any of these categories are ignored.
    """
    result = {'HEADER':[], 'SEQ_STRUCT':[]}
    
    for line in lines:
        if line.startswith('Filename') or line.startswith('Organism') or\
            line.startswith('Accession') or line.startswith('Citation') or\
            ":" in line:
            result['HEADER'].append(line.strip())
        elif len(line.split()) == 3:
            result['SEQ_STRUCT'].append(line.strip())
        else:
            continue #unknown
    return result

def BpseqParser(lines, num_base=1, unpaired_symbol='0'):
    """Return RnaSequence and structure (Pairs object) specified in file.

    lines -- filestream of bpseq file. File should contain a single record.
    num_base -- int, basis of the residue numbering. In bpseq files from
        the CRW website, the numbering starts at 1.
    unpaired_symbol -- string, symbol in the 'partner' column that indicates
        that a base is unpaired. In bpseq files from the CRW website, the
        unpaired_symbol is '0'. This parameter should be a string to allow
        other symbols that can't be casted to an integer to indicate
        unpaired bases.

    Bpseq file looks like this:
    
    Filename: d.16.b.E.coli.bpseq
    Organism: Escherichia coli
    Accession Number: J01695
    Citation and related information available at http://www....
    1 A 0
    2 A 0
    3 A 0
    4 U 0
    5 U 0
    6 G 0
    7 A 0
    8 A 0
    9 G 25
    10 A 24
    11 G 23
    12 U 22
    13 U 21
        
    So, 4 header lines, followed by a list of residues.
    Position (indexed to 1), residue, partner position
    """
    # separate header and residue lines
    grouped_lines = MinimalBpseqParser(lines)
    # parse header and seq/struct separately
    header_info = parse_header(grouped_lines['HEADER'])
    seq, struct = parse_residues(grouped_lines['SEQ_STRUCT'],\
        num_base, unpaired_symbol)
    #add header info to the sequence as Info object
    seq.Info = header_info
    
    return seq, struct

# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def bpseq_specify_output(lines, num_base=1, unpaired_symbol='0',
    return_vienna=False, remove_pseudo=False,\
    pseudoknot_function=opt_single_random):
    """Return Vienna structure of Pairs object with or without pseudoknots
    
    lines -- filestream of bpseq file. File should contain a single record.
    num_base -- int, basis of the residue numbering. In bpseq files from
        the CRW website, the numbering starts at 1.
    unpaired_symbol -- string, symbol in the 'partner' column that indicates
        that a base is unpaired. In bpseq files from the CRW website, the
        unpaired_symbol is '0'. This parameter should be a string to allow
        other symbols that can't be casted to an integer to indicate
        unpaired bases.
    return_vienna -- boolean, if True, a ViennaStructure object is returned,
        if False, a Pairs object is returned. If return_vienna is True,
        pseudoknots need to be removed from the structure.
    remove_pseudo -- boolean, if True, pseudoknots will be removed from
        the structure.
    pseudoknot_function -- function that takes a Pairs object as input and
        returns a nested version of the structure. The function should return
        a single nested structure, not a list of structures. Default is 
        opt_single_random, which retuns a single nested structure (it picks
        one at random in case of multiple structures with the maximum number of 
        base pairs). This default is chosen to assure the code always returns
        something. In case the experiment needs to be reproducible,
        a random choice isn't the best one to make, and one should
        use a different pseudoknot removal function. See struct/knots.py
        for documentation.
    """
    seq, pairs = BpseqParser(lines, num_base, unpaired_symbol)

    if pairs.hasPseudoknots() and (remove_pseudo or return_vienna):
        pairs = pseudoknot_function(pairs)
    if return_vienna:
        v = pairs.toVienna(len(seq))
        return seq, v
    else:
        return seq, pairs


if __name__ == "__main__":
    from sys import argv
    seq, struct = BpseqParser(open(argv[1]))
    print seq
    print seq.Info
    print struct

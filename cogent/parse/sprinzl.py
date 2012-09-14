#/usr/bin/env python
"""Parsers for the Sprinzl tRNA databases.
"""
from cogent.util.misc import InverseDict
from string import strip, maketrans
from cogent.core.sequence import RnaSequence
from cogent.core.info import Info as InfoClass

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Jeremy Widmann", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def Rna(x, Info=None):
    if isinstance(x, list):
        x = ''.join(x)
    if Info is None:
        Info = {}
    return RnaSequence(x.upper().replace('T','U'), Info=InfoClass(Info))

SprinzlFields =['Accession', 'AA', 'Anticodon', 'Species', 'Strain']
            
def OneLineSprinzlParser(infile):
    """Returns successive records from the tRNA database. First line labels.
   
    This was the first attempt at the parser, and requires quite a lot of
    preprocessing. Use SprinzlParser for something more general.
   
    Works on a file obtained by the following method:
    1. Do the default search.
    2. Show all the columns and autofit them.
    3. Delete the first column of numbers and all blank columns.
    4. Name the first 5 columns "Accession, AA, Anticodon, Species, Strain".
    5. Save the worksheet as plain text.
    """
    first = True
    for l in infile:
        line = l.strip()
        if not line:
            continue
        fields = line.split('\t')
        if first:   #label line
            label_fields = fields[5:]
            labels = InverseDict(enumerate(label_fields))
            first = False
        else:
            info = dict(zip(SprinzlFields, map(strip, fields[0:5])))
            info['Labels'] = labels
            yield Rna(map(strip, fields[5:]), Info=info)

GenomicFields = ['', 'Accession', 'AA', '', 'Anticodon', '', 'Species', \
    '', '', '', '', '', '', '', '', '', 'Strain', '', '', '', 'Taxonomy']

def _fix_structure(fields, seq):
    """Returns a string with correct # chars from db struct line.

    fields should be the result of line.split('\t')
    
    Implementation notes:
        Pairing line uses strange format: = is pair, * is GU pair, and
        nothing is unpaired. Cells are not padded out to the start or end of 
        the sequence length, presumably to infuriate the unwary.

        I don't _think_ it's possible to convert these into ViennaStructures
        since we don't know where each helix starts and ends, and the lengths
        of each piece can vary. I'd be happy to be proven wrong on this...

        For some reason, _sometimes_ spaces are inserted, and _sometimes_
        the cells are left entirely blank. Also, when there's a noncanonical
        pair in the helix, the helix is broken into two pieces, so counting
        pieces isn't going to work for figuring out the ViennaStructure.

    Expects as input the sequence and the raw structure line.
    """
    num_blanks = 4
    pieces = fields[num_blanks:]
    result = ['.'] * len(seq)
    for i, p in enumerate(pieces):
        if p and (p != ' '):
            result[i] = p
    return ''.join(result)

def _fix_sequence(seq):
    """Returns string where terminal gaps are replaced with terminal CCA.
    
        Some of the sequence in the Genomic tRNA Database have gaps where the
        acceptor stem (terminal CCA) should be.  This function checks the
        number of terminal gaps and replaces with appropriate part of terminal
        CCA.
    """
    if seq.endswith('---'):
        seq = seq[:-3]+'CCA'
    elif seq.endswith('--'):
        seq = seq[:-2]+'CA'
    elif seq.endswith('-'):
        seq = seq[:-1]+'A'
    return seq

def GenomicSprinzlParser(infile,fix_sequence=False):
    """Parser for the Genomic tRNA Database.

    Assumes the file has been prepared by the following method:
    1. Set all search fields to empty.
    2. Check all the results fields.
    3. Perform the search (this takes a while).
    4. Save the results worksheet as tab-delimited text.

    Note that the alignment length is supposed to be 99 bases, but not all the
    sequences have been padded out with the correct number of hyphens.
    """
    num_blanks = 4
    first = True
    for l in infile:
        #skip blank lines
        line = l.rstrip()
        if not line:
            continue
        fields = line.split('\t')
        if first:   #label line
            #for unknown reasons, some of the field headers have '.' instead
            #of '0', e.g. '7.' instead of '70'.
            line = line.replace('.', '0')
            fields = line.split('\t')
            labels = InverseDict(enumerate(fields[num_blanks:]))
            first = False
            offset = 0
        else:       #expect 3 record lines at a time
            if offset == 0:     #label line
                info = dict(zip(GenomicFields, map(strip, fields)))
                #add in the labels
                info['Labels'] = labels
                #convert the taxonomy from a string to a list
                info['Taxonomy'] = map(strip, info['Taxonomy'].split(';'))
                #convert the anticodon into RNA
                info['Anticodon'] = Rna(info['Anticodon'])
                #get rid of the empty fields
                del info['']
            elif offset == 1:   #sequence line
                raw_seq = ''.join(map(strip, fields))
                #for some reason, there are underscores in some sequences
                raw_seq = raw_seq.replace('_', '-')
                if fix_sequence:
                    raw_seq = _fix_sequence(raw_seq)
                seq = Rna(raw_seq, Info=info)
            elif offset == 2:   #structure line
                seq.Pairing = _fix_structure(fields, seq)
                yield seq
            #figure out which type of line we're expecting next
            offset += 1
            if offset > 2:
                offset = 0

def get_pieces(struct, splits):
    """Breaks up the structure at fixed positions, returns the pieces.
    
    struct: structure string in sprinzl format
    splits: list or tuple of positions to split on

    This is a helper function for the sprinzl_to_vienna function.
    
    struct = '...===...===.'
    splits = [0,3,7,-1,13]
    pieces -> ['...','===.','..===','.']
    """
    pieces = []
    for x in range(len(splits)-1):
        pieces.append(struct[splits[x]:splits[x+1]])
    return pieces

def get_counts(struct_piece):
    """Returns a list of the lengths or the paired regions in the structure.
    
    struct_pieces: string, piece of structure in sprinzl format

    This is a helper function for the sprinzl_to_vienna function

    struct_piece = '.===.=..'
    returns [3,1]
    """
    return map(len, filter(None, [i.strip('.') for i in \
        struct_piece.split('.')]))
        
def sprinzl_to_vienna(sprinzl_struct):
    """Constructs vienna structure from sprinzl sec. structure format
    
    sprinzl_struct: structure string in sprinzl format
    
    Many things are hardcoded in here, so if the format or the alignment
    changes, these values have to be adjusted!!!
    The correctness of the splits has been tested on the GenomicDB 
    database from Jan 2006, containing 8163 sequences.
    """
    assert len(sprinzl_struct) == 99
    gu='*'
    wc='='
    splits = [0,8,19,29,38,55,79,-11,len(sprinzl_struct)]
    direction = ['(','(',')','(',')','(',')',')']
    
    #get structural pieces
    s = sprinzl_struct.replace(gu,wc)
    pieces = get_pieces(s, splits)
    assert len(pieces) == len(splits)-1
    
    #get counts of structured regions in each piece, check validity
    counts = map(get_counts,pieces)
    pairs = [(0,-1),(1,2),(3,4),(5,6)]
    for i,j in pairs:
        assert sum(counts[i]) == sum(counts[j])
    #check counts matches directions 
    assert len(counts) == len(direction)
    
    #construct string string of brackets
    brackets = []
    for lengths, br in zip(counts,direction):
        for l in lengths:
            brackets.append(l*br)
    brackets = ''.join(brackets)
   
    #build vienna structure
    vienna = []
    x=0
    for sym in s:
        if sym == '.':
            vienna.append(sym)
        else:
            vienna.append(brackets[x])
            x += 1
    return ''.join(vienna)

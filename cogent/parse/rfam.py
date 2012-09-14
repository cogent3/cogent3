#!/usr/bin/env python
"""Provides a parser for Rfam format files.
"""
from string import strip
from cogent.parse.record import RecordError
from cogent.parse.record_finder import DelimitedRecordFinder
from cogent.parse.clustal import ClustalParser
from cogent.core.sequence import RnaSequence as Rna
from cogent.core.sequence import DnaSequence as Dna
from cogent.core.sequence import Sequence
from cogent.core.moltype import BYTES
from cogent.core.info import Info
from cogent.struct.rna2d import WussStructure
from cogent.util.transform import trans_all,keep_chars
from cogent.core.alignment import Alignment, DataError, SequenceCollection

__author__ = "Sandra Smit and Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Gavin Huttley",
                    "Rob Knight", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

def is_empty_or_html(line):
    """Return True for HTML line and empty (or whitespace only) line.

    line -- string

    The Rfam adaptor that retrieves records inlcudes two HTML tags in
    the record. These lines need to be ignored in addition to empty lines. 
    """
    if line.startswith('<pre') or line.startswith('</pre'):
        return True
    return (not line) or line.isspace()

Sequence = BYTES.Sequence
RfamFinder = DelimitedRecordFinder('//', ignore=is_empty_or_html)

def load_from_clustal(data, seq_constructor=Sequence, strict=True):
    recs = [(name, seq_constructor(seq, )) for name, seq in\
        ClustalParser(data, strict)]
    lengths = [len(i[1]) for i in recs]
    if lengths and max(lengths) == min(lengths):
        return Alignment(recs, MolType=BYTES)
    else:
        return SequenceCollection(recs, MolType=BYTES)

#all fields concerning the references are translated to None, except for 
# the MedLine ID, so that we can lookup the information if needed.
#RC = Reference comment
#RN = Reference Number
#RT = Reference Title
#RA = Reference Author
#RL = Reference Location
# The None fields are filtered out later
_field_names = {'AC':'Rfam',\
                'ID':'Identification',\
                'DE':'Description',\
                'AU':'Author',\
                'SE':'AlignmentSource',\
                'SS':'StructureSource',\
                'BM':'BuildCommands',\
                'GA':'GatheringThreshold',\
                'TC':'TrustedCutoff',\
                'NC':'NoiseCutoff',\
                'TP':'FamilyType',\
                'SQ':'Sequences',\
                'PI': 'PreviousIdentifications',\
                'DC': 'DatabaseComment',\
                'DR': 'DatabaseReference',\
                'RC': None,\
                'RN': None,\
                'RM': 'MedlineRef',\
                'RT': None,\
                'RA': None,\
                'RL': None,\
                'CC': 'Comment'}
                

def HeaderToInfo(header,strict=True):
    """Returns an Info object constructed from the header lines.

    Header is a list of lines that contain header information.
    Fields that can occur multiple times in a header are stored in a list.
    Fields that (should) occur only once are stored as a single value
    Comments are joined by ' ' to one field.
    Fields concerning the references are ignored, except for MedLine ID.
    """
    # construct temporary dictionary containing all original information
    initial_info = {}
    for line in header:
        line = line.strip()
        if not line:
            continue
        try:
            init,label,content = line.split(' ',2)
            if not init == '#=GF' or len(label) != 2:
                raise RecordError
        except:
            if strict:
                raise RecordError, "Failed to extract label and content " +\
                    "information from line %s"%(line)
            else:
                continue
        if label in ['BM','DR','RM','CC']:
            if label in initial_info:
                initial_info[label].append(content.strip())
            else:
                initial_info[label] = [content.strip()]
        else:
            initial_info[label] = content.strip()
            
    # transform initial dict into final one
    # throw away useless information; group information
    final_info={}
    for key in initial_info.keys():
        name = _field_names.get(key,key)
        if name == 'Comment':
            value = ' '.join(initial_info[key])
        else:
            value = initial_info[key]
        final_info[name] = value
    
    return Info(final_info)

def NameToInfo(sequence, strict=True):
    """Returns an Info object constructed from the sequence Name

    sequence: Sequence object with a Name attribute

    The label will be split on Genbank acc. no. and sequence coordinates.
    The coordinates will be shifted one position, since in Python the first
        position is 0.
    """
    #adjust label
    label = sequence.Name
    try:
        gb, pos = label.split('/',1) #split genbank label and pos
        if not gb:
            gb = None
        if not pos:
            pos = None
    except: #unable to split, so string doesn't contain '/'
        if strict:
            raise RecordError, "Failed to extract genbank id and positions" +\
            " from label %s"%label
        else:
            gb = None
            pos =None
    if pos:
        try:
            start, end = pos.split('-',1) #split start and end pos
        except:
            if strict:
                raise RecordError,\
                    "Failed to extract genbank id and positions from label %s"\
                    %label
            else:
                start = None
                end = None
    else:
        start = None
        end = None
    if start:
        # adjust start position to do the correct thing in python
        # see comment in docstring
        start = int(start)-1
    if end:
        end = int(end)
    info = Info({'GenBank':gb,'Start':start,'End':end})
    return info

def is_header_line(line):
    """Returns True if line is a header line"""
    return line.startswith('#=GF')

def is_seq_line(line):
    """Returns True if line is a sequence line"""
    return bool(line) and (not line[0].isspace()) and \
    (not line.startswith('#')) and (not line.startswith('//'))

def is_structure_line(line):
    """Returns True if line is a structure line"""
    return line.startswith('#=GC SS_cons')

def ChangedRnaSequence(data, seq_constructor=Rna):
    """Returns new RNA Sequence object, replaces dots with dashes in sequence.
    """
    return seq_constructor(str(data).replace('.','-'))

def ChangedDnaSequence(data, seq_constructor=Dna):
    """Returns new RNA Sequence object, replaces dots with dashes in sequence.
    """
    return seq_constructor(str(data).replace('.','-'))

def ChangedSequence(data, seq_constructor=Sequence):
    """Returns new Sequence object, replaces dots with dashes in sequence.
    """
    return seq_constructor(str(data).replace('.','-'))


def MinimalRfamParser(infile,strict=True,seq_constructor=ChangedRnaSequence):
    """Yield successive sequences as (header, sequences, structure) tuples.
    
    header is a list of header lines
    sequences is an Alignment object. Sequences are objects keyed by the
        original labels in the database.
    structure is a WussStructure
    """
    for record in RfamFinder(infile):
        header = []
        sequences = []
        structure = []
        for line in record:
            if is_header_line(line):
                header.append(line.strip())
            elif is_seq_line(line):
                sequences.append(line)
            elif is_structure_line(line):
                structure.append(line)
            else:
                continue
        #sequence and structure are required. 
        #for example when looking at the stockholm format of just one family
        if not sequences or not structure:
            if strict:
                error = 'Found record with missing element(s): '
                if not sequences:
                    error += 'sequences '
                if not structure:
                    error += 'structure '
                raise RecordError, error
            else:
                continue
        #join all sequence parts together, construct label
        try:
            new_seqs = load_from_clustal(sequences,strict=strict,
                seq_constructor=seq_constructor)
            sequences = new_seqs
        except (DataError, RecordError), e:
            if strict:
                raise RecordError, str(e)
            else:
                continue

        #construct the structure
        try:
            res = load_from_clustal(structure, strict=strict)
            assert len(res.NamedSeqs) == 1 #otherwise multiple keys
            structure = res.NamedSeqs['#=GC SS_cons']
        except (RecordError, KeyError, AssertionError), e:
            if strict:
                raise RecordError,\
                    "Can't parse structure of family: %s"%(str(header))
            else:
                structure = None
        yield header, sequences, structure
                
def RfamParser(lines, seq_constructor=ChangedRnaSequence, label_constructor=\
    HeaderToInfo,struct_constructor=WussStructure,strict=True,verbose=False):
    """Yields (family_info, sequences, structure).

    Treats lines as a stream of Rfam records.
    Family_info is the general information about the alignment.
    Sequences is an Alignment object. Each sequence has its own Info
        object with Genbank ID etc. Sequences are keyed by the original 
        label in the database.
    Structure is the consensus structure of the alignment, in Wuss format
    """
    for header, alignment, structure in MinimalRfamParser\
        (lines,strict=strict,seq_constructor=seq_constructor):
        if strict:
            try:
                family_info = label_constructor(header,strict=strict)
            except:
                raise RecordError,"Info construction failed on " +\
                    "record with header %s"%header
            try:
                for seq in alignment.Seqs:
                    _process_seq(seq, strict)
                structure = struct_constructor(structure)
                yield family_info, alignment, structure
            except Exception, e:
                raise RecordError,"Sequence construction failed on " +\
                    "record with reference %s"%(family_info.Refs)
        else:
            try:
                family_info = label_constructor(header,strict=strict)
                for seq in alignment.Seqs:
                    _process_seq(seq, strict)
                structure = struct_constructor(structure)
                yield family_info, alignment, structure
            except Exception, e:
                if verbose:
                    print Exception, e
                continue

def _process_seq(seq, strict):
    """Adds info to seq, and to Aligned object if seq is hidden."""
    if hasattr(seq, 'data'):
        real_seq = seq.data
    else:
        real_seq = seq
    seq.Info = NameToInfo(real_seq, strict=strict)
    if seq.Info and 'Name' in seq.Info:
        seq.Name = seq.Info.Name
    if seq is not real_seq:
        real_seq.Name = seq.Name
        real_seq.Info = seq.Info

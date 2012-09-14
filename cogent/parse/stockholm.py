#!/usr/bin/env python
"""Provides a parser for Stockholm format files.
"""
from string import strip
from cogent.parse.record import RecordError
from cogent.parse.record_finder import DelimitedRecordFinder
from cogent.parse.clustal import ClustalParser
from cogent.core.sequence import RnaSequence as Rna
from cogent.core.sequence import ProteinSequence as Protein
from cogent.core.moltype import BYTES
from cogent.core.info import Info
from cogent.struct.rna2d import WussStructure
from cogent.util.transform import trans_all,keep_chars
from cogent.core.alignment import Alignment, DataError, SequenceCollection
from collections import defaultdict

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def is_empty_or_html(line):
    """Return True for HTML line and empty (or whitespace only) line.

    line -- string

    The Stockholm adaptor that retrieves records inlcudes two HTML tags in
    the record. These lines need to be ignored in addition to empty lines. 
    """
    if line.startswith('<pre') or line.startswith('</pre'):
        return True
    return (not line) or line.isspace()

Sequence = BYTES.Sequence
StockholmFinder = DelimitedRecordFinder('//', ignore=is_empty_or_html)

def load_from_clustal(data, seq_constructor=Sequence, strict=True,gap_char='-'):
    recs=[(name, seq_constructor(seq.replace('.',gap_char), )) for name, seq in\
        ClustalParser(data, strict)]
    lengths = [len(i[1]) for i in recs]
    if lengths and max(lengths) == min(lengths):
        return Alignment(recs, MolType=BYTES)
    else:
        return SequenceCollection(recs, MolType=BYTES)

_gf_field_names = {'AC':'AccessionNumber',\
                'ID':'Identification',\
                'DE':'Description',\
                'AU':'Author',\
                'SE':'AlignmentSource',\
                'SS':'StructureSource',\
                'BM':'BuildMethod',\
                'SM':'SearchMethod',\
                'GA':'GatheringThreshold',\
                'TC':'TrustedCutoff',\
                'NC':'NoiseCutoff',\
                'TP':'FamilyType',\
                'SQ':'Sequences',\
                'DC':'DatabaseComment',\
                'DR':'DatabaseReference',\
                'RC':None,\
                'RN':'ReferenceNumber',\
                'RM':'MedlineRef',\
                'RT':None,\
                'RA':None,\
                'RL':None,\
                'PI':'PreviousIdentifications',\
                'KW':'Keywords',\
                'CC':'Comment',\
                'NE':'PfamAccession',\
                'NL':'Location',\
                'WK':'WikipediaLink',\
                'CL':'Clan',\
                'MB':'Membership',\
                'NH':'NewHampshire',\
                'TN':'TreeID',
                'FT':'Feature'}

_gs_field_names = {'AC':'AccessionNumber',\
                   'DE':'Description',\
                   'DR':'DatabaseReference',\
                   'OS':'Organism',\
                   'OC':'OrganismClassification',\
                   'BP':'BasePair'}
                   
_gr_field_names = {'SS':'SecondaryStructure',\
                   'SA':'SurfaceAccessibility',\
                   'TM':'TransMembrane',\
                   'PP':'PosteriorProbability',\
                   'LI':'LigandBinding',\
                   'AS':'ActiveSite',\
                   'pAS':'ASPfam',\
                   'sAS':'ASSwissProt',\
                   'IN':'Intron',\
                   'RF':'ReferenceAnnotation'}

_gc_field_names = {'SS_cons':'ConsensusSecondaryStructure',\
                   'SA':'SurfaceAccessibility',\
                   'TM':'TransMembrane',\
                   'PP':'PosteriorProbability',\
                   'LI':'LigandBinding',\
                   'AS':'ActiveSite',\
                   'pAS':'ASPfam',\
                   'sAS':'ASSwissProt',\
                   'IN':'Intron',\
                   'RF':'ReferenceAnnotation'}



def GfToInfo(gf_lines,strict=True):
    """Returns a dict constructed from the GF lines.

    gf_lines is a list of lines that contain per-file annotation.
    Fields that can occur multiple times in a header are stored in a list.
    Fields that (should) occur only once are stored as a single value
    Comments are joined by ' ' to one field.
    Fields concerning the references are ignored, except for MedLine ID.
    """
    # construct temporary dictionary containing all original information
    initial_info = {}
    for line in gf_lines:
        line = line.strip()
        if not line:
            continue
        try:
            init,feature,content = line.split(None,2)
            if not init == '#=GF':
                raise RecordError
        except:
            if strict:
                raise RecordError, "Failed to extract feature and content " +\
                    "information from line %s"%(line)
            else:
                continue
        if feature in ['BM','DR','RM','CC','FT']:
            if feature in initial_info:
                initial_info[feature].append(content.strip())
            else:
                initial_info[feature] = [content.strip()]
        else:
            initial_info[feature] = content.strip()
            
    # transform initial dict into final one
    # throw away useless information; group information
    final_info={}
    for key in initial_info.keys():
        name = _gf_field_names.get(key,key)
        if name == 'Comment':
            value = ' '.join(initial_info[key])
        else:
            value = initial_info[key]
        final_info[name] = value
    
    return final_info

def GcToInfo(gc_lines,strict=True):
    """Returns a dict constructed from the GC lines.

    gc_lines is a list of lines that contain per column annotation.
    Fields that (should) occur only once are stored as a single value
    """
    # construct temporary dictionary containing all original information
    initial_info = defaultdict(list)
    for line in gc_lines:
        line = line.strip()
        if not line:
            continue
        try:
            init,feature,content = line.split(None,2)
            if not init == '#=GC':
                raise RecordError
        except:
            if strict:
                raise RecordError, "Failed to extract feature and content " +\
                    "information from line %s"%(line)
            else:
                continue


        initial_info[feature].append(content.strip())
    # transform initial dict into final one
    # throw away useless information; group information
    final_info={}
    for key in initial_info.keys():
        name = _gc_field_names.get(key,key)
        value = initial_info[key]
        final_info[name] = ''.join(value)
    
    return final_info

def GsToInfo(gs_lines,strict=True):
    """Returns a dict constructed from the GS lines.

    gs_lines is a list of lines that contain per-sequence annotation.
    Fields that can occur multiple times in a header are stored in a list.
    Fields that (should) occur only once are stored as a single value
    """
    # construct temporary dictionary containing all original information
    initial_info = {}
    for line in gs_lines:
        line = line.strip()
        if not line:
            continue
        try:
            init,seqname,feature,content = line.split(None,3)
            if not init == '#=GS':
                raise RecordError
        except:
            if strict:
                raise RecordError, "Failed to extract feature and content " +\
                    "information from line %s"%(line)
            else:
                continue
        if feature in ['DE','DR','BP']:
            if feature in initial_info:
                initial_info[feature][seqname].append(content.strip())
            else:
                initial_info[feature]= {seqname:[content.strip()]}
        elif feature not in initial_info:
            initial_info[feature]= {seqname:content.strip()}
        else:
            initial_info[feature][seqname]=content.strip()
            
    # transform initial dict into final one
    # throw away useless information; group information
    final_info={}
    for key in initial_info.keys():
        name = _gs_field_names.get(key,key)
        value = initial_info[key]
        final_info[name] = value
    
    return final_info

def GrToInfo(gr_lines,strict=True):
    """Returns a dict constructed from the GR lines.

    gr_lines is a list of lines that contain per-sequence AND per-Column
    annotation.
    Fields that can occur multiple times in a header are stored in a list.
    Fields that (should) occur only once are stored as a single value
    """
    # construct temporary dictionary containing all original information
    initial_info = defaultdict(dict)
    for line in gr_lines:
        line = line.strip()
        if not line:
            continue
        try:
            init,seqname,feature,content = line.split(None,3)
            if not init == '#=GR':
                raise RecordError
        except:
            if strict:
                raise RecordError, "Failed to extract feature and content " +\
                    "information from line %s"%(line)
            else:
                continue
        if feature not in initial_info:
            initial_info[feature][seqname]=[]
        elif seqname not in initial_info[feature]:
            initial_info[feature][seqname]=[]
        initial_info[feature][seqname].append(content.strip())
            
    # transform initial dict into final one
    # throw away useless information; group information
    final_info={}
    for feature in initial_info.keys():
        name = _gr_field_names.get(feature,feature)
        value = initial_info[feature]
        for k,v in value.items():
            value[k]=''.join(v)
        final_info[name] = value
    
    return final_info

AllToInfo = {'GF':GfToInfo,'GC':GcToInfo,'GS':GsToInfo,'GR':GrToInfo}

def is_gf_line(line):
    """Returns True if line is a GF line"""
    return line.startswith('#=GF')

def is_gc_line(line):
    """Returns True if line is a GC line"""
    return line.startswith('#=GC')

def is_gs_line(line):
    """Returns True if line is a GS line"""
    return line.startswith('#=GS')

def is_gr_line(line):
    """Returns True if line is a GR line"""
    return line.startswith('#=GR')

def is_seq_line(line):
    """Returns True if line is a sequence line"""
    return bool(line) and (not line[0].isspace()) and \
    (not line.startswith('#')) and (not line.startswith('//'))

def is_structure_line(line):
    """Returns True if line is a structure line"""
    return line.startswith('#=GC SS_cons ')

def MinimalStockholmParser(infile,strict=True,seq_constructor=Rna):
    """Yield successive records as (gf, gc, gs, gr, sequences, structure).
    
    gf is a list of GF lines
    gc is a list of GC lines
    gs is a list of GS lines
    gr is a list of GR lines
    sequences is an Alignment object. Sequences are Rna objects keyed by the
        original labels in the database.
    structure is a WussStructure
    """
    for record in StockholmFinder(infile):
        gf = []
        gc = []
        gs = []
        gr = []
        sequences = []
        structure = []
        for line in record:
            if is_gf_line(line):
                gf.append(line.strip())
            elif is_gc_line(line):
                gc.append(line.strip())
                if is_structure_line(line):
                    structure.append(line)
            elif is_gs_line(line):
                gs.append(line.strip())
            elif is_gr_line(line):
                gr.append(line.strip())
            elif is_seq_line(line):
                sequences.append(line)
            
            else:
                continue
        #sequence and structure are required. 
        #for example when looking at the stockholm format of just one family
        if not sequences:
            if strict:
                error = 'Found record with missing element(s): '
                if not sequences:
                    error += 'sequences'
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
        if structure:
            try:
                res = load_from_clustal(structure, strict=strict, gap_char='.')
                assert len(res.NamedSeqs) == 1 #otherwise multiple keys
                structure = res.NamedSeqs['#=GC SS_cons']
            except (RecordError, KeyError, AssertionError), e:
                if strict:
                    raise RecordError,\
                        "Can't parse structure of family"
                structure = None
        yield {'GF':gf, 'GC':gc, 'GS':gs, 'GR':gr}, sequences, structure
                
def StockholmParser(lines, seq_constructor=Rna, info_constructor_dict=\
    AllToInfo,struct_constructor=WussStructure,strict=True):
    """Yields (family_info, sequences, structure).

    Treats lines as a stream of Stockholm records.
    Family_info is the general information about the alignment.
    Sequences is an Alignment object. Each sequence has its own Info
        object with Genbank ID etc. Sequences are keyed by the original 
        label in the database.
    Structure is the consensus structure of the alignment, in Wuss format
    """
    for annotation, alignment, structure in MinimalStockholmParser\
        (lines,strict=strict,seq_constructor=seq_constructor):
        family_info = {}
        if strict:
            for k,v in annotation.items():
                label_constructor = info_constructor_dict[k]
                try:
                    family_info[k] = label_constructor(v,strict=strict)
                except:
                    raise RecordError,"Info construction failed on " +\
                        "record on the %s annotation"%(k)
            try:
                for seq in alignment.Seqs:
                    _process_seq(seq, strict)
                structure = struct_constructor(structure)
                alignment.Info.update(family_info)
                alignment.Info.update({'Struct':structure})
                yield alignment
            except Exception, e:
                raise RecordError,"Sequence construction failed on " +\
                    "record with reference %s"%\
                        (family_info['GF'].get('AccessionNumber',None))
        else:
            try:
                for k,v in annotation.items():
                    label_constructor = info_constructor_dict[k]
                    family_info[k] = label_constructor(v,strict=strict)

                for seq in alignment.Seqs:
                    _process_seq(seq, strict)
                structure = struct_constructor(structure)
                alignment.Info.update(family_info)
                alignment.Info.update({'Struct':structure})
                yield alignment
            except Exception, e:
                continue


def _process_seq(seq, strict):
    """Adds info to seq, and to Aligned object if seq is hidden."""
    if hasattr(seq, 'data'):
        real_seq = seq.data
    else:
        real_seq = seq

    if seq.Info and 'Name' in seq.Info:
        seq.Name = seq.Info.Name
    if seq is not real_seq:
        real_seq.Name = seq.Name
        real_seq.Info = seq.Info

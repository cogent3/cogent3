#/usr/bin/env python
"""Parsers for FASTA and related formats.
"""
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.record import RecordError
from cogent.core.info import Info, DbRef
from cogent.core.moltype import BYTES, ASCII

from string import strip
import cogent
import re

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight","Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

Sequence = BYTES.Sequence

def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')

def is_gde_label(x):
    """Checks if x looks like a GDE label line."""
    return x and x[0] in '%#'

def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()

def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()

FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)

def MinimalFastaParser(infile, strict=True, \
    label_to_name=str, finder=FastaFinder, \
    is_label=None, label_characters='>'):
    """Yields successive sequences from infile as (label, seq) tuples.

    If strict is True (default), raises RecordError when label or seq missing.
    """
    
    for rec in finder(infile):
        #first line must be a label line
        if not rec[0][0] in label_characters:
            if strict:
                raise RecordError, "Found Fasta record without label line: %s"%\
                    rec
            else:
                continue
        #record must have at least one sequence
        if len(rec) < 2:
            if strict:
                raise RecordError, "Found label line without sequences: %s" % \
                    rec
            else:
                continue
            
        label = rec[0][1:].strip()
        label = label_to_name(label)
        seq = ''.join(rec[1:])

        yield label, seq

GdeFinder = LabeledRecordFinder(is_gde_label, ignore=is_blank) 

def MinimalGdeParser(infile, strict=True, label_to_name=str):
    return MinimalFastaParser(infile, strict, label_to_name, finder=GdeFinder,\
        label_characters='%#')

def xmfa_label_to_name(line):
    (loc, strand, contig) = line.split()
    (sp, loc) = loc.split(':')
    (lo, hi) = [int(x) for x in loc.split('-')]
    if strand == '-':
        (lo, hi) = (hi, lo)
    else:
        assert strand == '+'
    name = '%s:%s:%s-%s' % (sp, contig, lo, hi)
    return name
   
def is_xmfa_blank_or_comment(x):
    """Checks if x is blank or an XMFA comment line."""
    return (not x) or x.startswith('=') or x.isspace()

XmfaFinder = LabeledRecordFinder(is_fasta_label, \
    ignore=is_xmfa_blank_or_comment)

def MinimalXmfaParser(infile, strict=True):
    # Fasta-like but with header info like ">1:10-1000 + chr1"
    return MinimalFastaParser(infile, strict, label_to_name=xmfa_label_to_name,
        finder=XmfaFinder)

def MinimalInfo(label):
    """Minimal info data maker: returns Name, and empty dict for info{}."""
    return label, {}

def NameLabelInfo(label):
    """Returns name as label split on whitespace, and Label in Info."""
    return label.split()[0], {'Label':label}

def FastaParser(infile,seq_maker=None,info_maker=MinimalInfo,strict=True):
    """Yields successive sequences from infile as (name, sequence) tuples.

    Constructs the sequence using seq_maker(seq, info=Info(info_maker(label))).

    If strict is True (default), raises RecordError when label or seq missing.
    Also raises RecordError if seq_maker fails.

    It is info_maker's responsibility to raise the appropriate RecordError or
    FieldError on failure.

    Result of info_maker need not actually be an info object, but can just be
    a dict or other data that Info can use in its constructor.
    """
    if seq_maker is None:
        seq_maker = Sequence
    for label, seq in MinimalFastaParser(infile, strict=strict):
        if strict:
            #need to do error checking when constructing info and sequence
            try:
                name, info = info_maker(label) #will raise exception if bad
                yield name, seq_maker(seq, Name=name, Info=info)
            except Exception, e:
                raise RecordError, \
                "Sequence construction failed on record with label %s" % label
        else:
            #not strict: just skip any record that raises an exception
            try:
                name, info = info_maker(label)
                yield(name, seq_maker(seq, Name=name, Info=info))
            except Exception, e:
                continue

#labeled fields in the NCBI FASTA records
NcbiLabels = {
'dbj':'DDBJ',
'emb':'EMBL',
'gb':'GenBank',
'ref':'RefSeq',
}

def NcbiFastaLabelParser(line):
    """Creates an Info object and populates it with the line contents.
    
    As of 11/12/03, all records in genpept.fsa and the human RefSeq fasta
    files were consistent with this format.
    """
    info = Info()
    try:
        ignore, gi, db, db_ref, description = map(strip, line.split('|', 4))
    except ValueError:  #probably got wrong value
        raise RecordError, "Unable to parse label line %s" % line
    info.GI = gi
    info[NcbiLabels[db]] = db_ref
    info.Description = description
    return gi, info

def NcbiFastaParser(infile, seq_maker=None, strict=True):
    return FastaParser(infile, seq_maker=seq_maker, 
        info_maker=NcbiFastaLabelParser, strict=strict)

class RichLabel(str):
    """Object for overloaded Fasta labels. Holds an Info object storing keyed
    attributes from the fasta label. The str is created from a provided format
    template that uses the keys from the Info object."""
    
    def __new__(cls, info, template="%s"):
        """Arguments:
        
            - info: a cogent.core.info.Info instance
            - template: a string template, using a subset of the keys in info.
              Defaults to just '%s'.
        
        Example:
            label = RichLabel(Info(name='rat', species='Rattus norvegicus'),
                        '%(name)s')"""
        label = template % info
        new = str.__new__(cls, label)
        new.Info = info
        return new
    

def LabelParser(display_template, field_formatters, split_with=":", DEBUG=False):
    """returns a function for creating a RichLabel's from a string
    
    Arguments;
        - display_template: string format template
        - field_formatters: series of 
                (field index, field name, coverter function)
        - split_with: characters separating fields in the label.
          The display_template must use at least one of the assigned field
          names."""
    indexed = False
    for index, field, converter in field_formatters:
        if field in display_template:
            indexed = True
    assert indexed, "display_template [%s] does not use a field name"\
                    % display_template
    sep = re.compile("[%s]" % split_with)
    def call(label):
        label = [label, label[1:]][label[0] == ">"]
        label = sep.split(label)
        if DEBUG:
            print label
        info = Info()
        for index, name, converter in field_formatters:
            if callable(converter):
                try:
                    info[name] = converter(label[index])
                except IndexError:
                    print label, index, name
                    raise
            else:
                info[name] = label[index]
        return RichLabel(info, display_template)
    return call

def GroupFastaParser(data, label_to_name, group_key="Group", aligned=False,
        moltype=ASCII, done_groups=None, DEBUG=False):
    """yields related sequences as a separate seq collection
    
    Arguments:
        - data: line iterable data source
        - label_to_name: LabelParser callback
        - group_key: name of group key in RichLabel.Info object
        - aligned: whether sequences are to be considered aligned
        - moltype: default is ASCII
        - done_groups: series of group keys to be excluded
        """
    
    done_groups = [[], done_groups][done_groups is not None]
    parser = MinimalFastaParser(data, label_to_name=label_to_name, finder=XmfaFinder)
    group_ids = []
    current_collection = {}
    for label, seq in parser:
        seq = moltype.makeSequence(seq, Name=label, Info=label.Info)
        if DEBUG:
            print "str(label) ",str(label), "repr(label)", repr(label)
        if not group_ids or label.Info[group_key] in group_ids:
            current_collection[label] = seq
            if not group_ids:
                group_ids.append(label.Info[group_key])
        else:
            # we finish off check of current before creating a collection
            if group_ids[-1] not in done_groups:
                info = Info(Group=group_ids[-1])
                if DEBUG:
                    print "GroupParser collection keys", current_collection.keys()
                seqs = cogent.LoadSeqs(data=current_collection, moltype=moltype,
                                aligned=aligned)
                seqs.Info = info
                yield seqs
            current_collection = {label: seq}
            group_ids.append(label.Info[group_key])
    info = Info(Group=group_ids[-1])
    seqs = cogent.LoadSeqs(data=current_collection, moltype=moltype,
                    aligned=aligned)
    seqs.Info = info
    yield seqs
    

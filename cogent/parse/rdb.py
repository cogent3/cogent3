#!/usr/bin/env python
"""Provides a parser for Rdb format files.

Data in from the European rRNA database in distribution format.
"""
from string import strip, maketrans
from cogent.parse.record_finder import DelimitedRecordFinder
from cogent.parse.record import RecordError
from cogent.core.sequence import Sequence, RnaSequence
from cogent.core.info import Info
from cogent.core.alphabet import AlphabetError

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

RdbFinder = DelimitedRecordFinder('//')

_field_names = {'acc':'rRNA',\
                'src':'Source',\
                'str':'Strain',\
                'ta1':'Taxonomy1',\
                'ta2':'Taxonomy2',\
                'ta3':'Taxonomy3',\
                'ta4':'Taxonomy4',\
                'chg':'Changes',\
                'rem':'Remarks',\
                'aut':'Authors',\
                'ttl':'Title',\
                'jou':'Journal',\
                'dat':'JournalYear',\
                'vol':'JournalVolume',\
                'pgs':'JournalPages',\
                'mty':'Gene',\
                'del':'Deletions',\
                'seq':'Species'}


def InfoMaker(header_lines):
    """Returns an Info object constructed from the headerLines."""
    info = Info()
    for line in header_lines:
        all = line.strip().split(':',1)
        #strip out empty lines, lines without name, lines without colon
        if not all[0] or len(all) != 2: 
            continue
        try:
            name = _field_names[all[0]]
        except KeyError:
            name = all[0]
            
        value = all[1].strip()
        info[name] = value
    return info

def is_seq_label(x):
    "Check if x looks like a sequence label line."""
    return x.startswith('seq:')

def MinimalRdbParser(infile,strict=True):
    """Yield successive sequences as (headerLines, sequence) tuples.
    
    If strict is True (default) raises RecordError when 'seq' label is missing
    and if the record doesn't contain any sequences.
    """
    for rec in RdbFinder(infile):
        index = None
        for line in rec:
            if is_seq_label(line):
                index = rec.index(line) + 1 #index of first sequence line
        
        # if there is no line that starts with 'seq:' throw error or skip
        if not index:
            if strict:
                raise RecordError, "Found Rdb record without seq label "\
                    + "line: %s"%rec[0]
            else:
                continue
            
        headerLines = rec[:index]
        sequence = ''.join(rec[index:-1]) #strip off the delimiter
        if sequence.endswith('*'):
            sequence = sequence[:-1] #strip off '*'

        #if there are no sequences throw error or skip
        if not sequence:
            if strict:
                raise RecordError, "Found Rdb record without sequences: %s"\
                    %rec[0]
            else:
                continue

        yield headerLines, sequence

def create_acceptable_sequence(sequence):
    """Return clean sequence as string, 'o' -> '?', sec. structure deleted

    sequence: string of characters
    SeqConstructor: constructor function for sequence creation

    Will replace 'o' by '?'.
    Will strip out secondary structure annotation.
    """
    t = maketrans('o','?')
    # strip out secondary structure annotation {}[]()^
    return sequence.translate(t, "{}[]()^") #should be accepted by RnaSequence
    
    
def RdbParser(lines, SeqConstructor=RnaSequence, LabelConstructor=InfoMaker, \
    strict=True):
    """Yield sequences from the Rdb record.

    lines: a stream of Rdb records.
    SeqConstructor: constructor function to create the final sequence object
    LabelConstructor: function that creates Info dictionary from label lines
    strict: boolean, when True, an error is raised when one occurs, when False,
        the record is ignored when an error occurs.

    This function returns proper RnaSequence objects when possible. It strips
    out the secondary structure information, and it replaces 'o' by '?'. The
    original sequence is stored in the info dictionary under 'OriginalSeq'.
    If the original sequence is the desired end product, use MinimalRdbParser.
    """
    for header, sequence in MinimalRdbParser(lines,strict=strict):
        info = LabelConstructor(header)
        clean_seq = create_acceptable_sequence(sequence)
        # add original raw sequence to info
        info['OriginalSeq'] = sequence
        if strict:
            #need to do error checking while constructing info and sequence
            try:
                yield SeqConstructor(clean_seq, Info = info)
            except AlphabetError:
                raise RecordError(\
                "Sequence construction failed on record with reference %s."\
                %(info.Refs))
        else:
            #not strict: just skip any record that raises an exception
            try:
                yield SeqConstructor(clean_seq, Info=info)
            except:
                continue

if __name__ == '__main__':
    from sys import argv
    filename = argv[1]
    for sequence in RdbParser(open(filename)):
        print sequence.Info.Species
        print sequence

#!/usr/bin/env python
"""Parser for 454 Flowgram files"""

__author__ = "Jens Reeder, Julia Goodrich"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jens Reeder","Julia Goodrich", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jens Reeder"
__email__ = "jreeder@colorado.edu"
__status__ = "Development"

from string import strip
from random import sample
from itertools import izip

from cogent.parse.flowgram import Flowgram
from cogent.parse.record_finder import LabeledRecordFinder, is_fasta_label,\
     DelimitedRecordFinder, is_empty


def get_header_info(lines):
    """Returns the Information stored in the common header as a dictionary

    lines can be a a list or a file handle
    """
    header_dict = {}

    for line in lines:
        if line.startswith('Common Header'):
            continue
        if is_empty(line):
            break
       
        key, value = line.strip().split(':')
        header_dict[key] = value.strip()
    
    if isinstance(lines, file):
        lines.seek(0)

    return header_dict

def get_summaries(handle, number_list = None, name_list=None, all_sums = False):
    """Returns specified flowgrams and sequence summaries as generator
    handle can be a list of lines or a file handle
    number_list is a list of the summaries wanted by their index in the sff
        file, starts at 0
    name_list is a list of the summaries wanted by their name in the sff file
    all_sums if true will yield all the summaries in the order they appear in
        the file

    One and only one of the parameters must be set
    """
    sff_info = LabeledRecordFinder(is_fasta_label,constructor=strip)
    sum_gen = sff_info(handle)

    if number_list:
        assert not (name_list or all_sums)
        num = len(number_list)
        for i,s in enumerate(sum_gen):
            if i-1 in number_list:
                yield s
                num -= 1
            if num == 0:
                break
            
    elif name_list:
        assert not all_sums
        for s in sum_gen:
            if s[0].strip('>') in name_list:
                yield s

    elif all_sums:
        header = True
        for s in sum_gen:
            if header:
                header = False
                continue
            yield s
    else:
        raise ValueError, "number_list, name_list or all_sums must be specified"


def get_all_summaries(lines):
    """Returns all the flowgrams and sequence summaries in list of lists"""
    sff_info = LabeledRecordFinder(is_fasta_label,constructor=strip)

    return list(sff_info(lines))[1::]

def split_summary(summary):
    """Returns dictionary of one summary"""
    summary_dict = {}

    summary_dict["Name"] = summary[0].strip('>')
    for line in summary[1::]:
        key, value = line.strip().split(':')
        summary_dict[key] = value.strip()
        
    return summary_dict

def parse_sff(lines):
    """Creates list of flowgram objects from a SFF file
    """
    head = get_header_info(lines)
    summaries = get_all_summaries(lines)

    flows = []
    for s in summaries:
        t = split_summary(s)
        flowgram = t["Flowgram"]
        del t["Flowgram"]
        flows.append(Flowgram(flowgram, Name = t["Name"],
                              floworder =head["Flow Chars"], header_info = t))
    return flows, head


def lazy_parse_sff_handle(handle):
    """Returns one flowgram at a time 
    """
    sff_info = LabeledRecordFinder(is_fasta_label,constructor=strip)
    sff_gen = sff_info(handle)

    header_lines = sff_gen.next()
    header = get_header_info(header_lines)

    return (_sff_parser(sff_gen, header), header)

def _sff_parser(handle, header):
    for s in handle:
        t = split_summary(s)
        flowgram = t["Flowgram"]
        del t["Flowgram"]
        flowgram = Flowgram(flowgram, Name = t["Name"],
                            KeySeq=header["Key Sequence"],
                            floworder = header["Flow Chars"],
                            header_info = t)
        
        yield flowgram

def get_random_flows_from_sff(filename, num=100, size=None):
    """Reads size many flows from filename and return sample of num randoms.
    
    Note: size has to be the exact number of flowgrams in the file, otherwise 
    the result won't be random or less than num flowgrams will be returned

    filename: sff.txt input file

    num: number of flowgrams in returned sample

    size: number of flowgrams to sample from 
    """

    if(size==None):
        size = count_sff(open(filename))
    if (size<num):
        size = num
    
    (flowgrams, header) =  lazy_parse_sff_handle(open(filename))
    idxs = sample(xrange(size), num)
    idxs.sort()
    i = 0   
    for (j,f) in izip(xrange(size), flowgrams):
        if (idxs[i] == j):
            i += 1
            yield f
            if (i>=num):
                break

def count_sff(sff_fh):
    """Counts flowgrams in a sff file"""
    
    (flowgrams, header) = lazy_parse_sff_handle(sff_fh)
    i=0
    for f in flowgrams:
        i+=1
    return i


def sff_to_fasta(sff_fp, out_fp):
    """Transform an sff file to fasta"""
    (flowgrams, header) = lazy_parse_sff_handle(open(sff_fp))

    out_fh = open(out_fp, "w")
                                     
    for f in flowgrams:
        out_fh.write(f.toFasta()+"\n")


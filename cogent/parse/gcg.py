#!/usr/bin/env python

import logging
LOG = logging.getLogger('cogent.input')

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@anu.edu.au"
__status__ = "Production"

def MsfParser(f):
    """Read sequences from a msf format file"""
    alignmentdict = {}
    #parse optional header
    #parse optional text information
    #file header and sequence header are seperated by a line ending in '..'
    line = f.readline().strip()
    for line in f:
        line = line.strip()
        if line.endswith('..'):
            break
    #parse sequence info
    seqinfo = {}
    for line in f:
        line = line.strip()
        if line.startswith('//'):
            break
        line = line.split()
        if line and line[0] == 'Name:':
            seqinfo[line[1]] = int(line[3])
    #parse sequences
    sequences = {}
    for line in f:
        line = line.strip().split()
        if line and sequences.has_key(line[0]):
            sequences[line[0]] += ''.join(line[1:])
        elif line and seqinfo.has_key(line[0]):
            sequences[line[0]] = ''.join(line[1:])
    #consistency check
    if len(sequences) != len(seqinfo):
        LOG.warning("Number of loaded seqs[%s] not same as "\
            "expected[%s]." % (len(sequences), len(seqinfo)))
    for name in sequences:
        if len(sequences[name]) != seqinfo[name]:
            LOG.warning("Length of loaded seqs [%s] is [%s] not "\
            "[%s] as expected." % (name,len(sequences[name]),seqinfo[name]))
    
    #yield sequences
    for name in sequences:
        yield (name, sequences[name])

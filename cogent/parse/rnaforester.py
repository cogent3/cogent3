#!/usr/bin/env python

from string import strip,split
from cogent.struct.rna2d import Pairs,ViennaStructure

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def rnaforester_parser(lines):
    """Parser for RNAforester output format

    Returns a list containing: alignemnt,consensus sequence and consensus 
    structure
    Ex: [{alignment},consensus sequence,[consensus structure]]
    """
    
    result = []
    
    for block in cluster_parser(lines):
        for struct in line_parser(block):
            result.append(struct)

    return result


def cluster_parser(lines): 
    """To parse lines into rnaforester clluster blocks"""
    block = []
    first = True
    record = False
    for line in lines:
        if line.startswith('RNA Structure Cluster Nr:'): #new cluster block
            record = True
            if not first:
                yield block
                block = []
            first = False
        if record:
            block.append(line)
    yield block

def line_parser(block):
    """Parses the stdout output from RNAforester and return the concensus 
    structure prediction along with alignment and consensus sequence.
    """
    odd = True
    record = False
    first = True
    seq = ''
    con_seq = ''
    struct = ''
    alignment = {}
    for line in block:
        #find alignments
        if line.startswith('seq'):
            if line.__contains__(')') or line.__contains__('('):
                continue
            else:
                sline = line.strip().split()
                name = sline[0]
                tmp_seq = sline[-1]
                if alignment.__contains__(name):
                    seq = alignment[name]
                    seq = ''.join([seq,tmp_seq])
                    alignment[name] = seq
                else:
                    alignment[name] = tmp_seq

        if line.startswith('Consensus sequence/structure:'): #start
            record = True
            if not first:
                struct = to_pairs(struct)
                yield [alignment,con_seq,struct]
                result = []
            first = True
        elif record:
            if line.startswith('                         '):
                line = line.strip()
                if odd:
                    con_seq = ''.join([con_seq,line])
                    odd = False
                else:
                    struct = ''.join([struct,line])
                    odd = True

    struct = to_pairs(struct)
    yield [alignment,con_seq,struct]


def to_pairs(struct):
    """
    Takes a structure in dot-bracket notation and converts it to pairs notation

    functions tested in rna2d.py
    """
    struct = ViennaStructure(struct)
    struct = struct.toPairs()

    return struct


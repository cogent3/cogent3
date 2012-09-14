#!/usr/bin/env python
#file cogent/parse/gibbs.py
"""Parses Gibbs Sampler output file and creates MotifResults object."""

from cogent.parse.record_finder import LabeledRecordFinder
from cogent.motif.util import Location, ModuleInstance, Module, Motif,\
     MotifResults
from cogent.core.moltype import DNA, RNA, PROTEIN
from numpy import exp

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def get_sequence_and_motif_blocks(lines):
    """Returns main block of data as a list.
    """
    gibbs_map_maximization = \
        LabeledRecordFinder(lambda x: 'MAP MAXIMIZATION RESULTS' in x)
    seq_block, motif_block = list(gibbs_map_maximization(lines))
    return seq_block, motif_block

def get_sequence_map(lines):
    """Returns dict mapping Gibbs sequence number to sequence ID.
    
        - ex: sequence numbers mapping to gis:
            {'1':'1091044',
             '2':'11467494',
             '3':'11499727'}
    """
    sequence_map = {}
    sequence_finder = \
        LabeledRecordFinder(lambda x: x.startswith('Sequences to be Searched:'))
    sequence_block = list(sequence_finder(lines))[-1]
    for i in sequence_block[2:]:
        if i.startswith('#'):
            num,label = i.strip().split(' ',1)
            num = num.strip()
            label = label.strip()
            sequence_map[num[1:]] = label
        else:
            break
    return sequence_map
    

def get_motif_blocks(lines):
    """Returns list of motif blocks given main block as lines.
    """
    gibbs_motif = LabeledRecordFinder(lambda x: 'MOTIF' in x)
    return list(gibbs_motif(lines))[1:]
    
def get_motif_sequences(lines):
    """Returns list of tuples with motif sequence information given motif block.
    
        - result is list of tuples :
            [(seq_num, motif_start, motif_seq, motif_sig),]
    """
    motif_list = []
    motif_seq_finder = LabeledRecordFinder(lambda x: 'columns' in x)
    motifs = list(motif_seq_finder(lines))[-1]
    for m in motifs[2:]:
        if ',' in m:
            curr = m.strip().split()
            motif_num = curr[1]
            seq_num = curr[0].split(',')[0]
            motif_start = int(curr[2])-1
            #If motif does not start at beginning of sequence:
            if motif_start > 0:
                motif_seq = curr[4]
            #Motif starts at beginning of sequence, no context before motif
            else:
                motif_seq = curr[3]
            motif_sig = float(curr[-3])
            motif_list.append((seq_num,motif_start,motif_seq,motif_sig,motif_num))
        else:
            break
    return motif_list

def get_motif_p_value(lines):
    """Returns the motif p-value given motif block.
    """
    motif_p_finder = LabeledRecordFinder(lambda x: x.startswith('Log Motif'))
    motif_p_block = list(motif_p_finder(lines))[-1]
    log_motif_portion = float(motif_p_block[0].split()[-1])
    return exp(log_motif_portion)

def guess_alphabet(motif_list):
    """Returns alphabet given a motif_list from get_motif_sequences.
    
        - temp hack...should really think of a better way to get alphabet,
        but Gibbs sampler help tells nothing of how it guesses the alphabet
        or what it does with degenerate characters.
        - only allows for 2 degenerate nucleotide chars.
    """
    alpha_dict = {}
    for motif in motif_list:
        for char in motif[2]:
            if char not in alpha_dict:
                alpha_dict[char]=0
            alpha_dict[char]+=1
    if len(alpha_dict) > 6:
        alphabet=PROTEIN
    elif 'T' in alpha_dict:
        alphabet=DNA
    else:
        alphabet=RNA
    return alphabet
    

def build_module_objects(motif_block, sequence_map, truncate_len=None):
    """Returns module object given a motif_block and sequence_map.
    
        - motif_block is list of lines resulting from calling get_motif_blocks
        - sequence_map is the mapping between Gibbs sequence numbering and 
        sequence id from fasta file.
    """
    #Get motif id
    motif_id = motif_block[0].strip().split()[-1]

    #Get motif_list
    motif_list = get_motif_sequences(motif_block)
    #Get motif p-value
    motif_p = get_motif_p_value(motif_block)
    #Guess alphabet from motif sequences
    alphabet = guess_alphabet(motif_list)
    
    #Create Module object(s)
    gibbs_module = {}

    module_keys = ["1"]

    for motif in motif_list:

        seq_id = str(sequence_map[motif[0]])

        if truncate_len:
            seq_id = seq_id[:truncate_len]
        
        start = motif[1]
        seq = motif[2]
        sig = motif[3]
        motif_num = "1" 
        
        #Create Location object
        location = Location(seq_id, start, start + len(seq))
        #Create ModuleInstance
        mod_instance = ModuleInstance(seq,location,sig)
        cur_key = (seq_id,start)
        
        gibbs_module[(seq_id,start)]=mod_instance
    gibbs_mod = Module(gibbs_module,MolType=alphabet)
    gibbs_mod.Pvalue = motif_p
    gibbs_mod.ID = motif_id + module_keys[0]
    yield gibbs_mod

def module_ids_to_int(modules):
    """Given a list of modules changes alpha ids to int ids.
    """
    id_map = {}
    for m in modules:
        id_map[m.ID]=True
    for i,mid in enumerate(sorted(id_map.keys())):
        id_map[mid]=str(i)
    for m in modules:
        m.ID=id_map[m.ID]


def GibbsParser(lines, truncate_len=None, strict=True):
    """Returns a MotifResults object given a Gibbs Sampler results file.
    
        - only works with results from command line version of Gibbs Sampler
    """
    try:
        #Get sequence and motif blocks
        sequence_block, motif_block = get_sequence_and_motif_blocks(lines)
    except Exception, e:
        if strict:
            raise e
        else:
            return None
    #Create MotifResults object
    gibbs_motif_results = MotifResults()

    #Get sequence map
    sequence_map = get_sequence_map(sequence_block)
    #Get motif blocks
    motif_blocks = get_motif_blocks(motif_block)
    #Get modules
    for module in motif_blocks:
        if module[1] == 'No Motifs Detected':
            print "No Modules detcted!!", module[0]
            continue
        for cur_smod in build_module_objects(module, sequence_map, truncate_len=truncate_len):
            gibbs_motif_results.Modules.append(cur_smod)
    module_ids_to_int(gibbs_motif_results.Modules)
    
    for ct, module in enumerate(gibbs_motif_results.Modules):
        gibbs_motif_results.Motifs.append(Motif(module))
    gibbs_motif_results.Alphabet=module.Alphabet
    return gibbs_motif_results

if __name__ == "__main__":
    from sys import argv, exit
    print "Running..."

    if len(argv) != 2:
        print "Usage: gibbs.py <gibbs out file>"
        exit(1)
    mr = GibbsParser(open(argv[1]), 24)
    print mr
    for module in mr.Modules:
        print module.ID
        print str(module)


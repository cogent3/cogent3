#!/usr/bin/env python
from __future__ import division
from cogent.app.muscle_v38 import muscle_seqs
from cogent.app.util import get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from numpy import array

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

three_to_one = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H','ILE':'I','LYS':'K','LEU':'L','FME':'M','MET':'M','MSE':'M',
    'ASN':'N','PRO':'P',
    'GLN':'Q','ARG':'R','SER':'S','THR':'T','SEC':'U','VAL':'V','TRP':'W',
    'TYR':'Y'}
nucleotides = {'A':'A','C':'C','G':'G','U':'U','T':'T','I':'I',\
    '1MA':'A','2MA':'A','SRA':'A','MIA':'A','MAD':'A','A2M':'A','AVC':'A',
    'PPU':'A','AET':'A','A23':'A','5AA':'A','T6A':'A','MTU':'A','LCA':'A',
    '5MC':'C','CCC':'C','OMC':'C','CH':'C','10C':'C',
    '1MG':'G','QUO':'G','G7M':'G','GDP':'G','YG':'G','7MG':'G','OMG':'G',
    'M2G':'G','2MG':'G','GTP':'G','YYG':'G',
    'PSU':'U','H2U':'U','5MU':'U','4SU':'U','ONE':'U','+U':'U','DHU':'U',
    'IU':'U','SSU':'U','OMU':'U','UR3':'U','MNU':'U','5BU':'U','S4U':'U',
    'UMS':'U'
    }


def get_aligned_muscle(seq1,seq2):
    """Returns aligned sequences and frac_same using MUSCLE.

        This needs to be moved to the muscle app controller
    """
    outname = get_tmp_filename()
    res = muscle_seqs([seq1,seq2], add_seq_names=True, WorkingDir="/tmp", out_filename=outname)
    seq1_aligned,seq2_aligned =list(MinimalFastaParser(res['MuscleOut'].read()))
    res.cleanUp()
    del(res)
    seq1_aligned = seq1_aligned[1][1:]
    seq2_aligned = seq2_aligned[1][1:]
    frac_same = (array(seq1_aligned, 'c') == array(seq2_aligned, 'c')).sum(0)\
            / min(len(seq1), len(seq2))
    
    return seq1_aligned,seq2_aligned,frac_same

def get_chains(lines):
    """From list of lines in pdb records, returns dict {chain:[(pos,residue)]}.

    Keeps original 1-based numbering. All residues will be returned.
    """
    chains = {}
    last_resnum = None
    ter=False
    prev_chains=[]
    for line in lines:
        #skip if not an atom line
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            #If chain terminated, record chain.
            if line.startswith('TER'):
                ter=True
                prev_chains.append(chain)
            else:
                continue
        else:
            residue = line[17:20].strip()
            chain = line[21].strip()
            try:
                resnum = int(line[22:27].strip())
            #Some resnum columns have non-integer values
            except ValueError:
                resnum = line[22:27].strip()
            #Continue until next chain is found.
            if ter:
                if chain in prev_chains:
                    continue
                else:
                    ter=False
                    if chain not in chains:
                        chains[chain] = []
                    curr_chain = chains[chain]
            else:
                if chain not in chains:
                    chains[chain] = []
                curr_chain = chains[chain]

                if resnum != last_resnum:
                    curr_chain.append((resnum,residue))
                last_resnum = resnum
    return chains



def ungapped_to_pdb_numbers(chain_list):
    """From a chain list, map seq position -> res number."""
    return dict(enumerate([i[0]for i in chain_list]))

def chains_to_seqs(chains):
    """Returns sequences as an array of chars for each chain.
    
    Will concatenate all the residues that exist; it's the job of alignment
    or similar to align it with other seqs. Can get numbering 
    """
    result = {}
    chain_to_seq_type = {}
    for chain_id, residues in chains.items():
        seq_type = None
        curr_seq = []
        
        first_res = residues[0][1].strip()
        if first_res in three_to_one:
            seq_type = 'Protein'
        else:
            seq_type = 'Nucleotide'
            
        for res_id, seq in residues:
            seq = seq.strip()
            if seq_type == 'Nucleotide':
                curr_seq.append(nucleotides.get(seq, 'N'))
            else:
                curr_seq.append(three_to_one.get(seq, '?'))

        result[chain_id] = ''.join(curr_seq)
        chain_to_seq_type[chain_id]=seq_type
    return result, chain_to_seq_type

def get_best_muscle_hits(subject_seq, query_aln,threshold,use_shorter=True):
    """Returns subset of query_aln with alignment scores above threshold.
    
        - subject_seq is sequence aligned against query_aln seqs.
        - query_aln is dict or Alignment object with candidate seqs to be
        aligned with subject_seq.
        - threshold is an alignment score (fraction shared aligned length)
        which returned seqs must be above when aligned w/ subject_seq.
        - use_shorter (default=True) is to decide whether to use the length
        of the shorter sequence to calculate the alignment score.
    """
    keep={}
    #best = 0
    for query_label, query_seq in query_aln.items():
        subject_aligned, query_aligned, frac_same = \
            get_aligned_muscle(subject_seq,query_seq)
        #if frac_same > best:
        if frac_same > threshold: 
            keep[query_label]=query_seq
            #best=frac_same
    return keep

def get_matching_chains(subject_seq, pdb_lines,\
    subject_type='Protein',threshold=0.8):
    """Returns PDB chains that match subject_seq.
    
        - subject_seq must be a sequence string.
        - pdb_lines must be a list of lines from a PDB record.
        - subject_type must be the type of sequence that subject_seq is.  This
        is used to build blast database.
    """
    #Get PDB sequence info
    chains = get_chains(pdb_lines)
    #Get pdb numbering
    ungapped_to_pdb = {}
    for k,v in chains.items():
        ungapped_to_pdb[k]=ungapped_to_pdb_numbers(v)
    #Get pdb alignment and sequence types
    pdb_aln, pdb_types = chains_to_seqs(chains)
    
    #Get only sequences of subject_type
    pdb_matching = \
        dict([(k,pdb_aln[k]) for k,v in pdb_types.items() if v == subject_type])
    
    #if there is more than one chain in pdb_matching
    if len(pdb_matching) > 1:
        #Get best hits using MUSCLE.
        pdb_matching = \
            get_best_muscle_hits(subject_seq, pdb_matching, threshold)
    
    return pdb_matching, ungapped_to_pdb

def align_subject_to_pdb(subject_seq, pdb_matching):
    """Returns pairwise aligned subject_seq and pdb_matching alignment.
    
        - result will be a dict:
            {pdb_chain:(aligned_subject, aligned_pdb)}
    """
    result = {}
    for pdb_chain, pdb_seq in pdb_matching.items():
        subject_aligned, pdb_aligned,frac_same = \
            get_aligned_muscle(subject_seq, pdb_seq)
        result[pdb_chain]=(subject_aligned,pdb_aligned)
    
    return result


#####The following code must be in the .pml script:####
def iterate_blocks(seq, max_len=100):
    """Yields successive blocks up to max_len from seq."""
    curr = 0
    while curr < len(seq):
        yield seq[curr:curr+max_len]
        curr += max_len

def make_color_list(colors, prefix="color_"):
    """Makes list of colors, sequentially numbered after prefix."""
    return [(prefix+str(i+1), color) for i, color in enumerate(colors)]

def set_color_list(color_list):
    """Uses cmd to set all the items in a list of colors as named colors."""
    for name, color in color_list:
        cmd.set_color(name, color)

def set_seq_colors(colors, indices, chain_id):
    """Takes list of colors same length as seq, index mapping, and chain id."""
    for i, color in enumerate(colors):
        idx = indices[i].split('+')
        for block in iterate_blocks(idx):
            curr = '+'.join(block)
            cmd.color(color, "chain %s and resi %s" % (chain_id, curr))

def set_show_shapes(indices, chain_id, shape="sticks"):
    """Takes list of indices and a chain id and sets to shape.
    """
    for block in iterate_blocks(indices):
        str_indices = '+'.join(map(str,block))
        cmd.show(shape,"chain %s and resi %s" % (chain_id, str_indices))
    

#pymol coloring functions string:
PYMOL_FUNCTION_STRING = \
'''
def iterate_blocks(seq, max_len=100):
    """Yields successive blocks up to max_len from seq."""
    curr = 0
    while curr < len(seq):
        yield seq[curr:curr+max_len]
        curr += max_len

def make_color_list(colors, prefix="color_"):
    """Makes list of colors, sequentially numbered after prefix."""
    return [(prefix+str(i+1), color) for i, color in enumerate(colors)]

def set_color_list(color_list):
    """Uses cmd to set all the items in a list of colors as named colors."""
    for name, color in color_list:
        cmd.set_color(name, color)

def set_seq_colors(colors, indices, chain_id):
    """Takes list of colors same length as seq, index mapping, and chain id."""
    for i, color in enumerate(colors):
        idx = indices[i].split('+')
        for block in iterate_blocks(idx):
            curr = '+'.join(block)
            cmd.color(color, "chain %s and resi %s" % (chain_id, curr))

def set_show_shapes(indices, chain_id, shape="sticks"):
    """Takes list of indices and a chain id and sets to shape.
    """
    str_indices = "+".join(map(str,indices))
    cmd.show(shape,"chain %s and resi %s" % (chain_id, str_indices))

'''

MAIN_FUNCTION_STRING = \
'''
cmd.hide()
cmd.show("cartoon")
cmd.color("white")
    
'''
        

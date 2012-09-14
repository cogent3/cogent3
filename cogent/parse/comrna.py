#!/usr/bin/env python
#file comrna_parser.py
"""Parser for comRNA output format

To reduce number of structures that parser report set first to True, then parser
will only report structures from first block (Maximum stem similarity score block). 

The function common can be used to have the most common occuring structure be 
reported first, if all structure only occurs once the structure with the most pairs 
will be reported as the most common. 
"""
from cogent.util.transform import make_trans
from cogent.struct.rna2d   import Pairs
from cogent.struct.knots   import opt_single_random
from string                import index

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def comRNA_parser(lines=None,pseudo=True,first=False):
    """Parsed comRNA output.

    pseudo - if True, report results with pseudoknots; if flag
             is False pseudoknots will be removed.
    """
    names = get_names(lines)   
    result = []
    
    for block in minimalComrnaParser(lines):
        for structure in blockParser(block):
            for struct in structParser(structure):
                for pairs,seq in pairsParser(struct,names):
                    result.append([seq,pairs])
                if first:    
                    break
            if first:
                break
        if first:
            break

    if not pseudo:
        tmp = []
        for block in result:
            tmp.append([block[0],opt_single_random(block[-1])])
        result = tmp
    return result

def common(structs):
    """
    Will return a list of sequences and structures with the most common structure 
    first in the list. (rest or list unordered!)
    Don't care which of the sequences for the "winning" sequence that is reported
    since the are not ranked amongst them self.
    """
    
    frequency = {}
    v = 0
    indx = 0
    result = []
    tmp_list = [] #lookup the seq for the structures,dont care which winner seq 
    key = []
    for block in structs:
        tmp_list.extend(block)
        p = tuple(block[-1])
        if frequency.__contains__(p): #everytime struct p appears count up by 1
            frequency[p]+=1
        else:
            frequency[p]=1
        nr = frequency[p]
        if nr > v: #Which struct appears most times
            v = nr
            key = p

    #if winning structure has frequency == 1 all structure apper only once
    if frequency[key]==1:
        longest = 0
        for block in structs:
            l = len(block[-1])
            if l > longest: #pick longest sequence as the winner
                key = tuple(block[-1])
    
    winner = Pairs(key)
    indx = tmp_list.index(winner)-1
    result.append([tmp_list[indx],winner]) #adds the most common structure first
    del frequency[key]
    for i in frequency.keys(): #rest of structures added
        i = Pairs(i)
        indx = tmp_list.index(i)-1
        result.append([tmp_list[indx],i])
        
    return result

def get_names(lines):
    """
    Retrieves the names of the sequences in the output.
    """ 
    next = False 
    names = []
    for line in lines:
        if next:
            if len(line) == 1:
                break
            else:
                tmp = line.split()
                names.append(tmp[1])
        if line.startswith('Sequences loaded ...'):
            next = True
    return names

def minimalComrnaParser(lines):
    """
    Parses the output file in to blocks depending on the S score
    S score is the Maximum stem similarity score.
    """
    block = []
    first = True
    record = False
    for line in lines:
        if line.startswith('===========================  S ='):
            record = True
            if not first:
                yield block
                block = []
            first = False
        if record:    
            block.append(line)
    yield block

def blockParser(block):
    """
    Parses every block of S scores in to blocks of structures
    every S score block has 10 or less structures
    """
    struct = []
    first = True
    record = False
    for line in block:
        if line.startswith('Structure #'):
            record = True
            if not first:
                yield struct
                struct = []
            first = False
        if record:
            struct.append(line)
    yield struct

def structParser(lines):
    """
    Parses a structure block into a block containing the sequens and structures
    lines.
    """
    blc = 0 #blank line counter
    bc = 0 #block counter
    struct = []
    record  = False
    for line in lines:
        if len(line) == 1:
            blc +=1
            record = False
        if blc == 2:
            blc = 0
            bc +=1
            record = True
        if record and bc < 3:
            struct.append(line)

    yield struct

            
def pairsParser(seqBlock,names):
    """
    Takes a structure block and parse that into structures
    """
    for name in names:
        seq = []
        sIndx = [] #start index, where in the line the sequence start
        struct = [] #structure lines
        record = False
        for line in seqBlock:
            if line.startswith(name+' '):
                tmp = line.split()
                #if seq length is shorter then 80 for one seq and longer
                #for another seq the following block will be empty for the
                #shorter sequence. this if statement protects against that
                if len(tmp) == 4: 
                    try:
                        seq.append(tmp[2])#[name,start nr,seq,end nr]
                    except:
                        print 'LINE',line
                        print 'BLOCK', seqBlock
                    sIndx.append(index(line,tmp[2]))            
                    record = True
                else:
                    continue
            else:
                if record:
                    record = False
                    struct.append(line)

###############################################################################
# Construction of the full sequence and structure and then mapping each letter
#in structure to a position

        Fseq = '' #full sequence
        Fstruct = '' #full structure
        for i in range(len(seq)):
            # slice out corresponding structure to sequence
            #so you can get the same index for structure and sequence
            tmpStruct = struct[i][sIndx[i]:(sIndx[i]+len(seq[i]))]
            Fseq = ''.join([Fseq,seq[i]])
            Fstruct = ''.join([Fstruct,tmpStruct])
        #Applies a position to every letter in structure sequence    
        letterPos = zip(range(len(Fseq)),Fstruct)
        
###############################################################################
#Cunstruction of dictionary for where every letter in structure has a list of
#positions corresponding to that of that letter in respect to the sequence

        alphabet = {}
        for pos, letter in letterPos:
            indices = []
            #if the dict contains the letter you want to add to that list
            if alphabet.__contains__(letter): 
                indices = alphabet[letter]
                indices.append(pos)
                alphabet[letter] = indices
            #else you want to create a new list for that letter
            elif not letter==' ':
                indices.append(pos)
                alphabet[letter] = indices
                
###############################################################################
#Each list in alphabet needs to be split in two,
#oL and cL (open and close list), to be able to fold the positions into pairs

        pairs = []
        for value in alphabet.values():
            middle = len(value)/2
            oL = value[:middle]
            cL = value[middle:]
            #pairs are created by making a tuple of the first in oL to
            #the last in cl, second in oL to second last in cL and so on
            pairs.extend(zip(oL,cL.__reversed__()))

        yield Pairs(pairs),Fseq



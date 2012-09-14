#!/usr/bin/env python

"""Provides utility methods for sequence weighting
"""
from __future__ import division

from numpy import array, zeros, dot as matrixmultiply, ones, identity, take,\
    asarray, uint8 as UInt8, add, sqrt
from random import choice
from cogent.util.array import hamming_distance
from cogent.core.profile import Profile, CharMeaningProfile
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import Alignment

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight", "Gavin Huttley", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

PROTEIN_ORDER = ''.join(PROTEIN.Alphabet) + '-'
DNA_ORDER = ''.join(DNA.Alphabet) + '-'
RNA_ORDER = ''.join(RNA.Alphabet) + '-'


class Weights(dict):
    """Dictionary seq_ids: seq weight
    """

    def normalize(self):
        """Normalizes to one. Works in place!"""
        total = sum(self.values(), 0)
        for k, v in self.items():
            self[k] = v/total


def number_of_pseudo_seqs(alignment):
    """Returns the # of unique randomized sequences that can be generated 
    from the alignment.

    alignment: Alignment object

    A random sequence is created by randomly choosing at each position one 
    of the residues observed at that position (column) in the alighment. 
    A single occurrence of that residue type is sufficient to make it an 
    option and the choice is with equal likelihood from any of the observed
    characters. (See Implementation Notes for more details).
    """
    return reduce(lambda x, y: x*y, map(len,alignment.columnProbs()))

def pseudo_seqs_exact(alignment,n=None):
    """Returns all possible pseudo sequences (generated from the alignment)

    alignment: Alignment object
    n: has NO FUNCTION, except to make the API for all functions that generate
    pseudo sequences the same.
    
    This function is used by the VOR method (only when the number of unique
    pseudo sequences is lower than 1000). 

    The original sequences will be generated in this list of pseudo sequences.
    Duplications of sequences in the alignment don't cause duplications
    in the list of pseudo sequences.
    
    Example:
    AA, AA, BB in the alignment
    generates pseudo sequences: AA,BB,AB,BA

    ABC, BCC, BAC in the alignment
    generates pseudo sequences: AAC, ABC, ACC, BAC, BBC, and BCC
    """
    counts = alignment.columnFreqs()
    results = [[]]
    for col in counts:
        new_results = []
        for item in results:
            for option in col:
                new_results.append(item + [option])
        results = new_results
    return [''.join(i) for i in results]

def pseudo_seqs_monte_carlo(alignment,n=1000):
    """Yields sample of possible pseudo sequences (generated from alignment)

    alignment: Alignment object
    n = number of pseudo sequences generated (=sample size)

    This function is used by the VOR method, usually when the number of 
    possible pseudo sequences exceeds 1000.

    To see how the pseudo sequences are generated read the Implementation
    Notes at the top of this file.
    """
    freqs = alignment.columnFreqs()
    for x in range(n):
        seq = []
        for i in freqs:
            seq.append(choice(i.keys()))
        yield ''.join(seq)

def row_to_vote(row):
    """Changes distances to votes.

    There's one vote in total. The sequence with the lowest distance
    gets the vote. If there are multiple sequences equidistant at the 
    minimum distance, the vote is split over all the sequences at the 
    minimum distance.
    
    Example: 
    [4,2,3,7] -> [0,1,0,0]
    [1,3,2,1] -> [.5,0,0,.5]
    [1,1,1,1] -> [.25,.25,.25,.25]
    """
    result = array(row) - min(row) == 0
    return result/sum(result, 0)

def distance_matrix(alignment, distance_method=hamming_distance):
    """Returns distance matrix for seqs in the alignment.

    Order is either the Names in the alignment or the order in which
        is iterated over the rows.

    Distance is the Hamming distance between two sequences
    """
    #change sequences into arrays
    alignment = [array(str(alignment.NamedSeqs[k]), 'c') for k in\
        alignment.Names]

    nr_of_seqs = len(alignment)
    distances = zeros([nr_of_seqs,nr_of_seqs])
    for i, seq_a in enumerate(alignment):
        for j, seq_b in enumerate(alignment):
            if i < j:
                dist = distance_method(seq_a,seq_b)
                distances[i,j] = dist
                distances[j,i] = dist
            else: # i==j (diagonal) or i>j (lower left half)
                continue
    return distances

def eigenvector_for_largest_eigenvalue(matrix):
    """Returns eigenvector corresponding to largest eigenvalue of matrix.
    
    Implements a numerical method for finding an eigenvector by repeated 
    application of the matrix to a starting vector. For a matrix A the 
    process w(k) <-- A*w(k-1) converges to eigenvector w with the largest 
    eigenvalue. Because distance matrix D has all entries >= 0, a theorem 
    due to Perron and Frobenius on nonnegative matrices guarantees good 
    behavior of this method, excepting degenerate cases where the 
    iteration oscillates. For distance matrices a remedy is to add the 
    identity matrix to A, permitting the iteration to converge to the 
    eigenvector. (From Sander and Schneider (1991), and Vingron and 
    Sibbald (1993)) 
    
    Note: Only works on square matrices.
    """
    #always add the identity matrix to avoid oscillating behavior
    matrix = matrix + identity(len(matrix))

    #v is a random vector (chosen as the normalized vector of ones)
    v = ones(len(matrix))/len(matrix)

    #iterate until convergence
    for i in range(1000):
        new_v = matrixmultiply(matrix,v)
        new_v = new_v/sum(new_v, 0) #normalize
        if sum(map(abs,new_v-v), 0) > 1e-9:
            v = new_v #not converged yet
            continue
        else: #converged
            break
    
    return new_v


def distance_to_closest(alignment, distance_method=hamming_distance):
    """Returns vector of distances to closest neighbor for each s in alignment
    
    alignment: Alignment object
    distance_method: function used to calculate the distance between two seqs
    
    Function returns the closest distances according to the Names in the
    alignment

    example:
    Alignment({1:'ABCD',2:'ABCC',3:'CBDD',4:'ACAA'},Names=[3,2,1,4])
    [2,1,1,3]
    """
    names = alignment.Names
    seqs = [array(str(alignment.NamedSeqs[n]), 'c') for n in names]
    
    closest = []
    for i, key in enumerate(names):
        seq = seqs[i]
        dist = None
        for j, other_key in enumerate(names):
            if key == other_key:
                continue
            d = distance_method(seq,seqs[j])
            if dist:
                if d < dist:
                    dist = d
            else:
                dist = d
        closest.append(dist)
    return array(closest)

def SeqToProfile(seq, alphabet=None, char_order=None,\
    split_degenerates=False):
    """Generates a Profile object from a Sequence object.

    seq: Sequence object
    alphabet (optional): Alphabet object (if you want to split
        degenerate symbols, the alphabet object should have a 
        Degenerates property. Default is the Alphabet associated with 
        the Sequence object.
    char_order (optional): The order the characters occur in the Profile.
        Default is the list(alphabet)
    split_degenerates (optional): Whether you want the counts for the 
        degenerate symbols to be divided over the non-degenerate symbols they
        code for.
    
    A Profile is a position x character matrix describing which characters
    occur at each position. In a sequence (as opposed to an alignment) only
    one character occurs at each position. In general a sequence profile
    will only contain ones and zeros. However, you have the possibility of 
    splitting degenerate characters. For example, if a position is R, it 
    means that you have 50/50% chance of A and G. It's also possible to 
    ignore characters, which in a sequence profile will lead to positions
    (rows) containing only zeros.
    
    Example:
    Sequence = ACGU
    Profile(seq, CharOrder=UCAG):
    U   C   A   G
    0   0   1   0   first pos
    0   1   0   0   second pos
    0   0   0   1   third pos
    1   0   0   0   fourth pos

    Sequence= GURY
    Profile(seq,CharOrder=UCAG, split_degenerates=True)
    U   C   A   G
    0   0   0   1   first pos
    1   0   0   0   second pos
    0   0   .5  .5  third pos
    .5  .5  0   0   fourth pos

    Characters can also be ignored
    Sequence = ACN-
    Profile(seq, CharOrder=UCAGN, split_degenerates=True)
    U   C   A   G
    0   0   1   0   first pos
    0   1   0   0   second pos
    .25 .25 .25 .25 third pos
    0   0   0   0   fourth pos <--contains only zeros
    """

    if alphabet is None:
        alphabet = seq.MolType
    if char_order is None:
        char_order = list(alphabet)

    #Determine the meaning of each character based on the alphabet, the
    #character order, and the option to split degenerates
    char_meaning = CharMeaningProfile(alphabet, char_order,\
        split_degenerates)
    #construct profile data
    idxs = array(str(seq).upper(), 'c').view(UInt8)
    result_data = char_meaning.Data[idxs]
    #result_data = take(char_meaning.Data, asarray(str(seq).upper(), UInt8), axis=0)
    
    return Profile(result_data, alphabet, char_order)


def AlnToProfile(aln, alphabet=None, char_order=None, split_degenerates=False,\
    weights=None):
    """Generates a Profile object from an Alignment.

    aln: Alignment object
    alphabet (optional): an Alphabet object (or list of chars, but if you 
        want to split degenerate symbols, the alphabet must have a 
        Degenerates property. Default is the alphabet of the first seq in 
        the alignment.
    char_order (optional): order of the characters in the profile. Default
        is list(alphabet)
    split_degenerates (optional): Whether you want the counts for the 
        degenerate symbols to be divided over the non-degenerate symbols they
        code for.
    weights (optional): dictionary of seq_id: weight. If not entered all seqs
        are weighted equally

    A Profile is a position x character matrix describing which characters
    occur at each position of an alignment. The Profile is always normalized,
    so it gives the probabilities of each character at each position.
    
    Ignoring chars: you can ignore characters in the alignment by not putting
    the char in the CharOrder. If you ignore all characters at a particular
    position, an error will be raised, because the profile can't be normalized.

    Splitting degenerates: you can split degenerate characters over the 
    non-degenerate characters they code for. For example: R = A or G. So,
    an R at a position counts for 0.5 A and 0.5 G.
   
    Example:
    seq1    TCAG    weight: 0.5
    seq2    TAR-    weight: 0.25
    seq3    YAG-    weight: 0.25
    Profile(aln,alphabet=DNA,char_order="TACG",weights=w,
    split_degenerates=True)
    Profile:
       T      A      C      G
    [[ 0.875  0.     0.125  0.   ]
     [ 0.     0.5    0.5    0.   ]
     [ 0.     0.625  0.     0.375]
     [ 0.     0.     0.     1.   ]]
    """

    if alphabet is None:
        alphabet = aln.values()[0].MolType
    if char_order is None:
        char_order = list(alphabet)
    if weights is None:
        weights = dict.fromkeys(aln.keys(),1/len(aln))
    
    char_meaning = CharMeaningProfile(alphabet, char_order,\
        split_degenerates)

    profiles = []
    for k,v in aln.items():
        idxs = array(str(v).upper(), 'c').view(UInt8)
        profiles.append(char_meaning.Data[idxs] * weights[k])
    s = reduce(add,profiles)
    
    result = Profile(s,alphabet, char_order)
    try:
        result.normalizePositions()
    except Exception, e:
        raise ValueError,e
        #"Probably one of the rows in your profile adds up to zero,\n "+\
        #"because you are ignoring all of the characters in the "+\
        #"corresponding\n column in the alignment"
    return result


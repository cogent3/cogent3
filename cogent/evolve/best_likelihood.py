#!/usr/bin/env python
"""
Returns the likelihood resulting from a model in which motif probabilities 
are assumed to be equal to the observed motif frequencies, as described by 
Goldman (1993).  This model is not informative for inferring the evolutionary 
process, but its likelihood indicates the maximum possible likelihood value 
for a site-independent evolutionary process.
"""
from __future__ import division
from numpy import log

from cogent import LoadSeqs

__author__ = "Helen Lindsay, Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Helen Lindsay", "Gavin Huttley", "Daniel McDonald"]
cite = "Goldman, N. (1993).  Statistical tests of models of DNA substitution.  J Mol Evol, 36: 182-98"
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def _transpose(array):
    new_array = []
    num_cols = len(array[0])
    for col in range(num_cols):
        new_row = []
        for row in array:
            new_row += [row[col]]
        new_array.append(new_row)
    return new_array

def _take(array, indices):
    new_array = []
    for index in indices:
        new_array.append(array[index])
    return new_array

def aligned_columns_to_rows(aln, motif_len, exclude_chars = None, allowed_chars='ACGT'):
    """return alignment broken into motifs as a transposed list with
    sequences as columns and aligned columns as rows
    
    Arguments:
        exclude_chars: columns containing these characters will be excluded"""
    if exclude_chars:
        exclude_chars = set(exclude_chars)
        exclude_func = exclude_chars.intersection
    else:
        allowed_chars = set(allowed_chars)
        exclude_func = lambda x: not allowed_chars.issuperset(x)
    
    exclude_indices = set()
    array = []
    for name in aln.Names:
        motifs = list(aln.getGappedSeq(name).getInMotifSize(motif_len))
        array.append(motifs)
        for motif_index, motif in enumerate(motifs):
            if exclude_func(motif):
                    exclude_indices.update([motif_index])
    
    include_indices = set(range(len(array[0]))).difference(exclude_indices)
    include_indices = list(include_indices)
    include_indices.sort()
    array = _transpose(array)
    array = _take(array, include_indices)
    return array

def count_column_freqs(columns_list):
    """return the frequency of columns"""
    col_freq_dict = {}
    for column in columns_list:
        column = ' '.join(column)
        col_freq_dict[column] = col_freq_dict.get(column, 0) + 1
    return col_freq_dict

def get_ML_probs(columns_list, with_patterns=False):
    """returns the column log-likelihoods and frequencies
    
    Argument:
        - with_patterns: the site patterns are returned"""
    n = len(columns_list)
    col_freq_dict = count_column_freqs(columns_list)
    col_lnL_freqs = []
    for column_pattern, freq in col_freq_dict.items():
        # note, the behaviour of / is changed due to the __future__ import
        if with_patterns:
            row = [column_pattern, freq/n, freq]
        else:
            row = [freq/n, freq]
        col_lnL_freqs.append(row)
    return col_lnL_freqs

def get_G93_lnL_from_array(columns_list):
    """return the best log likelihood for a list of aligned columns"""
    col_stats = get_ML_probs(columns_list)
    log_likelihood = 0
    for freq, num in col_stats:
        pattern_lnL = log(freq)*num
        log_likelihood += pattern_lnL
    return log_likelihood

def BestLogLikelihood(aln, alphabet=None, exclude_chars = None,
    allowed_chars='ACGT', motif_length=None, return_length=False):
    """returns the best log-likelihood according to Goldman 1993.
    
    Arguments:
        - alphabet: a sequence alphabet object.
        - motif_length: 1 for nucleotide, 2 for dinucleotide, etc ..
        - exclude_chars: a series of characters used to exclude motifs
        - allowed_chars: only motifs that contain a subset of these are
          allowed
        - return_length: whether to also return the number of alignment columns
    """
    assert alphabet or motif_length, "Must provide either an alphabet or a"\
                                     " motif_length"
    # need to use the alphabet, so we can enforce character compliance
    if alphabet:
        kwargs = dict(moltype=alphabet.MolType)
        motif_length = alphabet.getMotifLen()
    else:
        kwargs = {}
    
    aln = LoadSeqs(data=aln.todict(), **kwargs)
    columns = aligned_columns_to_rows(aln, motif_length, exclude_chars,
                                        allowed_chars)
    num_cols = len(columns)
    log_likelihood = get_G93_lnL_from_array(columns)
    if return_length:
        return log_likelihood, num_cols
    
    return log_likelihood


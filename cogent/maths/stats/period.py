from __future__ import division
from random import shuffle, random, choice
import numpy

try:
    from math import factorial
except ImportError: # python version < 2.6
    from cogent.maths.stats.special import Gamma
    factorial = lambda x: Gamma(x+1)

from cogent.maths.stats.special import igam

__author__ = "Hua Ying, Julien Epps and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Julien Epps", "Hua Ying", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

def chi_square(x, p, df=1):
    """returns the chisquare statistic and it's probability"""
    N = len(x)
    end = N
    sim = numpy.logical_not(numpy.logical_xor(x[0:end-p], x[p:end]))*1
    s = ((numpy.ones((N-p,), float)-sim)**2).sum()
    D = s/(N-p)
    p_val = 1 - igam(df/2.0, D/2)
    return D, p_val

def g_statistic(X, p=None, idx=None):
    """
    return g statistic and p value
    arguments: 
        X - the periodicity profile (e.g. DFT magnitudes, autocorrelation etc)
            X needs to contain only those period values being considered, 
            i.e. only periods in the range [llim, ulim]
    """
    # X should be real
    X = abs(numpy.array(X))
    if p is None: 
        power = X.max(0)
        idx = X.argmax(0)
    else:
        assert idx is not None
        power = X[idx]
    g_obs = power/X.sum()
    M = numpy.floor(1/g_obs)
    pmax = len(X)
    
    result = numpy.zeros((int(M+1),), float)
    pmax_fact = factorial(pmax)
    for index in xrange(1, min(pmax, int(M))+1):
        v = (-1)**(index-1)*pmax_fact/factorial(pmax-index)/factorial(index)
        v *= (1-index*g_obs)**(pmax-1)
        result[index] = v
    
    p_val = result.sum()
    return g_obs, p_val

def _seq_to_symbols(seq, motifs, motif_length, result=None):
    """return symbolic represetation of the sequence
    Arguments: 
        - seq: a sequence
        - motifs: a list of sequence motifs
        - motif_length: length of first motif
    """
    if result is None:
        result = numpy.zeros(len(seq), numpy.uint8)
    else:
        result.fill(0)
    
    if motif_length is None:
        motif_length = len(motifs[0])
    
    for i in xrange(len(seq) - motif_length + 1):
        if seq[i: i + motif_length] in motifs:
            result[i] = 1
    
    return result

try:
    from cogent.maths._period import seq_to_symbols
    #raise ImportError
except ImportError:
    seq_to_symbols = _seq_to_symbols

class SeqToSymbols(object):
    """class for converting all occurrences of motifs in passed sequence
    to 1/0 otherwise"""
    def __init__(self, motifs, length=None, motif_length=None):
        super(SeqToSymbols, self).__init__()
        if type(motifs) == str:
            motifs = [motifs]
        self.motifs = motifs
        self.length = length
        self.motif_length = motif_length or len(motifs[0])
        self.working = None
        if length is not None:
            self.setResultArray(length)
    
    def setResultArray(self, length):
        """sets a result array for length"""
        self.working = numpy.zeros(length, numpy.uint8)
        self.length = length
    
    def __call__(self, seq, result=None):
        if result is None and self.working is None:
            self.setResultArray(len(seq))
        elif self.working is not None:
            if len(seq) != self.working.shape[0]:
                self.setResultArray(len(seq))
        
        result = self.working
        result.fill(0)
        if type(seq) != str:
            seq = ''.join(seq)
        
        return seq_to_symbols(seq, self.motifs, self.motif_length, result)
    

def circular_indices(vector, start, length, num):
    """docstring for circular_indices"""
    if start > length:
        start = start-length
        
    if start+num < length:
        return vector[start: start+num]
    # get all till end, then from beginning
    return vector[start:] + vector[:start+num-length]

def sampled_places(block_size, length):
    """returns randomly sampled positions with block_size to make a new vector
    with length
    """
    # Main condition is to identify when a draw would run off end, we want to
    # draw from beginning
    num_seg, remainder = divmod(length, block_size)
    vector = range(length)
    result = []
    for seg_num in xrange(num_seg):
        i = choice(vector)
        result += circular_indices(vector, i, length, block_size)
    
    if remainder:
        result += circular_indices(vector, i+block_size, length, remainder)
    
    assert len(result) == length, len(result)
    return result

def blockwise_bootstrap(signal, calc, block_size, num_reps, seq_to_symbols=None, num_stats=None):
    """returns observed statistic and the probability from the bootstrap
    test of observing more `power' by chance than that estimated from the
    observed signal
    
    Arguments:
        - signal: a series, can be a sequence object
        - calc: function to calculate the period power, e.g. ipdft, hybrid,
          auto_corr or any other statistic.
        - block_size: size of contiguous values for resampling
        - num_reps: number of randomly generated permutations
        - seq_to_symbols: function to convert a sequence to 1/0. If not
          provided, the raw data is used.
        - num_stats: the number of statistics being evaluated for each
          interation. Default to 1.
    """
    signal_length = len(signal)
    
    if seq_to_symbols is not None:
        dtype='c'
    else:
        dtype=None # let numpy guess
    
    signal = numpy.array(list(signal), dtype=dtype)
    
    if seq_to_symbols is not None:
        symbolic = seq_to_symbols(signal)
        data = symbolic
    else:
        data = signal
    
    obs_stat = calc(data)
    if seq_to_symbols is not None:
        if sum(symbolic) == 0:
            p = [numpy.array([1.0, 1.0, 1.0]), 1.0][num_stats == 1]
            
            return obs_stat, p
    
    if num_stats is None:
        try:
            num_stats = calc.getNumStats()
        except AttributeError:
            num_stats = 1
    
    if num_stats == 1:
        count = 0
    else:
        count = numpy.zeros(num_stats)
    
    for rep in range(num_reps):
        # get sample positions
        sampled_indices = sampled_places(block_size, signal_length)
        new_signal = signal.take(sampled_indices)
        if seq_to_symbols is not None:
            symbolic = seq_to_symbols(new_signal)
            data = symbolic
        else:
            data = new_signal
        sim_stat = calc(data)
        # count if > than observed
        if num_stats > 1:
            count[sim_stat >= obs_stat] += 1
        elif sim_stat >= obs_stat:
            count += 1
        
    return obs_stat, count / num_reps


# def percrb4():
#     """Return SNR and CRB for periodicity estimates from symbolic signals"""
#     # TODO: complete the function according to Julien's percrb4.m
#     pass
# 

from __future__ import division
from numpy import log, zeros, float64, int32, array, sqrt, dot, diag, where
from numpy.linalg import det, norm, inv

from cogent import DNA, RNA, LoadTable
from cogent.util.progress_display import display_wrap

__author__ = "Gavin Huttley and Yicheng Zhu"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha" # pending addition of protein distance metrics

def _same_moltype(ref, query):
    """if ref and query have the same states"""
    return set(ref) == set(query)

def get_pyrimidine_indices(moltype):
    """returns pyrimidine indices for the moltype"""
    states = list(moltype)
    if _same_moltype(RNA, moltype):
        return map(states.index, 'CU')
    elif _same_moltype(DNA, moltype):
        return map(states.index, 'CT')
    else:
        raise RuntimeError('Non-nucleic acid MolType')

def get_purine_indices(moltype):
    """returns purine indices for the moltype"""
    states = list(moltype)
    if not _same_moltype(RNA, moltype) and not _same_moltype(DNA, moltype):
        raise RuntimeError('Non-nucleic acid MolType')
    
    return map(states.index, 'AG')

def get_matrix_diff_coords(indices):
    """returns coordinates for off diagonal elements"""
    return [(i,j) for i in indices for j in indices if i != j]

def get_moltype_index_array(moltype, invalid=-9):
    """returns the index array for a molecular type"""
    canonical_chars = list(moltype)
    # maximum ordinal for an allowed character, this defines the length of
    # the required numpy array
    max_ord = max(map(ord, moltype.All.keys()))
    char_to_index = zeros(max_ord+1, int32)
    # all non canonical_chars are ``invalid''
    char_to_index.fill(invalid)
    
    for i in range(len(canonical_chars)):
        c = canonical_chars[i]
        o = ord(c)
        char_to_index[o] = i
    
    return char_to_index

def seq_to_indices(seq, char_to_index):
    """returns an array with sequence characters replaced by their index"""
    ords = map(ord, seq)
    indices = char_to_index.take(ords)
    return indices

def _fill_diversity_matrix(matrix, seq1, seq2):
    """fills the diversity matrix for valid positions.
    
    Assumes the provided sequences have been converted to indices with
    invalid characters being negative numbers (use get_moltype_index_array
    plus seq_to_indices)."""
    paired = array([seq1, seq2]).T
    paired = paired[paired.min(axis=1) >= 0]
    for i in range(len(paired)):
        matrix[paired[i][0], paired[i][1]] += 1
    

def _jc69_from_matrix(matrix):
    """computes JC69 stats from a diversity matrix"""
    invalid = None, None, None, None
    total = matrix.sum()
    diffs = total - sum(matrix[i,i] for i in range(matrix.shape[0]))
    if total == 0:
        return invalid
    
    p = diffs / total
    if p >= 0.75: # cannot take log
        return invalid
    
    factor = (1 - (4 / 3) * p)
    
    dist = -3.0 * log(factor) / 4
    var = p * (1 - p) / (factor * factor * total)
    return total, p, dist, var

def _tn93_from_matrix(matrix, freqs, pur_indices, pyr_indices, pur_coords, pyr_coords, tv_coords):
    invalid = None, None, None, None
    
    total = matrix.sum()
    freqs = matrix.sum(axis=0) + matrix.sum(axis=1)
    freqs /= (2*total)
    
    if total == 0:
        return invalid
    
    # 
    p = matrix.take(pur_coords + pyr_coords + tv_coords).sum() / total
    
    freq_purs = freqs.take(pur_indices).sum()
    prod_purs = freqs.take(pur_indices).prod()
    freq_pyrs = freqs.take(pyr_indices).sum()
    prod_pyrs = freqs.take(pyr_indices).prod()
    
    # purine transition diffs
    pur_ts_diffs = matrix.take(pur_coords).sum()
    pur_ts_diffs /= total
    # pyr transition  diffs
    pyr_ts_diffs = matrix.take(pyr_coords).sum()
    pyr_ts_diffs /= total
    # transversions
    tv_diffs = matrix.take(tv_coords).sum() / total
    
    coeff1 = 2 * prod_purs / freq_purs
    coeff2 = 2 * prod_pyrs / freq_pyrs
    coeff3 = 2 * (freq_purs * freq_pyrs - \
            (prod_purs * freq_pyrs / freq_purs) -\
            (prod_pyrs * freq_purs / freq_pyrs))
    
    
    term1 = 1 - pur_ts_diffs / coeff1 - tv_diffs / (2*freq_purs)
    term2 = 1 - pyr_ts_diffs / coeff2 - tv_diffs / (2*freq_pyrs)
    term3 = 1 - tv_diffs / (2 * freq_purs * freq_pyrs)
    
    if term1 <= 0 or term2 <= 0 or term3 <= 0: # log will fail
        return invalid
    
    dist = -coeff1 * log(term1) - coeff2 * log(term2) - coeff3 * log(term3)
    v1 = 1 / term1
    v2 = 1 / term2
    v3 = 1 / term3
    v4 = (coeff1 * v1 / (2 * freq_purs)) + \
         (coeff2 * v2 / (2 * freq_pyrs)) + \
         (coeff3 * v3 / (2 * freq_purs * freq_pyrs))
    var = v1**2 * pur_ts_diffs + v2**2 * pyr_ts_diffs + v4**2 * tv_diffs - \
          (v1 * pur_ts_diffs + v2 * pyr_ts_diffs + v4 * tv_diffs)**2
    var /= total
    
    return total, p, dist, var


def _logdet(matrix, use_tk_adjustment=True):
    """returns the LogDet from a diversity matrix
    Arguments:
        - use_tk_adjustment: when True, unequal state frequencies are allowed
    """
    
    invalid = None, None, None, None
    total = matrix.sum()
    diffs = total - sum(matrix[i,i] for i in range(matrix.shape[0]))
    if total == 0:
        return invalid
    
    p = diffs / total
    
    if diffs == 0: # seqs identical
        return total, p, 0.0, None
    
    # we replace missing diagonal states with a frequency of 0.5,
    # then normalise
    frequency = matrix.copy()
    unobserved = where(frequency.diagonal() == 0)[0]
    for index in unobserved:
        frequency[index, index] = 0.5
    
    frequency /= frequency.sum()
    # the inverse matrix of frequency, every element is squared
    M_matrix = inv(frequency)**2
    freqs_1 = frequency.sum(axis = 0)
    freqs_2 = frequency.sum(axis = 1)
    
    if use_tk_adjustment:
        mean_state_freqs = (freqs_1 + freqs_2) / 2
        coeff = (norm(mean_state_freqs)**2 - 1) / (matrix.shape[0] - 1)
    else:
        coeff = -1 / matrix.shape[0]
        
    
    FM_1 = diag(freqs_1)
    FM_2 = diag(freqs_2)
    
    try:
        d_xy = coeff * log(det(frequency) / sqrt(det(FM_1 * FM_2)))
    except FloatingPointError:
        return invalid
    
    if det(frequency) <= 0: #if the result is nan
        return invalid
    
    var_term = dot(M_matrix, frequency).transpose()[0].sum()
    var_denom = 16 * total
    if use_tk_adjustment:
        var = (var_term - (1 / sqrt(freqs_1 * freqs_2)).sum()) / var_denom
    else:
        # variance formula for TK adjustment is false
        var = (var_term - 1) / var_denom
    
    var = d_xy - 2 * var
    
    return total, p, d_xy, var


try:
    from _pairwise_distance import \
        _fill_diversity_matrix as fill_diversity_matrix
    # raise ImportError # for testing
except ImportError:
    fill_diversity_matrix = _fill_diversity_matrix

def _number_formatter(template):
    """flexible number formatter"""
    def call(val):
        try:
            result = template % val
        except TypeError:
            result = val
        return result
    return call

class _PairwiseDistance(object):
    """base class for computing pairwise distances"""
    
    def __init__(self, moltype, invalid=-9, alignment=None):
        super(_PairwiseDistance, self).__init__()
        self.moltype = moltype
        self.char_to_indices = get_moltype_index_array(moltype)
        self._dim = len(list(moltype))
        self._dists = None
        
        self.Names = None
        self.IndexedSeqs = None
        
        if alignment is not None:
            self._convert_seqs_to_indices(alignment)
        
        self._func_args = []
    
    def _convert_seqs_to_indices(self, alignment):
        assert type(alignment.MolType) == type(self.moltype), \
            'Alignment does not have correct MolType'
        
        self._dists = {}
        self.Names = alignment.Names[:]
        indexed_seqs = []
        for name in self.Names:
            seq = alignment.getGappedSeq(name)
            indexed = seq_to_indices(str(seq), self.char_to_indices)
            indexed_seqs.append(indexed)
        
        self.IndexedSeqs = array(indexed_seqs)
    
    @staticmethod
    def func():
        pass # over ride in subclasses
    
    @display_wrap
    def run(self, alignment=None, ui=None):
        """computes the pairwise distances"""
        if alignment is not None:
            self._convert_seqs_to_indices(alignment)
        
        matrix = zeros((self._dim, self._dim), float64)
        done = 0.0
        to_do = (len(self.Names) * len(self.Names) - 1) / 2
        for i in range(len(self.Names)-1):
            name_1 = self.Names[i]
            s1 = self.IndexedSeqs[i]
            for j in range(i+1, len(self.Names)):
                name_2 = self.Names[j]
                ui.display('%s vs %s' % (name_1, name_2), done / to_do )
                done += 1
                matrix.fill(0)
                s2 = self.IndexedSeqs[j]
                fill_diversity_matrix(matrix, s1, s2)
                total, p, dist, var = self.func(matrix, *self._func_args)
                self._dists[(name_1, name_2)] = (total, p, dist, var)
                self._dists[(name_2, name_1)] = (total, p, dist, var)
        
    
    def getPairwiseDistances(self):
        """returns a 2D dictionary of pairwise distances."""
        if self._dists is None:
            return None
        dists = {}
        for name_1 in self.Names:
            for name_2 in self.Names:
                if name_1 == name_2:
                    continue
                val = self._dists[(name_1, name_2)][2]
                dists[(name_1, name_2)] = val
                dists[(name_2, name_1)] = val
        
        return dists
    
    def _get_stats(self, stat, transform=None, **kwargs):
        """returns a table for the indicated statistics"""
        if self._dists is None:
            return None
        rows = []
        for row_name in self.Names:
            row = [row_name]
            for col_name in self.Names:
                if row_name == col_name:
                    row.append('')
                    continue
                val = self._dists[(row_name, col_name)][stat]
                if transform is not None:
                    val = transform(val)
                row.append(val)
            rows.append(row)
        header = [r'Seq1 \ Seq2'] + self.Names
        table = LoadTable(header=header, rows=rows, row_ids = True,
                missing_data='*', **kwargs)
        return table
    
    @property
    def Dists(self):
        kwargs = dict(title='Pairwise Distances', digits=4)
        return self._get_stats(2, **kwargs)
    
    @property
    def StdErr(self):
        stderr = lambda x: sqrt(x)
        kwargs = dict(title='Standard Error of Pairwise Distances', digits=4)
        return self._get_stats(3, transform=stderr, **kwargs)
    
    @property
    def Variances(self):
        kwargs = dict(title='Variances of Pairwise Distances', digits=4)
        table = self._get_stats(3, **kwargs)
        var_formatter = _number_formatter("%.2e")
        if table is not None:
            for name in self.Names:
                table.setColumnFormat(name, var_formatter)
        
        return table
    
    @property
    def Proportions(self):
        kwargs = dict(title='Proportion variable sites', digits=4)
        return self._get_stats(1, **kwargs)
    
    @property
    def Lengths(self):
        kwargs = dict(title='Pairwise Aligned Lengths', digits=0)
        return self._get_stats(0, **kwargs)
    

class _NucleicSeqPair(_PairwiseDistance):
    """docstring for _NucleicSeqPair"""
    def __init__(self, *args, **kwargs):
        super(_NucleicSeqPair, self).__init__(*args, **kwargs)
        if not _same_moltype(DNA, self.moltype) and \
            not _same_moltype(RNA, self.moltype):
            raise RuntimeError('Invalid MolType for this metric')
    

class JC69Pair(_NucleicSeqPair):
    """calculator for pairwise alignments"""
    def __init__(self, *args, **kwargs):
        """states: the valid sequence states"""
        super(JC69Pair, self).__init__(*args, **kwargs)
        self.func = _jc69_from_matrix
    

class TN93Pair(_NucleicSeqPair):
    """calculator for pairwise alignments"""
    def __init__(self, *args, **kwargs):
        """states: the valid sequence states"""
        super(TN93Pair, self).__init__(*args, **kwargs)
        self._freqs = zeros(self._dim, float64)
        self.pur_indices = get_purine_indices(self.moltype)
        self.pyr_indices = get_pyrimidine_indices(self.moltype)
        
        # matrix coordinates
        self.pyr_coords = get_matrix_diff_coords(self.pyr_indices)
        self.pur_coords = get_matrix_diff_coords(self.pur_indices)
        
        self.tv_coords = get_matrix_diff_coords(range(self._dim))
        for coord in self.pur_coords + self.pyr_coords:
            self.tv_coords.remove(coord)
        
        # flattened
        self.pyr_coords = [i * 4 + j for i, j in self.pyr_coords]
        self.pur_coords = [i * 4 + j for i, j in self.pur_coords]
        self.tv_coords = [i * 4 + j for i, j in self.tv_coords]
        
        self.func = _tn93_from_matrix
        self._func_args = [self._freqs, self.pur_indices,
            self.pyr_indices, self.pur_coords,
            self.pyr_coords, self.tv_coords]
        
    

class LogDetPair(_PairwiseDistance):
    """computes logdet distance between sequence pairs"""
    def __init__(self, use_tk_adjustment=True, *args, **kwargs):
        """Arguments:
            - use_tk_adjustment: use the correction of Tamura and Kumar 2002
        """
        super(LogDetPair, self).__init__(*args, **kwargs)
        self.func = _logdet
        self._func_args = [use_tk_adjustment]
    
    def run(self, use_tk_adjustment=None, *args, **kwargs):
        if use_tk_adjustment is not None:
            self._func_args = [use_tk_adjustment]
        
        super(LogDetPair, self).run(*args, **kwargs)

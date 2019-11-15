cimport numpy as np

__version__ = "('2019', '11', '15', 'a')"

# fills in a diversity matrix from sequences of integers
def _fill_diversity_matrix(np.ndarray[np.float64_t, ndim=2] matrix, np.ndarray[np.int32_t, ndim=1] seq1, np.ndarray[np.int32_t, ndim=1] seq2):
    """fills the diversity matrix for valid positions.
    
    Assumes the provided sequences have been converted to indices with
    invalid characters being negative numbers (use get_moltype_index_array
    plus seq_to_indices)."""
    cdef int i
    
    for i in range(len(seq1)):
        if seq1[i] < 0 or seq2[i] < 0:
            continue
        
        matrix[seq1[i], seq2[i]] += 1.0
    


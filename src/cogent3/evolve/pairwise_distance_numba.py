from numba import njit


__author__ = "Gavin Huttley, Yicheng Zhu and Ben Kaehler"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Yicheng Zhu", "Ben Kaehler", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


# fills in a diversity matrix from sequences of integers
@njit(cache=True)
def fill_diversity_matrix(matrix, seq1, seq2):
    """fills the diversity matrix for valid positions.

    Assumes the provided sequences have been converted to indices with
    invalid characters being negative numbers (use get_moltype_index_array
    plus seq_to_indices)."""

    for i in range(len(seq1)):
        if seq1[i] < 0 or seq2[i] < 0:
            continue
        matrix[seq1[i], seq2[i]] += 1.0

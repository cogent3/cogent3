from numba import njit

# turn off code coverage as njit-ted code not accessible to coverage


# fills in a diversity matrix from sequences of integers
@njit(cache=True)
def fill_diversity_matrix(matrix, seq1, seq2):  # pragma: no cover
    """fills the diversity matrix for valid positions.

    Assumes the provided sequences have been converted to indices with
    invalid characters being negative numbers (use get_moltype_index_array
    plus seq_to_indices)."""

    for i in range(len(seq1)):
        if seq1[i] < 0 or seq2[i] < 0:
            continue
        matrix[seq1[i], seq2[i]] += 1.0

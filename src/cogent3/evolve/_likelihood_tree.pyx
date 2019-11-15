#cython: boundscheck=False
#cython: wraparound=False

include "../../include/numerical_pyrex.pyx"
version_info = (2, 2)
__version__ = "('2019', '11', '15', 'a')"

cdef extern from "math.h":
    double log (double x)


def sum_input_likelihoods(child_indexes, Double2D result, likelihoods):
    cdef int M, S, U, C, motif, parent_col, child_col, child
    cdef Double2D plhs
    cdef Long1D index

    # S is parent seq length, U is unique columns in child seq
    # M is size of alphabet, C is number of children.
    C = len(child_indexes)
    M = S = 0
    checkArray2D(result, &S, &M)
    
    for child in range(C):
        index = child_indexes[child]
        plhs = likelihoods[child]
        U = 0
        checkArray1D(index, &S)
        checkArray2D(plhs, &U, &M)
        if child == 0:
            for parent_col in range(S):
                child_col = index[parent_col]
                for motif in range(M):
                    result[parent_col, motif] = plhs[child_col, motif]
        else:
            for parent_col in range(S):
                child_col = index[parent_col]
                for motif in range(M):
                    result[parent_col, motif] *= plhs[child_col, motif]
    return result
    
def get_total_log_likelihood(Double1D counts, Double2D input_likelihoods, Double1D mprobs):
    cdef int S, M, col, motif
    cdef double posn, total
    
    # M is size of alphabet, S is seq length
    S = M = 0
    checkArray1D(mprobs, &M)
    checkArray1D(counts, &S)
    checkArray2D(input_likelihoods, &S, &M)
    
    total = 0.0
    for col in range(S):
        posn = 0.0
        for motif in range(M):
            posn += input_likelihoods[col, motif] * mprobs[motif]
        total += log(posn)*counts[col]
    return total

def get_log_sum_across_sites(Double1D counts, Double1D input_likelihoods):
    cdef int S, col
    cdef double total
    
    S = 0
    checkArray1D(counts, &S)
    checkArray1D(input_likelihoods, &S)
    
    total = 0.0
    for col in range(S):
        total += log(input_likelihoods[col])*counts[col]
    return total

def log_dot_reduce(Long1D index, object patch_probs, Double2D switch_probs, Double2D plhs):
    cdef int i, j, col, site, N, U, S, most_probable_state
    cdef int exponent
    cdef double result, BASE
    cdef Double1D state, prev, tmp
    cdef object patch_probs1, patch_probs2
    BASE = 2.0 ** 1000
    patch_probs1 = patch_probs.copy()
    patch_probs2 = patch_probs.copy()
    state = patch_probs1
    prev = patch_probs2
    
    # S is seq length, U is unique columns in child seq
    # N is number of patch types
    N = U = S = 0
    checkArray1D(state, &N)
    checkArray1D(prev, &N)
    checkArray2D(switch_probs, &N, &N)
    checkArray2D(plhs, &U, &N)
    checkArray1D(index, &S)
    
    exponent = 0
    for site in range(S):
        col = index[site]
        if col >= U:
            raise ValueError((col, U))
        tmp = prev
        prev = state
        state = tmp
        most_probable_state = 0
        for i in range(N):
            state[i] = 0
            for j in range(N):
                state[i] += prev[j] * switch_probs[j, i]
            state[i] *= plhs[col, i]
            if state[i] > state[most_probable_state]:
                most_probable_state = i
        while state[most_probable_state] < 1.0:
            for i from 0 <= i < N:
                state[i] *= BASE
            exponent += -1
    result = 0.0
    for i in range(N):
        result += state[i]

    return log(result) + exponent * log(BASE)
        

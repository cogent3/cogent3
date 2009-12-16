include "../../include/numerical_pyrex.pyx"
version_info = (2, 1)
__version__ = "('1', '4')"

cdef extern from "math.h":
    double log (double x)

def sumInputLikelihoods(child_indexes, result, likelihoods):
    # M is dim of alphabet, S is non-redundandt parent seq length, 
    # U is length
    cdef int M, S, U, m, i, u
    cdef int c
    
    cdef double *values_data
    cdef long *index_data
    cdef double *target_data
    
    M = S = 0
    target_data = checkArrayDouble2D(result, &S, &M)
    first = 1
    c = 0
    for index in child_indexes:
        U = 0
        index_data = checkArrayLong1D(index, &S)
        values_data = checkArrayDouble2D(likelihoods[c], &U, &M)
        #if index_data[S-1] >= U:
        #    raise RangeError
        if c == 0:
            for i from 0 <= i < S:
                u = index_data[i]
                for m from 0 <= m < M:
                    target_data[M*i+m] = values_data[M*u+m]
        else:
            for i from 0 <= i < S: # col of parent data
                u = index_data[i] # col of childs data
                for m from 0 <= m < M:
                    target_data[M*i+m] *= values_data[M*u+m]
        c += 1
    return result
    
def getTotalLogLikelihood(counts, input_likelihoods, mprobs):
    cdef int S, M, i, m
    cdef double posn, total
    cdef double *likelihoods_data, *mprobs_data, *weights_data
    
    S = M = 0
    mprobs_data = checkArrayDouble1D(mprobs, &M)
    weights_data = checkArrayDouble1D(counts, &S)
    likelihoods_data = checkArrayDouble2D(input_likelihoods, &S, &M)
    total = 0.0
    for i from 0 <= i < S:
        posn = 0.0
        for m from 0 <= m < M:
            posn += likelihoods_data[i*M+m] * mprobs_data[m]
        total += log(posn)*weights_data[i]
    return total

def getLogSumAcrossSites(counts, input_likelihoods):
    cdef int S, i
    cdef double total
    cdef double *likelihoods_data,  *weights_data
    S = 0
    weights_data = checkArrayDouble1D(counts, &S)
    likelihoods_data = checkArrayDouble1D(input_likelihoods, &S)
    total = 0.0
    for i from 0 <= i < S:
        total += log(likelihoods_data[i])*weights_data[i]
    return total

def logDotReduce(index, patch_probs, switch_probs, plhs):
    cdef int site, i, j, k, n, uniq, exponent, length, most_probable_state
    cdef double result, BASE
    cdef double *sp, *pl, *state, *prev, *tmp
    cdef long *index_data
    cdef object patch_probs1, patch_probs2
    BASE = 2.0 ** 1000
    patch_probs1 = patch_probs.copy()
    patch_probs2 = patch_probs.copy()
    n = uniq = length = 0
    state = checkArrayDouble1D(patch_probs1, &n)
    prev = checkArrayDouble1D(patch_probs2, &n)
    sp = checkArrayDouble2D(switch_probs, &n, &n)
    pl = checkArrayDouble2D(plhs, &uniq, &n)
    index_data = checkArrayLong1D(index, &length)
    exponent = 0
    for site from 0 <= site < length:
        k = index_data[site]
        if k >= uniq:
            raise ValueError((k, uniq))
        tmp = prev
        prev = state
        state = tmp
        most_probable_state = 0
        for i from 0 <= i < n:
            state[i] = 0
            for j from 0 <= j < n:
                state[i] += prev[j] * sp[j*n+i]
            state[i] *= pl[k*n+i]
            if state[i] > state[most_probable_state]:
                most_probable_state = i
        while state[most_probable_state] < 1.0:
            for i from 0 <= i < n:
                state[i] *= BASE
            exponent += -1
    result = 0.0
    for i from 0 <= i < n:
        result += state[i]

    return log(result) + exponent * log(BASE)
        

from numpy import pi, exp, sqrt, cos
cimport numpy as np

version_info = (3, 2)
__version__ = "('2019', '11', '15', 'a')"

# TODO intro the version idea of peter's see email from him on Wednesday, 26 May 2010

def goertzel_inner(np.ndarray[np.float64_t, ndim=1] x, int N, int period):
    """returns the power from series x for period"""
    cdef int n
    cdef np.float64_t coeff, s, s_prev, s_prev2, power
    
    coeff = 2.0 * cos(2 * pi / period)
    s_prev = 0.0
    s_prev2 = 0.0
    for n in range(N):
        s = x[n] + coeff * s_prev - s_prev2
        s_prev2 = s_prev
        s_prev = s
    
    power = sqrt(s_prev2**2 + s_prev**2 - coeff * s_prev2 * s_prev)
    return power

def ipdft_inner(np.ndarray[np.float64_t, ndim=1] x,
                       np.ndarray[np.complex128_t, ndim=1] X,
                       np.ndarray[np.complex128_t, ndim=1] W,
                       int ulim, int N):
    """use this when repeated calls for window of same length are to be
    made"""
    cdef int n, p
    cdef np.complex128_t w
    
    for p in range(ulim):
        w = 1.0
        for n in range(N):
            if n != 0:
                w = w * W[p]
            X[p] = X[p] + x[n]*w
    return X

def autocorr_inner(np.ndarray[np.float64_t, ndim=1] x, np.ndarray[np.float64_t, ndim=1] xc, int N):
    cdef int m, n
    
    for m in range(-N+1, N):
        for n in range(N):
            if 0 <= n-m < N:
                xc[m+N-1] += (x[n]*x[n-m])

def seq_to_symbols(char* seq, list motifs, int motif_length,
    np.ndarray[np.uint8_t, ndim=1] result):
    cdef int i, j, N, num_motifs
    cdef bytes got
    
    N = len(seq)
    num_motifs = len(motifs)
    motif_length = len(motifs[0])
    for i in range(N - motif_length + 1):
        got = seq[i: i+motif_length]
        for j in range(num_motifs):
            if got == motifs[j]:
                result[i] = 1
    
    return result
    
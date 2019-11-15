#cython: boundscheck=False
#cython: wraparound=False

include "../../include/numerical_pyrex.pyx"

# The weird indentation in this file is to match the POG version.

cdef extern from "Python.h":
    PyErr_Occurred()
    int PyErr_CheckSignals()
    double Py_HUGE_VAL

cdef extern from "math.h":
    double log (double x)

version_info = (3, 2)
__version__ = "('2019', '11', '15', 'a')"

cdef double SCALE_STEP, MIN_FLOAT_VALUE
SCALE_STEP = 2.0**50
MIN_FLOAT_VALUE = 1.0 / SCALE_STEP

cdef int MAX_XCOUNT
MAX_XCOUNT = 256

cdef long MIN_SCALE, MAX_SCALE
MIN_SCALE = -10000
MAX_SCALE = +10000  # or 0 if all numbers should be probabilities

def fmpt(mantissa, exponent, msg=''):
    return "%s * SCALE_STEP ** %s %s" % (mantissa, exponent, msg)

ctypedef unsigned char [:,:,::1] UChar3D

def calc_rows(Long1D plan, Long1D seq1_index, Long1D seq2_index, 
        int i_low, int i_high, int j_low, int j_high, preds, 
        Long2D state_directions, Double2D T, 
        Double2D xgap_scores, Double2D ygap_scores, Double3D match_scores, 
        rows, UChar3D track, track_enc, int viterbi, int use_logs=0, int local=False, 
        int use_scaling=True):
    
    """The faster, sequence only (no POG) version.  Forward or Viterbi 
    algorithm, with doubles or with slower but practically unoverflowable 
    (double, long) GMP-like numbers.  Viterbi is also available in the ever 
    popular addition-of-logs version.  All this with any possible pair HMM 
    transition matrix.
        
    Limitations
       - HMM states must be in a sensible order: M and X, then Y, then END.
    """
    
    # These are array lengths/indicies and so could be Py_ssize_t
    cdef int prev_i, prev_j, state, prev_state, min_prev_state, N
    cdef int row_length, row_length1, row_count, row_count1, tmp_rows
    cdef int a_count, b_count, a, b, a_low, a_high, b_low, b_high
    cdef int dest_states, dest_state, d4, j, i
    cdef int last_i, last_j, last_state
    cdef int bin, x, y, bin_count, max_x, max_y
    cdef int current_row_index, source_row_index
    cdef int source_i

    cdef int dx, dy
    cdef int tcode_x, tcode_y, tcode_s
    cdef double d_score, mantissa, partial_sum, sub_partial_sum, max_mantissa, overall_max_mantissa
    cdef long exponent, max_exponent, overall_max_exponent
    cdef Double3D mantissas
    cdef Long3D exponents
    cdef long pointer_a, pointer_b, pointer_state
    
    assert not (use_logs and not viterbi)
    assert not (use_logs and use_scaling)
    assert not (local and not viterbi)
    
    N = 0
    checkArray2D(T, &N, &N)
    row_length = 0
    row_count = 0
    checkArray1D(plan, &row_count)
    
    dest_states = 0
    d4 = 4
    # Array of (state, bin, dx, dy) tuples describing the HMM states.
    checkArray2D(state_directions, &dest_states, &d4)
        
    checkArray1D(seq1_index, &row_count)
    checkArray1D(seq2_index, &row_length)
    
    max_x = max_y = bin_count = 0
    checkArray3D(match_scores, &bin_count, &max_x, &max_y)
    checkArray2D(xgap_scores, &bin_count, &max_x)
    checkArray2D(ygap_scores, &bin_count, &max_y)
    
    for i from 0 <= i < row_count:
        assert 0 <= seq1_index[i] < max_x
    for j from 0 <= j < row_length:
        assert 0 <= seq2_index[j] < max_y
    
    assert j_low >= 0 and j_high > j_low and j_high <= row_length
    
    (mantissas, exponents) = rows
    tmp_rows = 0
    checkArray3D(mantissas, &tmp_rows, &row_length, &N)
    if use_scaling:
        checkArray3D(exponents, &tmp_rows, &row_length, &N)
    
    cdef double impossible                
    if use_logs:
        impossible = log(0.0) # -inf
    else:
        impossible = 0.0

    if viterbi and track is not None and track_enc is not None:
        checkArray3D(track, &row_count, &row_length, &N)
        (tcode_x, tcode_y, tcode_s) = track_enc
    else:
        track = None
        tcode_x = tcode_y = tcode_s = 0
    
    # For local
    overall_max_exponent = MIN_SCALE
    overall_max_mantissa = impossible
    last_i = last_j = last_state = -1
    
    for i from i_low <= i < i_high:
        x = seq1_index[i]
        
        if PyErr_CheckSignals():
            raise PyErr_Occurred()
        
        current_row_index = plan[i]
        
        #for prev_state from 1 <= prev_state < N:
        #    current_row_data[0, prev_state] = impossible
        
        for j from j_low <= j < j_high:
            
            for dest_state from 0 <= dest_state < dest_states:
                state = state_directions[dest_state, 0]
                bin = state_directions[dest_state, 1]
                dx = state_directions[dest_state, 2]
                dy = state_directions[dest_state, 3]
                
                max_mantissa = impossible
                max_exponent = MIN_SCALE
                partial_sum = 0.0
                pointer_state = N  # ie ERROR
                
                source_i = i - dx
                if source_i < 0:
                    continue
                source_row_index = plan[source_i]
                
                prev_j = j - dy
                if prev_j < 0:
                    continue
                
                min_prev_state = 1
                a = dx
                b = dy
                
                if (local and dx and dy) or (prev_j == 0 and source_i == 0):
                    partial_sum = max_mantissa = T[0, state]
                    max_exponent = 0
                    pointer_state = 0
                    pointer_a = a
                    pointer_b = b
                    
                if use_scaling:
                            sub_partial_sum = 0.0
                            for prev_state from min_prev_state <= prev_state < N:
                                exponent = exponents[source_row_index, prev_j, prev_state]
                                if exponent == MIN_SCALE:
                                    continue
                                
                                mantissa = mantissas[source_row_index, prev_j, prev_state] 
                                mantissa = mantissa * T[prev_state, state]
                                
                                if mantissa < MIN_FLOAT_VALUE:
                                    if mantissa == 0.0:
                                        continue
                                    if mantissa < 0.0:
                                        if T[prev_state, state] < 0.0:
                                            raise ArithmeticError(fmpt(mantissa, exponent, 
                                                    "transition is a negative probability"))
                                        raise ArithmeticError(fmpt(mantissa, exponent, 
                                                "product is a negative probability"))
                                    while mantissa < MIN_FLOAT_VALUE:
                                        mantissa *= SCALE_STEP
                                        exponent += -1
                                        if exponent <= MIN_SCALE:
                                          raise ArithmeticError(fmpt(mantissa, exponent,
                                                "underflows"))
                                                
                                elif mantissa > 1.0:
                                    mantissa *= MIN_FLOAT_VALUE
                                    exponent += 1
                                    if exponent > MAX_SCALE:
                                        raise ArithmeticError(fmpt(mantissa, exponent, 
                                            "is unexpectedly large"))
                                                
                                if exponent > max_exponent:
                                    if exponent == max_exponent + 1:
                                        sub_partial_sum = partial_sum
                                    else:
                                        sub_partial_sum = 0.0
                                    partial_sum = 0.0
                                    max_mantissa = 0.0
                                    max_exponent = exponent
                                
                                if exponent == max_exponent:
                                    partial_sum += mantissa
                                    if viterbi and mantissa > max_mantissa:
                                        max_mantissa = mantissa
                                        pointer_state = prev_state
                                        pointer_a = a
                                        pointer_b = b
                                
                                elif exponent == max_exponent - 1:
                                    sub_partial_sum += mantissa
                                    
                            partial_sum += sub_partial_sum * MIN_FLOAT_VALUE
                else:
                            for prev_state from min_prev_state <= prev_state < N:
                                mantissa = mantissas[source_row_index, prev_j, prev_state] 
                                if use_logs:
                                    mantissa = mantissa + T[prev_state, state]
                                else:
                                    mantissa = mantissa * T[prev_state, state]
                                    partial_sum += mantissa
                                if viterbi and mantissa > max_mantissa:
                                    max_mantissa = mantissa
                                    pointer_state = prev_state
                                    pointer_a = a
                                    pointer_b = b
                        
                if viterbi:
                    mantissa = max_mantissa
                    if track is not None:
                        track[i, j, state] = (        
                            (pointer_a << tcode_x) |  
                            (pointer_b << tcode_y) |
                            (pointer_state << tcode_s))
                else:
                    mantissa = partial_sum
    
                if dy:
                    y = seq2_index[j]
                    if dx:
                        d_score = match_scores[bin, x, y]
                    else:
                        d_score = ygap_scores[bin, y]
                elif dx:
                    d_score = xgap_scores[bin, x]
                elif use_logs:
                    d_score = 0.0
                else:
                    d_score = 1.0
                
                if use_logs:
                    mantissa += d_score
                else:
                    mantissa *= d_score
                
                mantissas[current_row_index, j, state] = mantissa
                if use_scaling:
                    exponents[current_row_index, j, state] = max_exponent
                    
                if local and dx and dy:
                    if (use_scaling and max_exponent > overall_max_exponent) or (
                            (not use_scaling or max_exponent == overall_max_exponent) and (
                            mantissa > overall_max_mantissa)):
                        overall_max_exponent = max_exponent
                        overall_max_mantissa = mantissa
                        last_i = i
                        last_j = j
                        last_state = state
    if not local:
        last_i = i_high - 1
        last_j = j_high - 1
        last_state = state
    else:
        mantissa = overall_max_mantissa
        max_exponent = overall_max_exponent
        
    if use_scaling:
        score = log(mantissa) + log(SCALE_STEP) * max_exponent
    elif use_logs:
        score = mantissa
    else:
        score = log(mantissa)
    return ((last_i, last_j), last_state, score)

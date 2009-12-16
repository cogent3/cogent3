include "numerical_pyrex.pyx"

cdef extern from "Python.h":
    PyErr_Occurred()
    int PyErr_CheckSignals()
    double Py_HUGE_VAL

cdef extern from "math.h":
    double log (double x)

version_info = (3, 0)
__version__ = "('1', '4')"

cdef double SCALE_STEP, MIN_FLOAT_VALUE
SCALE_STEP = 2.0**50
MIN_FLOAT_VALUE = 1.0 / SCALE_STEP

cdef int MAX_XCOUNT
MAX_XCOUNT = 256

cdef long MIN_SCALE, MAX_SCALE
MIN_SCALE = -10000
MAX_SCALE = +10000  # or 0 if all numbers should be probabilities

#cdef unsigned long * checkArrayULong3D(ArrayType a, int *x, int *y, int *z) except NULL:
#    return <unsigned long *> checkArray3D(a, c'i', sizeof(long), x, y, z)

#cdef unsigned int * checkArrayUInt3D(ArrayType a, int *x, int *y, int *z) except NULL:
#    return <unsigned int *> checkArray3D(a, c'i', sizeof(int), x, y, z)

cdef unsigned char * checkArrayUChar3D(ArrayType a, int *x, int *y, int *z) except NULL:
    return <unsigned char *> checkArray3D(a, c'i', sizeof(char), x, y, z)

def fmpt(mantissa, exponent, msg=''):
    return "%s * SCALE_STEP ** %s %s" % (mantissa, exponent, msg)

def calc_rows(ArrayType plan, ArrayType seq1_index, ArrayType seq2_index, 
        int i_low, int i_high, int j_low, int j_high, preds, 
        ArrayType state_directions, ArrayType T, 
        ArrayType xgap_scores, ArrayType ygap_scores, ArrayType match_scores, 
        rows, track, track_enc, int viterbi, int use_logs=0, int local=False, 
        int use_scaling=True):
    
    """The ultimate in 2D Pyrex dynamic programming - Forward or Viterbi 
    algorithm, with doubles or with slower but practically unoverflowable 
    (double, long) GMP-like numbers.  Viterbi is also available in the ever 
    popular addition-of-logs version.  All this with any possible pair HMM 
    transition matrix.
    
    One time to use something faster than this is when the inputs are sequences
    rather than alignments.  This code expects alignments (which can be single
    sequences) represented as POGs (ie: DAGs).
    
    Limitations
       - HMM states must be in a sensible order: M and X, then Y, then END.
    """
        
    cdef int dx, dy, prev_i, prev_j, state, prev_state, N, row_length
    cdef int a_count, b_count, a, b, a_low, a_high, b_low, b_high
    cdef int dest_states, dest_state, d4, j, i, source_i
    cdef int last_i, last_j, last_state, overall_max_exponent
    cdef int tcode_x, tcode_y, tcode_s
    cdef unsigned char *track_data
    cdef double overall_max_mantissa
    cdef double d_score, mantissa, partial_sum, sub_partial_sum, max_mantissa
    cdef long exponent, index, max_exponent, *source_row_ex_data
    cdef double *T_data, *source_row_data
    cdef double *match_score_data
    cdef long *dest_states_data, *plan_data
    cdef ArrayType i_sources, j_sources
    cdef double *current_row_data, *source_row_data_cache[256] # MAX_XCOUNT
    cdef long *current_row_ex_data, *source_row_ex_data_cache[256] # MAX_XCOUNT
    cdef long pointer_a, pointer_b, pointer_state
    cdef long *x_index, *y_index
    cdef int x, y, max_x, max_y
    
    cdef long *i_sources_data, *i_sources_offsets_data
    cdef long *j_sources_data, *j_sources_offsets_data
    cdef long i_sources_start, i_sources_end
    cdef long j_sources_start, j_sources_end
    cdef int i_link_count, j_link_count
    cdef int row_count, row_length1, plan_index, row_count1
    
    (mantissas, exponents) = rows

    assert not (use_logs and not viterbi)
    assert not (use_logs and use_scaling)
    assert not (local and not viterbi)
    
    N = 0
    T_data = checkArrayDouble2D(T, &N, &N)
    row_length = 0
    row_count = 0
    plan_data = checkArrayLong1D(plan, &row_count)
    
    dest_states = 0
    d4 = 4
    # Array of (state, bin, dx, dy) tuples describing the HMM states.
    dest_states_data = checkArrayLong2D(state_directions, &dest_states, &d4)
    
    cdef int bin_count, bin
    cdef double *xgap_score_data, *ygap_score_data
    
    x_index = checkArrayLong1D(seq1_index, &row_count)
    y_index = checkArrayLong1D(seq2_index, &row_length)
    
    max_x = max_y = bin_count = 0
    match_score_data = checkArrayDouble3D(
            match_scores, &bin_count, &max_x, &max_y)
    xgap_score_data = checkArrayDouble2D(xgap_scores, &bin_count, &max_x)
    ygap_score_data = checkArrayDouble2D(ygap_scores, &bin_count, &max_y)
    
    for i from 0 <= i < row_count:
        assert 0 <= x_index[i] < max_x
    for j from 0 <= j < row_length:
        assert 0 <= y_index[j] < max_y
    
    assert j_low >= 0 and j_high > j_low and j_high <= row_length
    
    (pog1, pog2) = preds
    (j_sources, j_sources_offsets) = pog2.asCombinedArray()
    j_link_count = 0
    j_sources_data = checkArrayLong1D(j_sources, &j_link_count)
    row_length1 = row_length + 1
    j_sources_offsets_data = checkArrayLong1D(j_sources_offsets, &row_length1)
    
    (i_sources, i_sources_offsets) = pog1.asCombinedArray()
    i_link_count = 0
    i_sources_data = checkArrayLong1D(i_sources, &i_link_count)
    row_count1 = row_count + 1 # ???
    i_sources_offsets_data = checkArrayLong1D(i_sources_offsets, &row_count1)

    cdef double impossible                
    if use_logs:
        impossible = log(0.0) # -inf
    else:
        impossible = 0.0

    if viterbi and track is not None and track_enc is not None:
        track_data = checkArrayUChar3D(track, &row_count, &row_length, &N)
        (tcode_x, tcode_y, tcode_s) = track_enc
    else:
        track_data = NULL
        tcode_x = tcode_y = tcode_s = 0
    
    # For local
    overall_max_exponent = MIN_SCALE
    overall_max_mantissa = impossible
    last_i = last_j = last_state = -1
    
    for i from i_low <= i < i_high:
        x = x_index[i]
        
        if PyErr_CheckSignals():
            raise PyErr_Occurred()

        plan_index = plan_data[i]
        current_row_data = checkArrayDouble2D(mantissas[plan_index], &row_length, &N)
        if use_scaling:
            current_row_ex_data = checkArrayLong2D(exponents[plan_index], &row_length, &N)
        else:
            current_row_ex_data = NULL
            
        i_sources_start = i_sources_offsets[i]
        i_sources_end = i_sources_offsets[i+1]
    
        source_row_data_cache[0] = current_row_data
        source_row_ex_data_cache[0] = current_row_ex_data
        
        a_count = i_sources_end-i_sources_start
        for a from 0 <= a < a_count:
            prev_i = i_sources_data[a+i_sources_start]
            plan_index = plan_data[prev_i]
            source_row_data_cache[a+1] = checkArrayDouble2D(
                    mantissas[plan_index], &row_length, &N)
            if use_scaling:
                source_row_ex_data_cache[a+1] = checkArrayLong2D(
                        exponents[plan_index], &row_length, &N)
            else:
                source_row_ex_data_cache[a+1] = NULL
        
        if i == 0:
            if use_logs:
                current_row_data[0] = 0.0       
            else:
                current_row_data[0] = 1.0       
            if use_scaling:
                current_row_ex_data[0] = 0             
        else:
            current_row_data[0*N+0] = impossible       
            if use_scaling:
                current_row_ex_data[0*N+0] = MIN_SCALE             
     
        j_sources_end = j_sources_offsets_data[j_low]
        for j from j_low <= j < j_high:                
            j_sources_start = j_sources_end
            j_sources_end = j_sources_offsets_data[j+1]

            for dest_state from 0 <= dest_state < dest_states:
                state = dest_states_data[dest_state*4+0]
                bin = dest_states_data[dest_state*4+1]
                dx = dest_states_data[dest_state*4+2]
                dy = dest_states_data[dest_state*4+3]
                
                max_mantissa = impossible
                max_exponent = MIN_SCALE
                partial_sum = 0.0
                pointer_state = N  # ie ERROR
                
                if dx:
                    a_low = 1
                    a_high = a_count + 1
                else:
                    a_low = 0
                    a_high = 1
                
                if dy:
                    b_low = 1
                    b_high = j_sources_end - j_sources_start + 1
                else:
                    b_low = 0
                    b_high = 1
    
                if use_scaling:
                    sub_partial_sum = 0.0
                    # keep these next 8 lines same as below, plus source_row_ex_data
                    for a from a_low <= a < a_high:
                        source_row_data = source_row_data_cache[a]
                        source_row_ex_data = source_row_ex_data_cache[a]
                        for b from b_low <= b < b_high:
                            if dy:
                                prev_j = j_sources_data[b-1+j_sources_start]
                            else:
                                prev_j = j
                            for prev_state from (prev_j>0) <= prev_state < N:
                                index = prev_j*N + prev_state
                                
                                exponent = source_row_ex_data[index]
                                if exponent == MIN_SCALE:
                                    continue
                                
                                mantissa = (source_row_data[index] 
                                     * T_data[prev_state*N+state])
                                
                                if mantissa < MIN_FLOAT_VALUE:
                                    if mantissa == 0.0:
                                        continue
                                    if mantissa < 0.0:
                                        if T_data[prev_state*N+state] < 0.0:
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
                    # keep these next 7 lines same as above/below, 
                    # less source_row_ex_data
                    for a from a_low <= a < a_high:
                        source_row_data = source_row_data_cache[a]
                        for b from b_low <= b < b_high:
                            if dy:
                                prev_j = j_sources_data[b-1+j_sources_start]
                            else:
                                prev_j = j
                            for prev_state from (prev_j>0) <= prev_state < N:
                                index = prev_j*N + prev_state
                                if use_logs:
                                    mantissa = (source_row_data[index] 
                                         + T_data[prev_state*N+state])
                                else:
                                    mantissa = (source_row_data[index] 
                                         * T_data[prev_state*N+state])
                                    partial_sum += mantissa
                                if viterbi and mantissa > max_mantissa:
                                    max_mantissa = mantissa
                                    pointer_state = prev_state
                                    pointer_a = a
                                    pointer_b = b
                        
                if viterbi:
                    mantissa = max_mantissa
                    if track_data:
                        track_data[(i*row_length+j)*N+state] = (        
                            (pointer_a << tcode_x) |  
                            (pointer_b << tcode_y) |
                            (pointer_state << tcode_s))
                else:
                    mantissa = partial_sum
    
                if dy:
                    y = y_index[j]
                    if dx:
                        d_score = match_score_data[((bin*max_x+x)*max_y)+y]
                    else:
                        d_score = ygap_score_data[bin*max_y+y]
                elif dx:
                    d_score = xgap_score_data[bin*max_x+x]
                elif use_logs:
                    d_score = 0.0
                else:
                    d_score = 1.0
                
                if use_logs:
                    mantissa += d_score
                else:
                    mantissa *= d_score
                
                current_row_data[j*N+state] = mantissa
                if use_scaling:
                    current_row_ex_data[j*N+state] = max_exponent
                    
                if local and dx and dy:
                    if (use_scaling and max_exponent > overall_max_exponent) or (
                            (not use_scaling or max_exponent == overall_max_exponent) and (
                            mantissa >= overall_max_mantissa)):
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

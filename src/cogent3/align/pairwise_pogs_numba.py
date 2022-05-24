import numpy as np

from numba import boolean, float64, int64, njit, optional, uint8
from numba.core.types.containers import Tuple


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


@njit(
    Tuple(types=(Tuple(types=(int64, int64)), int64, float64,))(
        int64[::1],
        int64[::1],
        int64[::1],
        int64,
        int64,
        int64,
        int64,
        optional(int64[::1]),
        optional(int64[::1]),
        optional(int64[::1]),
        optional(int64[::1]),
        int64[:, ::1],
        float64[:, ::1],
        float64[:, ::1],
        float64[:, ::1],
        float64[:, :, ::1],
        optional(float64[:, :, ::1]),
        float64,
        optional(int64[:, :, ::1]),
        optional(uint8[:, :, ::1]),
        optional(int64[::1]),
        boolean,
        boolean,
        boolean,
        boolean,
    ),
    cache=True,
)
def calc_rows(
    plan,
    x_index,
    y_index,
    i_low,
    i_high,
    j_low,
    j_high,
    i_sources,
    i_sources_offsets,
    j_sources,
    j_sources_offsets,
    state_directions,
    T,
    xgap_scores,
    ygap_scores,
    match_scores,
    mantissas,
    mantissa,
    exponents,
    track,
    track_enc,
    viterbi,
    local=False,
    use_scaling=False,
    use_logs=False,
):
    assert not (use_logs and not viterbi)
    assert not (use_logs and use_scaling)
    assert not (local and not viterbi)

    MIN_SCALE = -10000
    MAX_SCALE = +10000
    SCALE_STEP = 2.0 ** 50
    MIN_FLOAT_VALUE = 1.0 / SCALE_STEP
    source_row_index_cache = np.zeros(256)

    N = max(T.shape[0], T.shape[1])

    dest_states = max(0, state_directions.shape[0])

    row_count = x_index.shape[0]
    row_length = y_index.shape[0]

    max_x = match_scores.shape[1]
    max_y = match_scores.shape[2]

    max_x = max(xgap_scores.shape[1], max_x)
    max_y = max(ygap_scores.shape[1], max_y)

    for i in range(row_count):
        assert 0 <= x_index[i] <= max_x

    for j in range(row_length):
        assert 0 <= y_index[j] <= max_y

    assert j_low >= 0 and j_high > j_low and j_high <= row_length

    row_length = max(mantissas.shape[1], row_length)
    N = max(mantissas.shape[2], N)

    if use_scaling:
        row_length = max(exponents.shape[1], row_length)
        N = max(exponents.shape[2], N)

    if use_logs:
        impossible = -np.inf
    else:
        impossible = 0.0

    if viterbi and track is not None and track_enc is not None:
        N = max(track.shape[2], N)
        (tcode_x, tcode_y, tcode_s) = track_enc
    else:
        track = None
        tcode_x = tcode_y = tcode_s = 0

    overall_max_exponent = MIN_SCALE
    overall_max_mantissa = impossible
    last_i = last_j = last_state = -1

    max_exponent = MIN_SCALE

    for i in range(i_low, i_high):
        x = x_index[i]

        i_sources_start = i_sources_offsets[i]
        i_sources_end = i_sources_offsets[i + 1]

        current_row_index = plan[i]
        source_row_index_cache[0] = current_row_index

        a_count = i_sources_end - i_sources_start
        for a in range(a_count):
            prev_i = i_sources[a + i_sources_start]
            source_row_index_cache[a + 1] = plan[prev_i]

        if i == 0:
            if use_logs:
                mantissas[current_row_index, 0, 0] = 0.0
            else:
                mantissas[current_row_index, 0, 0] = 1.0
            if use_scaling:
                exponents[current_row_index, 0, 0] = 0
        else:
            mantissas[current_row_index, 0, 0] = impossible
            if use_scaling:
                exponents[current_row_index, 0, 0] = MIN_SCALE

        j_sources_end = j_sources_offsets[j_low]
        for j in range(j_low, j_high):
            j_sources_start = j_sources_end
            j_sources_end = j_sources_offsets[j + 1]

            for dest_state in range(dest_states):
                state = state_directions[dest_state, 0]
                bin = state_directions[dest_state, 1]
                dx = state_directions[dest_state, 2]
                dy = state_directions[dest_state, 3]

                max_mantissa = impossible
                max_exponent = MIN_SCALE
                partial_sum = 0.0
                pointer_state = N

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

                pointer_a = 0
                pointer_b = 0

                if use_scaling:
                    sub_partial_sum = 0.0

                    for a in range(a_low, a_high):
                        source_row_index = int(source_row_index_cache[a])
                        for b in range(b_low, b_high):
                            if dy:
                                prev_j = j_sources[b - 1 + j_sources_start]
                            else:
                                prev_j = j
                            min_prev_state = prev_j > 0

                            for prev_state in range(min_prev_state, N):
                                exponent = exponents[
                                    source_row_index, prev_j, prev_state
                                ]
                                if exponent == MIN_SCALE:
                                    continue

                                transition = T[prev_state, state]

                                mantissa = mantissas[
                                    source_row_index, prev_j, prev_state
                                ]
                                mantissa *= transition

                                if mantissa < MIN_FLOAT_VALUE:
                                    if mantissa == 0.0:
                                        continue
                                    assert mantissa >= 0.0 and transition >= 0.0

                                    while mantissa < MIN_FLOAT_VALUE:
                                        mantissa *= SCALE_STEP
                                        exponent += -1
                                        assert exponent > MIN_SCALE

                                elif mantissa > 1.0:
                                    mantissa *= MIN_FLOAT_VALUE
                                    exponent += 1
                                    assert exponent <= MAX_SCALE

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
                    for a in range(a_low, a_high):
                        source_row_index = int(source_row_index_cache[a])
                        for b in range(b_low, b_high):
                            if dy:
                                prev_j = j_sources[b - 1 + j_sources_start]
                            else:
                                prev_j = j
                            min_prev_state = prev_j > 0

                            for prev_state in range(min_prev_state, N):
                                mantissa = mantissas[
                                    source_row_index, prev_j, prev_state
                                ]
                                transition = T[prev_state, state]
                                if use_logs:
                                    mantissa += transition
                                else:
                                    mantissa *= transition
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
                            (pointer_a << tcode_x)
                            | (pointer_b << tcode_y)
                            | (pointer_state << tcode_s)
                        )
                else:
                    mantissa = partial_sum

                if dy:
                    y = y_index[j]
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
                        (not use_scaling or max_exponent == overall_max_exponent)
                        and (mantissa > overall_max_mantissa)
                    ):
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
        score = np.log(mantissa) + np.log(SCALE_STEP) * max_exponent
    elif use_logs:
        score = mantissa
    else:
        score = np.log(mantissa)
    return ((last_i, last_j), last_state, score)

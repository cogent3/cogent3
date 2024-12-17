import enum
import itertools
import typing

import numba
import numpy
from numpy import (
    nan,
)

from cogent3.core import new_alphabet
from cogent3.core.moltype import IUPAC_gap, IUPAC_missing
from cogent3.util import dict_array
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI

DEFAULT_EXCLUDES = f"{IUPAC_gap}{IUPAC_missing}"
DEFAULT_NULL_VALUE = nan


class MI_METHODS(enum.Enum):
    mi = "mi"
    nmi = "nmi"
    rmi = "rmi"


# Comments on design
# the revised mutual information calculations are based on a sequence alignment
# represented as a numpy uint8 array.

# The resampled mutual information calculation uses a cache of all entropy terms
# from the independent and joint position, produced by _calc_entropy_components().
# For each combination of alternate possible states, there are two values in the
# entropy terms that need to be modified. These calculations are done by
# _calc_updated_entropy() and _calc_temp_entropy(). This caching, plus the
# numba.jit compilation, makes rmi competitive performance-wise with the other
# methods.


@numba.jit
def _count_states(
    state_vector: numpy.ndarray,
    num_states: int,
    counts: numpy.ndarray | None = None,
) -> numpy.ndarray:  # pragma: no cover
    """computes counts from a single vector of states"""
    if counts is None:
        counts = numpy.empty(num_states, dtype=numpy.int64)

    counts.fill(0)
    for state in state_vector:
        if state < num_states:
            counts[state] += 1

    return counts


@numba.jit
def _vector_entropy(counts: numpy.ndarray) -> float:  # pragma: no cover
    """computes entropy for a single vector of integers"""
    total = counts.sum()
    if total <= 1:
        return 0.0 if total else numpy.nan

    entropy = 0.0
    for count in counts:
        if count > 0:
            prob = count / total
            entropy += prob * -numpy.log2(prob)

    return entropy


@numba.jit
def _count_joint_states(
    joint_states: numpy.ndarray,
    num_states: int,
    counts: numpy.ndarray | None = None,
) -> numpy.ndarray:  # pragma: no cover
    if counts is None:
        counts = numpy.empty((num_states, num_states), dtype=numpy.int64)

    counts.fill(0)
    for joint_state in joint_states:
        i, j = joint_state
        if i >= num_states or j >= num_states:
            continue

        counts[i, j] += 1
    return counts


@numba.jit
def _calc_joint_entropy(counts: numpy.ndarray) -> float:  # pragma: no cover
    entropy = 0.0
    total_counts = counts.sum()
    for count in counts.flatten():
        if count > 0:
            prob = count / total_counts
            entropy += prob * -numpy.log2(prob)
    return entropy


@numba.jit
def _calc_column_entropies(
    columns: numpy.ndarray,
    num_states: int,
) -> numpy.ndarray:  # pragma: no cover
    """
    Calculate the entropy for each column in the input array.

    Parameters:
    array (numpy.ndarray): Input array of unsigned 8-bit integers.

    Returns:
    numpy.ndarray: Array of entropy values for each column.
    """
    n_cols = columns.shape[1]
    entropies = numpy.zeros(n_cols, dtype=numpy.float64)
    counts = numpy.empty(num_states, dtype=numpy.int64)
    for col in range(n_cols):
        counts = _count_states(columns[:, col], num_states, counts)
        entropies[col] = _vector_entropy(counts)

    return entropies


@numba.jit
def _make_weights(
    counts: numpy.ndarray,
    weights: numpy.ndarray | None = None,
) -> numpy.ndarray:  # pragma: no cover
    """Return the weights for replacement states for each possible character.
    We compute the weight as the normalized frequency of the replacement state
    divided by 2*n."""
    zeroes = counts == 0
    total = counts.sum()
    char_prob = counts.astype(numpy.float64) / total
    if weights is None:
        weights = numpy.empty((counts.shape[0], counts.shape[0]), dtype=numpy.float64)

    weights.fill(0.0)
    denom = 2 * total
    for i in range(counts.shape[0]):
        if zeroes[i]:
            continue
        diag = char_prob[i]
        weights[i, :] = char_prob / (1 - diag) / denom
        weights[i, i] = 0.0

    return weights


@numba.jit
def _calc_entropy_components(
    counts: numpy.ndarray,
    total: int,
) -> tuple[float, numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    """Return the entropy and arrays of entropy components and non-zero status"""
    non_zero = counts != 0
    assert counts.ndim == 1, "designed for 1D arrays"
    freqs = counts.astype(numpy.float64) / total
    log2 = numpy.zeros(counts.shape, dtype=numpy.float64)
    et = numpy.zeros(counts.shape, dtype=numpy.float64)
    log2[non_zero] = numpy.log2(freqs[non_zero])
    et[non_zero] = (
        -log2[non_zero] * freqs[non_zero]
    )  # the terms in the entropy equation
    h = numpy.sum(et)  # entropy of the original data
    return h, et, non_zero


@numba.jit
def _calc_temp_entropy(
    entropy: float,
    counts: numpy.ndarray,
    entropy_terms: numpy.ndarray,
    total: int,
    index: int,
) -> float:  # pragma: no cover
    # compute the intermediate column 1 entropy term
    new_count = counts[index] - 1
    orig_term = entropy_terms[index]
    if new_count > 0:
        freq = new_count / total
        log2 = -numpy.log2(freq)
        new_term = freq * log2
    else:
        new_term = 0.0

    return entropy - orig_term + new_term


@numba.jit
def _calc_updated_entropy(
    temp_entropy: float,
    counts: numpy.ndarray,
    entropy_terms: numpy.ndarray,
    total: int,
    index: int,
) -> float:  # pragma: no cover
    new_count = counts[index] + 1
    orig_term = entropy_terms[index]
    freq = new_count / total
    log2 = -numpy.log2(freq)
    new_term = freq * log2
    return temp_entropy - orig_term + new_term


@numba.jit
def _calc_pair_scale(
    counts_12: numpy.ndarray,
    counts_1: numpy.ndarray,
    counts_2: numpy.ndarray,
    weights_1: numpy.ndarray,
    weights_2: numpy.ndarray,
    states_12: numpy.ndarray,
    states_1: numpy.ndarray,
    states_2: numpy.ndarray,
    coeffs: numpy.ndarray,
) -> tuple[float, numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    """Return entropies and weights for comparable alignment.
    A comparable alignment is one in which, for each paired state ij, all
    alternate observable paired symbols are created. For instance, let the
    symbols {A,C} be observed at position i and {A,C} at position j. If we
    observe the paired types {AC, AA}. A comparable alignment would involve
    replacing an AC pair with a CC pair."""

    # break down the joint entropy into individual terms so we can easily adjust
    # for the different combinations
    total = counts_1.sum()
    counts_12 = counts_12.flatten()
    je_orig, jet, j_nz = _calc_entropy_components(counts_12, total)

    # individual entropy components for column 1
    orig_e_1, et_1, nz_1 = _calc_entropy_components(counts_1, total)

    # individual entropy components for column 2
    orig_e_2, et_2, nz_2 = _calc_entropy_components(counts_2, total)
    orig_mi = orig_e_1 + orig_e_2 - je_orig

    num_scales = j_nz.sum() * nz_1.sum() + j_nz.sum() * nz_2.sum()
    scales = numpy.zeros((num_scales, 2), dtype=numpy.float64)
    pairs = numpy.zeros((num_scales, 2), dtype=numpy.uint64)
    new_coord = numpy.zeros(2, dtype=numpy.int64)
    if orig_e_1 == 0.0 or orig_e_2 == 0.0 or je_orig == 0.0:
        return 0.0, pairs, scales

    n = 0
    for pair in states_12.flatten()[j_nz]:
        i, j = new_alphabet.index_to_coord(pair, coeffs=coeffs)
        if counts_12[pair] == 0:
            continue

        # compute the intermediate column 1 entropy term
        tmp_e_1 = _calc_temp_entropy(orig_e_1, counts_1, et_1, total, i)

        # compute the intermediate joint entropy term
        tmp_je = _calc_temp_entropy(je_orig, counts_12, jet, total, pair)

        for k in states_1:
            if k == i or counts_1[k] == 0:
                continue
            # compute the new entropy for column 1
            new_e_1 = _calc_updated_entropy(tmp_e_1, counts_1, et_1, total, k)

            # compute the new joint-entropy
            new_coord[:] = k, j
            new_index = new_alphabet.coord_to_index(new_coord, coeffs=coeffs)
            n_je = _calc_updated_entropy(tmp_je, counts_12, jet, total, new_index)

            # the weight
            w = weights_1[i, k]
            new_mi = new_e_1 + orig_e_2 - n_je
            scales[n][0] = new_mi
            scales[n][1] = w
            pairs[n][0] = i
            pairs[n][1] = j
            n += 1

        # compute the intermediate column 2 entropy term
        tmp_e_2 = _calc_temp_entropy(orig_e_2, counts_2, et_2, total, j)
        for k in states_2:
            if k == j or counts_2[k] == 0:
                continue

            # compute the new entropy for column 1
            new_e_2 = _calc_updated_entropy(tmp_e_2, counts_2, et_2, total, k)

            # compute the new joint-entropy
            new_coord[:] = i, k
            new_index = new_alphabet.coord_to_index(new_coord, coeffs=coeffs)
            n_je = _calc_updated_entropy(tmp_je, counts_12, jet, total, new_index)

            # the weight
            w = weights_2[j, k]
            new_mi = orig_e_1 + new_e_2 - n_je
            scales[n][0] = new_mi
            scales[n][1] = w
            pairs[n][0] = i
            pairs[n][1] = j
            n += 1

    return orig_mi, pairs, scales


@numba.jit
def _count_col_joint(
    pos_12: numpy.ndarray,
    counts_12: numpy.ndarray,
    counts_1: numpy.ndarray,
    counts_2: numpy.ndarray,
    num_states: int,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    counts_12.fill(0)
    counts_1.fill(0)
    counts_2.fill(0)
    for pos_1, pos_2 in pos_12:
        if pos_1 >= num_states or pos_2 >= num_states:
            continue
        counts_12[pos_1, pos_2] += 1
        counts_1[pos_1] += 1
        counts_2[pos_2] += 1

    return counts_12, counts_1, counts_2


@numba.jit
def _rmi_calc(
    positions: numpy.ndarray,
    alignment: numpy.ndarray,
    num_states: numpy.uint8,
) -> tuple[numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    coeffs = new_alphabet.coord_conversion_coeffs(num_states, 2, dtype=numpy.int64)
    weights_1 = numpy.empty((num_states, num_states), dtype=float)
    weights_2 = numpy.empty((num_states, num_states), dtype=float)

    n_seqs = alignment.shape[0]
    stats = numpy.empty(len(positions), dtype=numpy.float64)
    stats.fill(numpy.nan)

    counts_1 = numpy.empty(num_states, dtype=numpy.int64)
    counts_2 = numpy.empty(num_states, dtype=numpy.int64)
    counts_12 = numpy.zeros((num_states, num_states), dtype=numpy.int64)
    joint_states = numpy.empty((n_seqs, 2), dtype=numpy.uint8)

    # the state indices
    states_12 = numpy.arange(counts_12.size).reshape(counts_12.shape)
    states_1 = numpy.arange(counts_1.size)
    states_2 = numpy.arange(counts_2.size)

    for pair in range(len(positions)):
        i, j = positions[pair]
        joint_states[:, 0] = alignment[:, i]
        joint_states[:, 1] = alignment[:, j]

        counts_12, counts_1, counts_2 = _count_col_joint(
            joint_states,
            counts_12,
            counts_1,
            counts_2,
            num_states,
        )

        weights_1 = _make_weights(counts_1, weights_1)
        weights_2 = _make_weights(counts_2, weights_2)
        entropy, pairs, scales = _calc_pair_scale(
            counts_12,
            counts_1,
            counts_2,
            weights_1,
            weights_2,
            states_12,
            states_1,
            states_2,
            coeffs,
        )
        if entropy == 0.0:
            stats[pair] = 0.0
        else:
            stat = 0.0
            for i in range(scales.shape[0]):
                e, w = scales[i]
                # we round the revised entropy to avoid floating point errors
                # this is in effect a more stringent condition
                if entropy > numpy.round(e, 10):
                    continue

                p1, p2 = pairs[i]
                stat += w * counts_12[p1, p2]

            stats[pair] = 1 - stat

    return positions, stats


@numba.jit
def _calc_all_entropies(
    joint_states: numpy.ndarray,
    joint_counts: numpy.ndarray,
    counts_1: numpy.ndarray,
    counts_2: numpy.ndarray,
    num_states: int,
) -> tuple[float, float, float]:  # pragma: no cover
    joint_counts, counts_1, counts_2 = _count_col_joint(
        joint_states,
        joint_counts,
        counts_1,
        counts_2,
        num_states,
    )
    entropy_1 = _vector_entropy(counts_1)
    entropy_2 = _vector_entropy(counts_2)
    if entropy_1 == 0.0 or entropy_2 == 0.0:
        return entropy_1, entropy_2, 0.0

    joint_entropy = _calc_joint_entropy(joint_counts)
    return entropy_1, entropy_2, joint_entropy


@numba.jit
def _general_mi_calc(
    positions: numpy.ndarray,
    canonical_pos: numpy.ndarray,
    alignment: numpy.ndarray,
    entropies: numpy.ndarray,
    num_states: numpy.uint8,
    metric_id: int = MI_METHODS.mi,
) -> tuple[numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    n_seqs = alignment.shape[0]
    stats = numpy.empty(len(positions), dtype=numpy.float64)
    stats.fill(numpy.nan)

    counts_1 = numpy.empty(num_states, dtype=numpy.int64)
    counts_2 = numpy.empty(num_states, dtype=numpy.int64)
    joint_counts = numpy.zeros((num_states, num_states), dtype=numpy.int64)
    joint_states = numpy.empty((n_seqs, 2), dtype=numpy.uint8)
    for pair in range(len(positions)):
        i, j = positions[pair]
        joint_states[:, 0] = alignment[:, i]
        joint_states[:, 1] = alignment[:, j]
        if canonical_pos[i] and canonical_pos[j]:
            entropy_i = entropies[i]
            entropy_j = entropies[j]
            joint_counts = _count_joint_states(joint_states, num_states, joint_counts)
            joint_entropy = _calc_joint_entropy(
                joint_counts,
            )
        else:
            # we need to compute all entropies
            # cases where either position has a non-canonical
            # state are omitted
            entropy_i, entropy_j, joint_entropy = _calc_all_entropies(
                joint_states,
                joint_counts,
                counts_1,
                counts_2,
                num_states,
            )

        # MI
        stat = entropy_i + entropy_j - joint_entropy
        if metric_id == 2 and joint_entropy != 0.0:
            # normalised MI
            stat /= joint_entropy

        stats[pair] = stat
    return positions, stats


def _gen_combinations(num_pos: int, chunk_size: int) -> typing.Iterator[numpy.ndarray]:
    combs = itertools.combinations(range(num_pos), 2)

    while True:
        if chunk := list(itertools.islice(combs, chunk_size)):
            yield numpy.array(chunk)
        else:
            break


class calc_mi:
    """calculator for mutual information or normalised mutual information

    callable with positions to calculate the statistic for.
    """

    def __init__(
        self,
        data: numpy.ndarray,
        num_states: int,
        metric_id: enum.Enum,
    ) -> None:
        """
        Parameters
        ----------
        data
            a 2D numpy array of uint8 values representing a multiple sequence
            alignment where sequences are the first dimension and positions are
            the second dimension.
        num_states
            the number of canonical states in the moltype. Sequence elements that
            exceed this value are not included in the calculation.
        metric_id
            either 1 (for MI) or 2 (for NMI)
        """
        self._data = data
        self._num_states = num_states
        self._metric_id = metric_id
        self._canonical_pos = numpy.all(self._data < self._num_states, axis=0)
        self._entropies = _calc_column_entropies(self._data, self._num_states)

    def __call__(self, positions: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
        return _general_mi_calc(
            positions,
            self._canonical_pos,
            self._data,
            self._entropies,
            self._num_states,
            metric_id=self._metric_id,
        )


class calc_rmi:
    """calculator for resampled mutual information

    Callable with positions to calculate the statistic for.
    When called, it returns the positions and their corresponding statistic.
    """

    def __init__(self, data: numpy.ndarray, num_states: int) -> None:
        """
        Parameters
        ----------
        data
            a 2D numpy array of uint8 values representing a multiple sequence
            alignment where sequences are the first dimension and positions are
            the second dimension.
        num_states
            the number of canonical states in the moltype. Sequence elements that
            exceed this value are not included in the calculation.
        """
        self._data = data
        self._num_states = num_states

    def __call__(self, positions: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
        return _rmi_calc(positions, self._data, self._num_states)


@UI.display_wrap
def coevolution_matrix(
    *,
    alignment: "Alignment",
    positions: list[int] | None = None,
    stat: str = "nmi",
    parallel: bool = False,
    par_kw: dict | None = None,
    show_progress: bool = False,
    ui=None,
) -> dict_array.DictArray:
    """measure pairwise coevolution

    Parameters
    ----------
    aln
        sequence alignment
    stat
        either 'mi' (mutual information), 'nmi' (normalised MI) or 'rmi' (resampled MI)
    parallel
        run in parallel on your machine
    par_kw
        providing {'max_workers': 6} defines the number of workers to use, see
        arguments for cogent3.util.parallel.as_completed()
    show_progress
        displays a progress bar

    Returns
    -------
    Returns a DictArray with the pairwise coevolution values as a lower triangle. The other
    values are nan.
    """
    stat = {"mi": 1, "nmi": 2, "rmi": 3}[MI_METHODS(stat).name]
    num_states = len(alignment.moltype.alphabet)

    if hasattr(alignment, "__array__"):
        # new_type Alignment classes will support direct conversion
        # to numpy.uint8 arrays
        data = numpy.array(alignment)
    else:
        data = alignment.to_type(array_align=True).array_seqs

    if positions:
        positions = list(itertools.chain(*positions))
        data = data[:, tuple(positions)]
    else:
        positions = range(data.shape[1])

    num_pos = data.shape[1]

    calc = calc_rmi(data, num_states) if stat == 3 else calc_mi(data, num_states, stat)

    mutual_info = numpy.empty((num_pos, num_pos), dtype=numpy.float64)
    mutual_info.fill(numpy.nan)

    # we generate the positions as a numpy.array of tuples
    chunk_size = 10_000
    num_chunks = num_pos * (num_pos - 1) // 2 // chunk_size

    position_combinations = _gen_combinations(num_pos, chunk_size)
    if parallel:
        par_kw = par_kw or {}
        to_do = PAR.as_completed(calc, position_combinations, **par_kw)
    else:
        to_do = map(calc, position_combinations)

    for pos_pairs, stats in ui.series(
        to_do,
        noun="Sets of pairwise positions",
        count=num_chunks + 1,
    ):
        indices = numpy.ravel_multi_index(pos_pairs.T[::-1], (num_pos, num_pos))
        mutual_info.put(indices, stats)

    positions = list(positions)
    return dict_array.DictArray.from_array_names(mutual_info, positions, positions)

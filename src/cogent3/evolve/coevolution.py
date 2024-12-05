import enum
import itertools
import typing
from os.path import basename
from pickle import Pickler, Unpickler
from random import shuffle

import numba
import numpy
from numpy import (
    array,
    e,
    greater_equal,
    isnan,
    less_equal,
    log,
    nan,
    nonzero,
    ones,
    ravel,
    zeros,
)
from numpy.linalg import norm

from cogent3 import PROTEIN, make_aligned_seqs
from cogent3.core import new_alphabet
from cogent3.core.alphabet import CharAlphabet
from cogent3.core.moltype import IUPAC_gap, IUPAC_missing
from cogent3.core.sequence import Sequence
from cogent3.evolve.substitution_model import (
    EmpiricalProteinMatrix,
    Parametric,
)
from cogent3.maths.stats.distribution import binomial_exact
from cogent3.maths.stats.number import CategoryCounter, CategoryFreqs
from cogent3.maths.stats.special import ROUND_ERROR
from cogent3.util import dict_array
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI
from cogent3.util import warning as c3warn

DEFAULT_EXCLUDES = "".join([IUPAC_gap, IUPAC_missing])
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

    def __init__(self, data: numpy.ndarray, num_states: int, metric_id: enum.Enum):
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

    def __init__(self, data: numpy.ndarray, num_states: int):
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

    if stat == 3:
        calc = calc_rmi(data, num_states)
    else:
        calc = calc_mi(data, num_states, stat)

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


_reason_1 = (
    "discontinuing support for non-MI based calculations, use coevolution_matrix()"
)
_reason_2 = "discontinuing support for one-off calculations, use coevolution_matrix()"
_reason_3 = "discontinuing support due to poor performance, use coevolution_matrix()"


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def build_rate_matrix(
    count_matrix,
    freqs,
    aa_order="ACDEFGHIKLMNPQRSTVWY",
):  # pragma: no cover
    epm = EmpiricalProteinMatrix(count_matrix, freqs)
    word_probs = array([freqs[aa] for aa in aa_order])
    num = word_probs.shape[0]
    mprobs_matrix = ones((num, num), float) * word_probs

    return epm.calcQ(word_probs, mprobs_matrix)


# Mutual Information Analysis
# Mutual Information Calculators


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def mi(h1, h2, joint_h):  # pragma: no cover
    """Calc Mutual Information given two entropies and their joint entropy"""
    return h1 + h2 - joint_h


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def normalized_mi(h1, h2, joint_h):  # pragma: no cover
    """MI normalized by joint entropy, as described in Martin 2005"""
    return mi(h1, h2, joint_h) / joint_h


nmi = normalized_mi

# Other functions used in MI calculations


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def join_positions(pos1, pos2):  # pragma: no cover
    """Merge two positions and return as a list of strings

    pos1: iterable object containing the first positions data
    pos2: iterable object containing the second positions data

    Example:
        >>> join_positions("ABCD", "1234")
            ['A1', 'B2', 'C3', 'D4']
    """
    return ["".join([r1, r2]) for r1, r2 in zip(pos1, pos2, strict=False)]


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def joint_entropy(pos1, pos2):  # pragma: no cover
    """Calculate the joint entroy of a pair of positions"""
    return CategoryCounter(join_positions(pos1, pos2)).entropy


# Exclude handlers (functions for processing position strings with exclude
# characters)


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def ignore_excludes(pos, excludes=DEFAULT_EXCLUDES):  # pragma: no cover
    """Return position data as-is (results in excludes treated as other chars)"""
    return pos


# Functions for scoring coevolution on the basis of Mutual Information


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def mi_pair(
    alignment,
    pos1,
    pos2,
    h1=None,
    h2=None,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calculate mutual information of a pair of alignment positions

    alignment: the full alignment object
    pos1: index of 1st position in alignment to be compared
     (zero-based, not one-based)
    pos2: index of 2nd position in alignment to be compared
     (zero-based, not one-based)
    h1: entropy of pos1, if already calculated (to avoid time to recalc)
    h2: entropy of pos2, if already calculated (to avoid time to recalc)
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes position data and returns it
     with exclude characters processed in someway. Position data should be
     an iterable object containing the characters present at each position.
     f(position_data,excludes=gDefaultExcludes) -> position_data

    """
    col1 = list(alignment[pos1].positions)[0]
    col2 = list(alignment[pos2].positions)[0]
    # Detect and process exclude characters.
    # This bit of code is slow, and not necessary if
    # exclude_hanlder == ignore_excludes, so I explicitly
    # check, and bypass this block if possible.
    if exclude_handler != ignore_excludes:
        for col in (col1, col2):
            states = set(col)
            for exclude in excludes:
                if exclude in states:
                    try:
                        _ = exclude_handler(col, excludes)
                        break
                    except TypeError:
                        return null_value

    # Calculate entropy of pos1 & pos2, if they weren't passed in.
    if not h1:
        h1 = CategoryCounter(col1).entropy
    if not h2:
        h2 = CategoryCounter(col2).entropy
    # Calculate the joint entropy of pos1 & pos2
    joint_h = joint_entropy(col1, col2)
    # Calculate MI using the specified method -- return null_value when
    # the specified MI cannot be calculated
    # (e.g., mi_calculator=nmi and joint_h=0.0)
    try:
        result = mi_calculator(h1, h2, joint_h)
        if result <= ROUND_ERROR:
            result = 0.0
    except ZeroDivisionError:
        result = null_value
    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def mi_position(
    alignment,
    position,
    positional_entropies=None,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calc mi b/w position and all other positions in an alignment

    alignment: the full alignment object
    position: the position number of interest -- NOTE: this is the
     position index, not the sequenece position (so zero-indexed, not
    one-indexed)
    positional_entropies: a list containing the entropy of each position in
     the alignment -- these can be passed in to avoid recalculating if
     calling this function over more than one position (e.g., in
     mi_alignment)
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    aln_length = len(alignment)
    # Create result vector
    result = zeros(aln_length, float)

    # compile positional entropies if not passed in
    if positional_entropies is None:
        positional_entropies = [CategoryCounter(p).entropy for p in alignment.positions]

    # Will want to make a change here so that we don't need to recalculate
    # all values when calling from mi_alignment
    for i in range(aln_length):
        result[i] = mi_pair(
            alignment,
            pos1=position,
            pos2=i,
            h1=positional_entropies[position],
            h2=positional_entropies[i],
            mi_calculator=mi_calculator,
            null_value=null_value,
            excludes=excludes,
            exclude_handler=exclude_handler,
        )
    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def mi_alignment(
    alignment,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calc mi over all position pairs in an alignment

    alignment: the full alignment object
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    aln_length = len(alignment)
    # Create result matrix
    result = zeros((aln_length, aln_length), float)

    # Compile postional entropies for each position in the alignment
    # I believe I started using this rather than alignment.uncertainties
    # b/c the latter relies on converting a ArrayAlignment to an Alignment --
    # need to check into this.
    positional_entropies = alignment.entropy_per_pos()

    # Calculate pairwise MI between position_number and all alignment
    # positions, and return the results in a vector.
    for i in range(aln_length):
        for j in range(i + 1):
            result[i, j] = mi_pair(
                alignment,
                pos1=i,
                pos2=j,
                h1=positional_entropies[i],
                h2=positional_entropies[j],
                mi_calculator=mi_calculator,
                null_value=null_value,
                excludes=excludes,
                exclude_handler=exclude_handler,
            )
    # copy the lower triangle to the upper triangle to make
    # the matrix symmetric
    ltm_to_symmetric(result)
    return result


# End Mutual Information Analysis

# Start Normalized Mutual Information Analysis (Martin 2005)


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def normalized_mi_pair(
    alignment,
    pos1,
    pos2,
    h1=None,
    h2=None,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calc normalized mutual information of a pair of alignment positions

    alignment: the full alignment object
    pos1: index of 1st position in alignment to be compared
     (zero-based, not one-based)
    pos2: index of 2nd position in alignment to be compared
     (zero-based, not one-based)
    h1: entropy of pos1, if already calculated (to avoid time to recalc)
    h2: entropy of pos2, if already calculated (to avoid time to recalc)
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    return mi_pair(
        alignment,
        pos1,
        pos2,
        h1=h1,
        h2=h2,
        mi_calculator=nmi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_pair = normalized_mi_pair


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def normalized_mi_position(
    alignment,
    position,
    positional_entropies=None,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calc normalized mi b/w position and all other positions in an alignment

    alignment: the full alignment object
    position: the position number of interest -- NOTE: this is the
     position index, not the sequenece position (so zero-indexed, not
    one-indexed)
    positional_entropies: a list containing the entropy of each position in
     the alignment -- these can be passed in to avoid recalculating if
     calling this function over more than one position (e.g., in
     mi_alignment)
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    return mi_position(
        alignment,
        position,
        positional_entropies=positional_entropies,
        mi_calculator=nmi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_position = normalized_mi_position


@c3warn.deprecated_callable("2024.12", reason=_reason_2, is_discontinued=True)
def normalized_mi_alignment(
    alignment,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):  # pragma: no cover
    """Calc normalized mi over all position pairs in an alignment

    alignment: the full alignment object
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.
    """
    return mi_alignment(
        alignment=alignment,
        mi_calculator=normalized_mi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_alignment = normalized_mi_alignment
# End Normalized Mutual Information Analysis


# Start Statistical coupling analysis (SCA) (Suel 2003)
class SCAError(Exception):
    pass


# PROTEIN's alphabet contains U, so redefining the alphabet for now
# rather than use PROTEIN.alphabet. May want to revist this decision...
AAGapless = CharAlphabet("ACDEFGHIKLMNPQRSTVWY")
default_sca_alphabet = AAGapless
# AAGapless = PROTEIN.alphabet

# Dictionary of mean AA-frequencies in all natural proteins
# Compiled by Rama Ranganathan from 36,498 unique eukaryotic proteins
# from the Swiss-Prot database
protein_dict = {
    "A": 0.072658,
    "C": 0.024692,
    "D": 0.050007,
    "E": 0.061087,
    "F": 0.041774,
    "G": 0.071589,
    "H": 0.023392,
    "I": 0.052691,
    "K": 0.063923,
    "L": 0.089093,
    "M": 0.02315,
    "N": 0.042931,
    "P": 0.052228,
    "Q": 0.039871,
    "R": 0.052012,
    "S": 0.073087,
    "T": 0.055606,
    "V": 0.063321,
    "W": 0.01272,
    "Y": 0.032955,
}
default_sca_freqs = protein_dict


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def freqs_to_array(f, alphabet):  # pragma: no cover
    """Takes data in freqs object and turns it into array.

    f = dict or CategoryCounter object
    alphabet = Alphabet object or just a list that specifies the order
        of things to appear in the resulting array
    """
    return array([f.get(i, 0) for i in alphabet])


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def get_allowed_perturbations(
    counts,
    cutoff,
    alphabet,
    num_seqs=100,
):  # pragma: no cover
    """Returns list of allowed perturbations as characters

    count: Profile object of raw character counts at each position
    num_seqs: number of sequences in the alignment
    cutoff: minimum number of sequences in the subalignment (as fraction
    of the total number of seqs in the alignment.

    A perturbation is allowed if the subalignment of sequences that
    contain the specified char at the specified position is larger
    that the cutoff value * the total number of sequences in the alignment.

    """
    result = []
    abs_cutoff = cutoff * num_seqs

    for char, count in zip(alphabet, counts, strict=False):
        if count >= abs_cutoff:
            result.append(char)
    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def probs_from_dict(d, alphabet):  # pragma: no cover
    """Convert dict of alphabet char probabilities to list in alphabet's order

    d: probabilities of observing each character in alphabet (dict indexed
     by char)
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)

    """
    return array([d[c] for c in alphabet])


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def freqs_from_aln(aln, alphabet, scaled_aln_size=100):  # pragma: no cover
    """Return the frequencies in aln of chars in alphabet's order

    aln: the alignment object
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    *Warning: characters in aln that are not in alphabet are silently
        ignored. Is this the desired behavior?

    Need to combine this function with get_position_frequences (and renamed
     that one to be more generic) since they're doing the same thing now.

    """
    alphabet_as_indices = array([aln.alphabet.to_indices(alphabet)]).transpose()
    aln_data = ravel(aln.array_positions)
    return (alphabet_as_indices == aln_data).sum(1) * (scaled_aln_size / len(aln_data))


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def get_positional_frequencies(
    aln,
    position_number,
    alphabet,
    scaled_aln_size=100,
):  # pragma: no cover
    """Return the freqs in aln[position_number] of chars in alphabet's order

    aln: the alignment object
    position_number: the index of the position of interest in aln
     (note: zero-based alignment indexing)
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    *Warning: characters in aln that are not in alphabet are silently
        ignored. Is this the desired behavior?

    """
    alphabet_as_indices = array([aln.alphabet.to_indices(alphabet)]).transpose()
    position_data = aln.array_positions[position_number]
    return (alphabet_as_indices == position_data).sum(1) * (
        scaled_aln_size / len(position_data)
    )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def get_positional_probabilities(
    pos_freqs,
    natural_probs,
    scaled_aln_size=100,
):  # pragma: no cover
    """Get probs of observering the freq of each char given it's natural freq
    In Suel 2003 supplementary material, this step is defined as:
     "... each element is the binomial probability of observing each
      amino acid residue at position j given its mean frequency in
      all natural proteins."
    This function performs the calculate for a single position.

    pos_freqs: the frequencies of each char in the alphabet at a
     position-of-interest in the alignment (list of floats, typically
     output of get_positional_frequencies)
    natural_probs: the natural probabilities of observing each char
     in the alphabet (list of floats: typically output of probs_from_dict)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    Note: It is critical that the values in pos_freqs and natural_probs are
     in the same order, which should be the order of chars in the alphabet.

    """
    results = []
    for pos_freq, natural_prob in zip(pos_freqs, natural_probs, strict=False):
        try:
            results.append(binomial_exact(pos_freq, scaled_aln_size, natural_prob))
        # Because of the scaling of alignments to scaled_aln_size, pos_freq is
        # a float rather than an int. So, if a position is perfectly conserved,
        # pos_freq as a float could be greater than scaled_aln_size.
        # In this case I cast it to an int. I don't like this alignment
        # scaling stuff though.
        except ValueError:
            results.append(binomial_exact(int(pos_freq), scaled_aln_size, natural_prob))
    return array(results)


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def get_subalignments(aln, position, selections):  # pragma: no cover
    """returns subalns w/ seq[pos] == selection for each in selections
    aln: an alignment object
    position: int in alignment to be checked for each perturbation
    selections: characters which must be present at seq[pos] for
        seq to be in subalignment

    Note: This method returns a list of subalignments corresponding
        to the list of selections. So, if you specify selections as
        ['A','G'], you would get two subalignments back -- the first
        containing sequences with 'A' at position, and the second
        containing sequences with 'G' at position. If you want all
        sequences containing either 'A' or 'G', merge the resulting
        subalignments.

    """
    result = []
    for s in aln.alphabet.to_indices(selections):
        seqs_to_keep = nonzero(aln.array_seqs[:, position] == s)[0]
        sub_align = aln.get_sub_alignment(seqs=seqs_to_keep)
        if sub_align is not None:
            result.append(sub_align)
    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def get_dg(position_probs, aln_probs):  # pragma: no cover
    """Return delta_g vector

    position_probs: the prob of observing each alphabet chars frequency in
     the alignment position-of-interest, given it's background frequency
     in all proteins (list of floats, typically the output of
     get_positional_probabilities)
    aln_probs: the prob of observing each alphabet chars frequency in the
     full alignment, given it's background frequency (list of floats)

    """
    results = []
    for position_prob, aln_prob in zip(position_probs, aln_probs, strict=False):
        results.append(log(position_prob / aln_prob))
    return array(results)


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def get_dgg(all_dgs, subaln_dgs, scaled_aln_size=100):  # pragma: no cover
    """Return delta_delta_g value

    all_dgs: the dg vector for a position-of-interest in the alignment
     (list of floats, typically the output of get_dg)
    subaln_dgs: the dg vector for a sub-alignment of the position-of-
     interest in the alignment (list of floats, typically the output
     of get_dg applied to a sub-alignment)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    * There are two weird issues in this function with respect to the
    desciption of the algorithm in the Suel 2003 supplementary material.
    In order to get the values presented in their GPCR paper, we need to
    (1) divide the euclidian norm by the scaled_aln_size, and then (2)
    multiply the result by e.
    ** IT IS CRITICAL TO UNDERSTAND WHY
    WE NEED TO APPLY THESE STEPS BEFORE PUBLISHING ANYTHING THAT USES
    THIS CODE.**

    * A possible reason for the mysterious e scaling is that we are
    misinterpreting what they mean when they say ddg is 'the magnitude of
    this difference vector.' We are assuming they are referring to the
    Euclidian norm, but until I see their code, I can't be sure about
    this.
    """
    return norm(all_dgs - subaln_dgs) / scaled_aln_size * e


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def sca_pair(
    alignment,
    pos1,
    pos2,
    cutoff,
    position_freqs=None,
    position_probs=None,
    dgs=None,
    perturbations=None,
    scaled_aln_size=100,
    null_value=DEFAULT_NULL_VALUE,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):  # pragma: no cover
    """ Calculate statistical coupling b/w a pair of alignment columns

        alignment: full alignment object
        pos1: the first position used to probe for statistical coupling
         (subalignments will be generated based on allowed perturbations
         at this position) -- int, zero-based indexing into alignment
        pos2: the second position used to probe for statistical coupling
         -- int, zero-based indexing into alignment
        cutoff: the percentage of sequences that must contain a specific
         char at a specific pos1 to result in an allowed sub-alignment.
         (According to the Ranganathan papers, this should be the value
         determined by their 3rd criteria.)
        position_freqs: if precalculated, a matrix containing the output
         of get_positional_frequencies for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        position_probs: if precalculated, a matrix containing the output
         of get_positional_probabilities for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        dgs: if precalculated, a matrix containing the output
         of get_dg for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        perturbations: if precalculated, a matrix containing the output
         of get_allowed_perturbations for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        scaled_aln_size: the scaled number of sequences in the alignment. The
         original SCA implementation treats all alignments as if they contained
         100 sequences when calculating frequencies and probabilities. 100 is
         therefore the default value.
        null_value: the value which should be returned if SCA cannot or
         should not be calculated (e.g., no allowed perturbations or
         pos1==pos2, respectively).
        return_all: if cutoff <= 0.50, it is possible that there will be more
         than one allowed_perturbation per position. In these cases, either all
         of the values could be returned (return_all=True) or the max of the
         values can be returned (return_all=False, default). If you'd like one
         value, but not the max, wrap this function with return_all=True, and
         handle the return value as desired.
        alphabet: an ordered iterable object containing the characters in the
         alphabet. For example, this can be a CharAlphabet object, a list,
         or a string.

        **IMPORTANT NOTE: SCA, unlike (all?) other methods implemented here,
         requires the full alignment, even to calculate coupling between just
         a pair of positions. Because frequencies of characters in the full
         alignment are compared with frequencies at each position, you cannot
         simply pull out two columns of the alignment, and pass them to this
         function as a subalignment. Your results would differ from calculating
         coupling of the same positions with the full alignment. For example:
            sca_pair(aln,10,20,0.85) != \
            sca_pair(aln.take_positions([10,20]),0,1,0.85)
    """

    # Calculate frequency distributions
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)

    # get positional frequencies
    if position_freqs:
        pos1_freqs = position_freqs[pos1]
        pos2_freqs = position_freqs[pos2]
    else:
        pos1_freqs = get_positional_frequencies(
            alignment,
            pos1,
            alphabet,
            scaled_aln_size,
        )
        pos2_freqs = get_positional_frequencies(
            alignment,
            pos2,
            alphabet,
            scaled_aln_size,
        )
    # get positional probability vectors ("... each element is the binomial
    # probability of observing each amino acid residue at position j given its
    # mean frequency in all natural proteins." Suel 2003 supplementary
    # material)
    if position_probs:
        pos2_probs = position_probs[pos2]
    else:
        pos2_probs = get_positional_probabilities(
            pos2_freqs,
            natural_probs,
            scaled_aln_size,
        )

    # get statistical energies for pos2 in full alignment
    if dgs:
        pos2_dg = dgs[pos2]
    else:
        pos2_dg = get_dg(pos2_probs, aln_probs)

    # determine allowed perturbations
    if perturbations:
        allowed_perturbations = perturbations[pos1]
    else:
        allowed_perturbations = get_allowed_perturbations(
            pos1_freqs,
            cutoff,
            alphabet,
            scaled_aln_size,
        )
    # should we do something different here on return_all == True?
    if not allowed_perturbations:
        return null_value

    # generate the subalignments which contain each allowed
    # perturbation residue at pos1
    subalignments = get_subalignments(alignment, pos1, allowed_perturbations)

    # calculate ddg for each allowed perturbation
    ddg_values = []
    for subalignment in subalignments:
        # Calculate dg for the subalignment
        subaln_freqs = freqs_from_aln(subalignment, alphabet, scaled_aln_size)
        subaln_probs = get_positional_probabilities(
            subaln_freqs,
            natural_probs,
            scaled_aln_size,
        )
        subaln_pos2_freqs = get_positional_frequencies(
            subalignment,
            pos2,
            alphabet,
            scaled_aln_size,
        )
        subaln_pos2_probs = get_positional_probabilities(
            subaln_pos2_freqs,
            natural_probs,
            scaled_aln_size,
        )
        subaln_dg = get_dg(subaln_pos2_probs, subaln_probs)
        ddg_values.append(get_dgg(pos2_dg, subaln_dg, scaled_aln_size))

    if return_all:
        return list(zip(allowed_perturbations, ddg_values, strict=False))
    return max(ddg_values)


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def sca_position(
    alignment,
    position,
    cutoff,
    position_freqs=None,
    position_probs=None,
    dgs=None,
    perturbations=None,
    scaled_aln_size=100,
    null_value=DEFAULT_NULL_VALUE,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):  # pragma: no cover
    """Calculate statistical coupling b/w a column and all other columns

    alignment: full alignment object
    position: the position of interest to probe for statistical coupling
     (subalignments will be generated based on allowed perturbations
     at this position) -- int, zero-based indexing into alignment
    cutoff: the percentage of sequences that must contain a specific
     char at a specific pos1 to result in an allowed sub-alignment.
     (According to the Ranganathan papers, this should be the value
     determined by their 3rd criteria.)
    position_freqs: if precalculated, a matrix containing the output
     of get_positional_frequencies for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    position_probs: if precalculated, a matrix containing the output
     of get_positional_probabilities for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    dgs: if precalculated, a matrix containing the output
     of get_dg for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    perturbations: if precalculated, a matrix containing the output
     of get_allowed_perturbations for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.
    null_value: the value which should be returned if SCA cannot or
     should not be calculated (e.g., no allowed perturbations or
    pos1==pos2, respectively).
    return_all: if cutoff <= 0.50, it is possible that there will be more
     than one allowed_perturbation per position. In these cases, either all
     of the values could be returned (return_all=True) or the max of the
     values can be returned (return_all=False, default). If you'd like one
     value, but not the max, wrap this function with return_all=True, and
     handle the return value as desired.
    alphabet: an ordered iterable object containing the characters in the
     alphabet. For example, this can be a CharAlphabet object, a list,
     or a string.

    """
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)
    if not position_freqs:
        position_freqs = []
        for i in range(len(alignment)):
            position_freqs.append(
                get_positional_frequencies(alignment, i, alphabet, scaled_aln_size),
            )

    if not position_probs:
        position_probs = []
        for i in range(len(alignment)):
            position_probs.append(
                get_positional_probabilities(
                    position_freqs[i],
                    natural_probs,
                    scaled_aln_size,
                ),
            )
    if not dgs:
        dgs = []
        for i in range(len(alignment)):
            dgs.append(get_dg(position_probs[i], aln_probs))

    if not perturbations:
        perturbations = []
        for i in range(len(alignment)):
            perturbations.append(
                get_allowed_perturbations(
                    position_freqs[i],
                    cutoff,
                    alphabet,
                    scaled_aln_size,
                ),
            )

    result = []
    for i in range(len(alignment)):
        result.append(
            sca_pair(
                alignment,
                position,
                i,
                cutoff,
                position_freqs=position_freqs,
                position_probs=position_probs,
                dgs=dgs,
                perturbations=perturbations,
                scaled_aln_size=scaled_aln_size,
                null_value=null_value,
                return_all=return_all,
                alphabet=alphabet,
                background_freqs=background_freqs,
            ),
        )
    return array(result)


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def sca_alignment(
    alignment,
    cutoff,
    null_value=DEFAULT_NULL_VALUE,
    scaled_aln_size=100,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):  # pragma: no cover
    """Calculate statistical coupling b/w all columns in alignment

    alignment: full alignment object
    cutoff: the percentage of sequences that must contain a specific
     char at a specific pos1 to result in an allowed sub-alignment.
     (According to the Ranganathan papers, this should be the value
     determined by their 3rd criteria.)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.
    null_value: the value which should be returned if SCA cannot or
     should not be calculated (e.g., no allowed perturbations or
    pos1==pos2, respectively).
    return_all: if cutoff <= 0.50, it is possible that there will be more
     than one allowed_perturbation per position. In these cases, either all
     of the values could be returned (return_all=True) or the max of the
     values can be returned (return_all=False, default). If you'd like one
     value, but not the max, wrap this function with return_all=True, and
     handle the return value as desired.
    alphabet: an ordered iterable object containing the characters in the
     alphabet. For example, this can be a CharAlphabet object, a list,
     or a string.

    """
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)
    # get all positional frequencies
    position_freqs = []
    for i in range(len(alignment)):
        position_freqs.append(
            get_positional_frequencies(alignment, i, alphabet, scaled_aln_size),
        )

    # get all positional probabilities
    position_probs = []
    for i in range(len(alignment)):
        position_probs.append(
            get_positional_probabilities(
                position_freqs[i],
                natural_probs,
                scaled_aln_size,
            ),
        )

    # get all delta_g vectors
    dgs = []
    for i in range(len(alignment)):
        dgs.append(get_dg(position_probs[i], aln_probs))

    # get all allowed perturbations
    perturbations = []
    for i in range(len(alignment)):
        perturbations.append(
            get_allowed_perturbations(
                position_freqs[i],
                cutoff,
                alphabet,
                scaled_aln_size,
            ),
        )

    result = []
    for i in range(len(alignment)):
        result.append(
            sca_position(
                alignment,
                i,
                cutoff,
                position_freqs=position_freqs,
                position_probs=position_probs,
                dgs=dgs,
                perturbations=perturbations,
                scaled_aln_size=scaled_aln_size,
                null_value=null_value,
                return_all=return_all,
                alphabet=alphabet,
                background_freqs=background_freqs,
            ),
        )
    return array(result)


# End statistical coupling analysis

# Start Resampled Mutual Information Analysis
# (developed by Hutley and Easton, and first published in
# Caporaso et al., 2008)


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def make_weights(counts, n):  # pragma: no cover
    """Return the weights for replacement states for each possible character.
    We compute the weight as the normalized frequency of the replacement state
    divided by 2*n."""
    char_prob = list(counts.to_freqs().items())
    weights = []
    for C, P in char_prob:
        alts = CategoryFreqs({c: p for c, p in char_prob if c != C})
        alts = alts.to_normalized()
        alts = CategoryCounter({c: w / (2 * n) for c, w in list(alts.items())})
        weights += [(C, alts)]
    return weights


def calc_pair_scale(seqs, obs1, obs2, weights1, weights2):  # pragma: no cover
    """Return entropies and weights for comparable alignment.
    A comparable alignment is one in which, for each paired state ij, all
    alternate observable paired symbols are created. For instance, let the
    symbols {A,C} be observed at position i and {A,C} at position j. If we
    observe the paired types {AC, AA}. A comparable alignment would involve
    replacing an AC pair with a CC pair."""
    # scale is calculated as the product of mi from col1 with alternate
    # characters. This means the number of states is changed by swapping
    # between the original and selected alternate, calculating the new mi

    pair_freqs = CategoryCounter(seqs)
    weights1 = dict(weights1)
    weights2 = dict(weights2)
    scales = []
    for a, b in list(pair_freqs.keys()):
        weights = weights1[a]

        pr = a + b
        pair_freqs -= pr
        obs1 -= a

        # make comparable alignments by mods to col 1
        for c, w in list(weights.items()):
            new_pr = c + b
            pair_freqs += new_pr
            obs1 += c

            entropy = mi(obs1.entropy, obs2.entropy, pair_freqs.entropy)
            scales += [(pr, entropy, w)]

            pair_freqs -= new_pr
            obs1 -= c

        obs1 += a
        # make comparable alignments by mods to col 2
        weights = weights2[b]
        obs2 -= b
        for c, w in list(weights.items()):
            new_pr = a + c
            pair_freqs += new_pr
            obs2 += c

            entropy = mi(obs1.entropy, obs2.entropy, pair_freqs.entropy)
            scales += [(pr, entropy, w)]

            obs2 -= c
            pair_freqs -= new_pr

        obs2 += b

        pair_freqs += pr
    return scales


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def resampled_mi_pair(
    alignment,
    pos1,
    pos2,
    weights=None,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    """returns scaled mutual information for a pair.

    Parameters
    ----------
    alignment
        Alignment instance
    pos1, pos2
        alignment positions to be assessed
    weights
        Freq objects of weights for pos1, pos2
    excludes
        states to be excluded.

    """
    positions = list(alignment.positions)
    col1 = positions[pos1]
    col2 = positions[pos2]
    seqs = ["".join(p) for p in zip(col1, col2, strict=False)]
    for col in (col1, col2):
        states = {}.fromkeys(col)
        for exclude in excludes:
            if exclude in states:
                try:
                    _ = exclude_handler(col, excludes)
                    break
                except TypeError:
                    return null_value

    excludes = excludes or []
    num = len(seqs)
    col1 = CategoryCounter(col1)
    col2 = CategoryCounter(col2)
    seq_freqs = CategoryCounter(seqs)
    if weights:
        weights1, weights2 = weights
    else:
        weights1 = make_weights(col1, num)
        weights2 = make_weights(col2, num)

    entropy = mi(col1.entropy, col2.entropy, seq_freqs.entropy)
    scales = calc_pair_scale(seqs, col1, col2, weights1, weights2)
    scaled_mi = 1 - sum([w * seq_freqs[pr] for pr, e, w in scales if entropy <= e])

    return scaled_mi


rmi_pair = resampled_mi_pair


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def resampled_mi_position(
    alignment,
    position,
    positional_entropies=None,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    aln_length = len(alignment)
    result = zeros(aln_length, float)

    if positional_entropies is None:
        positional_entropies = alignment.entropy_per_pos()

    for i in range(aln_length):
        result[i] = resampled_mi_pair(
            alignment,
            pos1=position,
            pos2=i,
            excludes=excludes,
            exclude_handler=exclude_handler,
            null_value=null_value,
        )
    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def resampled_mi_alignment(
    alignment,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    """returns scaled mutual information for all possible pairs."""
    aln_length = len(alignment)
    result = zeros((aln_length, aln_length), float)
    positional_entropies = alignment.entropy_per_pos()

    for i in range(aln_length):
        result[i] = resampled_mi_position(
            alignment=alignment,
            position=i,
            positional_entropies=positional_entropies,
            excludes=excludes,
            exclude_handler=exclude_handler,
            null_value=null_value,
        )
    return result


# End Resampled Mutual Information Analysis

# Begin ancestral_states analysis


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def get_ancestral_seqs(
    aln,
    tree,
    sm=None,
    pseudocount=1e-6,
    optimise=True,
):  # pragma: no cover
    """Calculates ancestral sequences by maximum likelihood

    Parameters
    ----------
    sm
        a Parametric instance. If not provided, one is
        constructed from the alignment alphabet
    pseudocount
        unobserved sequence states must not be zero, this value
        is assigned to sequence states not observed in the alignment.
    optimise
        whether to optimise the likelihood function.

        Note: for the sake of reduced alphabets, we calculate the
         substitution model from the alignment. This also appears
         to be what what described in Tuffery 2000, although they're
         not perfectly clear about it.
    """
    from cogent3.core.alignment import ArrayAlignment

    sm = sm or Parametric(aln.alphabet, recode_gaps=True)
    lf = sm.make_likelihood_function(tree, sm.motif_probs)
    lf.set_alignment(aln, motif_pseudocount=pseudocount)
    if optimise:
        lf.optimise(local=True, show_progress=False)
    return ArrayAlignment(lf.likely_ancestral_seqs(), moltype=aln.moltype)


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def ancestral_state_alignment(
    aln,
    tree,
    ancestral_seqs=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    result = []
    for i in range(len(aln)):
        row = [null_value] * len(aln)
        for j in range(i + 1):
            row[j] = ancestral_state_pair(aln, tree, i, j, ancestral_seqs, null_value)
        result.append(row)
    return ltm_to_symmetric(array(result))


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def ancestral_state_position(
    aln,
    tree,
    position,
    ancestral_seqs=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    result = []
    for i in range(len(aln)):
        result.append(
            ancestral_state_pair(aln, tree, position, i, ancestral_seqs, null_value),
        )
    return array(result)


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def ancestral_state_pair(
    aln,
    tree,
    pos1,
    pos2,
    ancestral_seqs=None,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    """ """
    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    ancestral_names_to_seqs = dict(
        list(zip(ancestral_seqs.names, ancestral_seqs.array_seqs, strict=False)),
    )
    distances = tree.get_distances()
    tips = tree.get_node_names(tipsonly=True)
    # map names to nodes (there has to be a built-in way to do this
    # -- what is it?)
    nodes = dict([(n, tree.get_node_matching_name(n)) for n in tips])
    # add tip branch lengths as distance b/w identical tips -- this is
    # necessary for my weighting step, where we want correlated changes
    # occuring on a single branch to be given the most weight
    distances.update(dict([((n, n), nodes[n].length) for n in nodes]))
    result = 0
    names_to_seqs = dict(list(zip(aln.names, aln.array_seqs, strict=False)))
    for i in range(len(tips)):
        org1 = tips[i]
        seq1 = names_to_seqs[org1]
        for j in range(i, len(tips)):
            org2 = tips[j]
            seq2 = names_to_seqs[org2]
            ancestor = nodes[org1].last_common_ancestor(nodes[org2]).name
            if ancestor == org1 == org2:
                # we're looking for correlated change along a
                # single branch
                ancestral_seq = ancestral_names_to_seqs[nodes[org1].ancestors()[0].name]
            else:
                # we're looking for correlated change along different
                # branches (most cases)
                ancestral_seq = ancestral_names_to_seqs[ancestor]

            # get state of pos1 in org1, org2, and ancestor
            org1_p1 = seq1[pos1]
            org2_p1 = seq2[pos1]
            ancestor_p1 = ancestral_seq[pos1]

            # if pos1 has changed in both organisms since their lca,
            # this is a position of interest
            if org1_p1 != ancestor_p1 and org2_p1 != ancestor_p1:
                # get state of pos2 in org1, org2, and ancestor
                org1_p2 = seq1[pos2]
                org2_p2 = seq2[pos2]
                ancestor_p2 = ancestral_seq[pos2]
                # if pos2 has also changed in both organisms since their lca,
                # then we add a count for a correlated change
                if org1_p2 != ancestor_p2 and org2_p2 != ancestor_p2:
                    # There are a variety of ways to score. The simplest is
                    # to increment by one, which seems to be what was done
                    # in other papers.) This works well, but in a quick test
                    # (alpha helices/myoglobin with several generally
                    # high scoring alphabets) weighting works better. A more
                    # detailed analysis is in order.
                    # result += 1
                    # Now I weight based on distance so
                    # changes in shorter time are scored higher than
                    # in longer time. (More ancient changes
                    # are more likely to be random than more recent changes,
                    # b/c more time has passed for the changes to occur in.)
                    # This gives results
                    # that appear to be better under some circumstances,
                    # and at worst, about the same as simply incrementing
                    # by 1.
                    result += 1 / distances[(org1, org2)]
                    # Another one to try might involve discounting the score
                    # for a pair when one changes and the other doesn't.
    return result


# End ancestral_states analysis


# Methods for running coevolutionary analyses on sequence data.
method_abbrevs_to_names = {
    "mi": "Mutual Information",
    "nmi": "Normalized Mutual Information",
    "sca": "Statistical Coupling Analysis",
    "an": "Ancestral States",
    "rmi": "Resampled Mutual Information",
}

# Method-specific error checking functions
# Some of the coevolution algorithms require method-specific input validation,
# but that code isn't included in the alrogithm-specific functions (e.g.
# sca_alignment,
# sca_pair) because those are sometimes run many times. For example,
# sca_alignment makes many calls to sca_pair, so we can't have sca_pair
# perform validation every time it's called. My solution is to have the
# coevolve_* functions perform the input validation, and recommend that
# users always perform analyses via these functions. So, in the above example,
# the user would access sca_alignment via coevolve_alignment('sca', ...). Since
# sca_alignment makes calls to sca_pair, not coevolve_pair, the input
# validation
# is only performed once by coevolve_alignment.


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def sca_input_validation(alignment, **kwargs):  # pragma: no cover
    """SCA specific validations steps"""

    # check that all required parameters are present in kwargs
    required_parameters = ["cutoff"]
    # users must provide background frequencies for MolTypes other
    # than PROTEIN -- by default, protein background frequencies are used.
    if alignment.moltype != PROTEIN:
        required_parameters.append("background_freqs")
    for rp in required_parameters:
        if rp not in kwargs:
            raise ValueError("Required parameter was not provided: " + rp)

    # check that the value provided for cutoff is valid (ie. between 0 and 1)
    if not 0.0 <= kwargs["cutoff"] <= 1.0:
        raise ValueError("Cutoff must be between zero and one.")

    # check that the set of chars in alphabet and background_freqs are
    # identical
    try:
        alphabet = kwargs["alphabet"]
    except KeyError:
        # We want to use the PROTEIN alphabet minus the U character for
        # proteins since we don't have a background frequency for U
        if alignment.moltype == PROTEIN:
            alphabet = AAGapless
        else:
            alphabet = alignment.moltype.alphabet
    try:
        background_freqs = kwargs["background_freqs"]
    except KeyError:
        background_freqs = default_sca_freqs
    validate_alphabet(alphabet, background_freqs)


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def validate_alphabet(alphabet, freqs):  # pragma: no cover
    """SCA validation: ValueError if set(alphabet) != set(freqs.keys())"""
    alphabet_chars = set(alphabet)
    freq_chars = set(freqs.keys())
    if alphabet_chars != freq_chars:
        raise ValueError(
            "Alphabet and background freqs must contain identical sets of chars.",
        )


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def ancestral_states_input_validation(alignment, **kwargs):  # pragma: no cover
    """Ancestral States (AS) specific validations steps"""
    # check that all required parameters are present in kwargs
    required_parameters = ["tree"]
    for rp in required_parameters:
        if rp not in kwargs:
            raise ValueError("Required parameter was not provided: " + rp)

    # validate the tree
    validate_tree(alignment, kwargs["tree"])

    # if ancestral seqs are provided, validate them. (If calculated on the fly,
    # we trust them.)
    if "ancestral_seqs" in kwargs:
        validate_ancestral_seqs(alignment, kwargs["tree"], kwargs["ancestral_seqs"])


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def validate_ancestral_seqs(alignment, tree, ancestral_seqs):  # pragma: no cover
    """AS validation: ValueError if incompatible aln, tree, & ancestral seqs

    Incompatibility between the alignment and the ancestral_seqs is
        different sequence lengths. Incompatbility between the tree and
        the ancestral seqs is imperfect overlap between the names of the
        ancestors in the tree and the ancestral sequence names.
    """
    if len(alignment) != len(ancestral_seqs):
        raise ValueError("Alignment and ancestral seqs are different lengths.")
    # is there a better way to get all the ancestor names? why doesn't
    # tree.ancestors() do this?
    edges = set(tree.get_node_names()) - set(tree.get_tip_names())
    seqs = set(ancestral_seqs.names)
    if edges != seqs:
        raise ValueError(
            "Must be ancestral seqs for all edges and root in tree, and no more.",
        )


@c3warn.deprecated_callable("2024.12", reason=_reason_1, is_discontinued=True)
def validate_tree(alignment, tree):  # pragma: no cover
    """AS validation: ValueError if tip and seq names aren't same"""
    if set(tree.get_tip_names()) != set(alignment.names):
        raise ValueError("Tree tips and seqs must have perfectly overlapping names.")


# End method-specific error checking functions

# General (opposed to algorithm-specific) validation functions


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def validate_position(alignment, position):  # pragma: no cover
    """ValueError if position is outside the range of the alignment"""
    if not 0 <= position < len(alignment):
        raise ValueError(
            "Position is outside the range of the alignment: " + str(position),
        )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def validate_alignment(alignment):  # pragma: no cover
    """ValueError on ambiguous alignment characters"""
    bad_seqs = []
    for name, ambiguous_pos in list(alignment.get_ambiguous_positions().items()):
        if ambiguous_pos:
            bad_seqs.append(name)
    if bad_seqs:
        raise ValueError(
            f"Ambiguous characters in sequences: {'; '.join(map(str, bad_seqs))}",
        )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolve_alignments_validation(
    method,
    alignment1,
    alignment2,
    min_num_seqs,
    max_num_seqs,
    **kwargs,
):  # pragma: no cover
    """Validation steps required for intermolecular coevolution analyses"""
    valid_methods_for_different_moltypes = {}.fromkeys(
        [mi_alignment, nmi_alignment, resampled_mi_alignment],
    )
    if (
        alignment1.moltype != alignment2.moltype
    ) and method not in valid_methods_for_different_moltypes:
        raise AssertionError(
            "Different MolTypes only supported for %s"
            % " ".join(map(str, list(valid_methods_for_different_moltypes.keys()))),
        )

    alignment1_names = set([n.split("+")[0].strip() for n in alignment1.names])
    alignment2_names = set([n.split("+")[0].strip() for n in alignment2.names])

    if "tree" in kwargs:
        tip_names = set(
            [n.split("+")[0].strip() for n in kwargs["tree"].get_tip_names()],
        )
        assert (
            alignment1_names == alignment2_names == tip_names
        ), "Alignment and tree sequence names must perfectly overlap"
    else:
        # no tree passed in
        assert (
            alignment1_names == alignment2_names
        ), "Alignment sequence names must perfectly overlap"

    # Determine if the alignments have enough sequences to proceed.
    if alignment1.num_seqs < min_num_seqs:
        raise ValueError(
            "Too few sequences in merged alignment: %d < %d"
            % (alignment1.num_seqs, min_num_seqs),
        )

    # Confirm that min_num_seqs <= max_num_seqs
    if max_num_seqs and min_num_seqs > max_num_seqs:
        raise ValueError(
            "min_num_seqs (%d) cannot be greater than max_num_seqs (%d)."
            % (min_num_seqs, max_num_seqs),
        )


# End general validation functions

# Start alignment-wide intramolecular coevolution analysis


# coevolve alignment functions: f(alignment,**kwargs) -> 2D array
coevolve_alignment_functions = {
    "mi": mi_alignment,
    "nmi": normalized_mi_alignment,
    "rmi": resampled_mi_alignment,
    "sca": sca_alignment,
    "an": ancestral_state_alignment,
}


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolve_alignment(method, alignment, **kwargs):  # pragma: no cover
    """Apply coevolution method to alignment (for intramolecular coevolution)

    method: f(alignment,**kwargs) -> 2D array of coevolution scores
    alignment: alignment object for which coevolve scores should be
        calculated
    **kwargs: parameters to be passed to method()
    """
    # Perform method specific validation steps
    if method == sca_alignment:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_alignment:
        ancestral_states_input_validation(alignment, **kwargs)
    validate_alignment(alignment)
    return method(alignment, **kwargs)


# End alignment-wide intramolecular coevolution analysis

# Start intermolecular coevolution analysis


# Mapping between coevolve_alignment functions and coevolve_pair functions.
# These are used in coevolve_alignments, b/c under some circumstance the
# alignment function is used, and under other circumstance the pair
# function is used, but the user shouldn't have to know anything about
# that.
coevolve_alignment_to_coevolve_pair = {
    mi_alignment: mi_pair,
    normalized_mi_alignment: normalized_mi_pair,
    resampled_mi_alignment: resampled_mi_pair,
    sca_alignment: sca_pair,
    ancestral_state_alignment: ancestral_state_pair,
}


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def merge_alignments(alignment1, alignment2):  # pragma: no cover
    """Append alignment 2 to the end of alignment 1

    This function is used by coevolve_alignments to merge two alignments
     so they can be evaluated by coevolve_alignment.
    """
    result = {}
    # Created maps from the final seq ids (i.e., seq id before plus) to the
    # seq ids in the original alignments
    aln1_name_map = dict([(n.split("+")[0].strip(), n) for n in alignment1.names])
    aln2_name_map = dict([(n.split("+")[0].strip(), n) for n in alignment2.names])

    try:
        for merged_name, orig_name in list(aln1_name_map.items()):
            result[merged_name] = alignment1.get_gapped_seq(
                orig_name,
            ) + alignment2.get_gapped_seq(aln2_name_map[merged_name])
    except ValueError:  # Differing MolTypes
        for merged_name, orig_name in list(aln1_name_map.items()):
            result[merged_name] = Sequence(
                alignment1.get_gapped_seq(orig_name),
            ) + Sequence(alignment2.get_gapped_seq(aln2_name_map[merged_name]))
    except KeyError:
        raise KeyError(
            "A sequence identifier is in alignment2 "
            + "but not alignment1 -- did you filter out sequences identifiers"
            + " not common to both alignments?",
        )
    return make_aligned_seqs(result, array_align=True)


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def n_random_seqs(alignment, n):  # pragma: no cover
    """Given alignment, return n random seqs in a new alignment object.

    This function is used by coevolve_alignments.

    """
    seq_names = alignment.names
    shuffle(seq_names)
    return alignment.take_seqs(seq_names[:n])


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolve_alignments(
    method,
    alignment1,
    alignment2,
    return_full=False,
    merged_aln_filepath=None,
    min_num_seqs=2,
    max_num_seqs=None,
    sequence_filter=n_random_seqs,
    **kwargs,
):  # pragma: no cover
    """Apply method to a pair of alignments (for intermolecular coevolution)

    method: the *_alignment function to be applied
    alignment1: alignment of first molecule (ArrayAlignment)
    alignment2: alignment of second molecule (ArrayAlignment)
    return_full: if True, returns intra- and inter-molecular
     coevolution data in a square matrix (default: False)
    merged_aln_filepath: if provided, will write the merged
     alignment to file (useful for running post-processing filters)
    min_num_seqs: the minimum number of sequences that should be
     present in the merged alignment to perform the analysis
     (default: 2)
    max_num_seqs: the maximum number of sequences to include
     in an analysis - if the number of sequences exceeds
     max_num_seqs, a random selection of max_num_seqs will be
     used. This is a time-saving step as too many sequences can
     slow things down a lot. (default: None, any number of
     sequences is allowed)
    sequence_filter: function which takes an alignment and an int
     and returns the int number of sequences from the alignment in
     a new alignment object (defualt: util.n_random_seqs(alignment,n))
     if None, a ValueError will be raised if there are more than
     max_num_seqs

    This function allows for calculation of coevolve scores between
     pairs of alignments. The results are returned in a rectangular
     len(alignment1) x len(alignment2) matrix.

    There are some complications involved in preparing alignments for
     this function, because it needs to be obvious how to associate the
     putative interacting sequences. For example, if looking for
     interactions between mammalian proteins A and B, sequences are
     required from the same sets of species, and it must be apparant how
     to match the sequences that are most likely to be involved in
     biologically meaningful interactions. This typically means matching
     the sequences of proteins A&B that come from the same species. In
     other words, interaction of T. aculeatus proteinA and
     H. sapien proteinB likely don't form a biologically relevant
     interaction, because the species are so diverged.

     Matching of sequences is performed via the identifiers, but it is
     the responsibility of the user to correctly construct the sequence
     identifiers before passing the alignments (and tree, if applicable)
     to this function. To faciliate matching sequence identifiers, but not
     having to discard the important information already present in a
     sequence identifier obtained from a database such as KEGG or RefSeq,
     sequence identifiers may contain a plus symbol (+). The characters
     before the + are used to match sequences between the alignments and
     tree. The characters after the + are ignored by this function. So, a
     good strategy is to make the text before the '+' a taxonomic
     identifier and leave the text after the '+' as the original sequence
     identifier. For example, your sequence/tip names could look like:

     alignment1: 'H. sapien+gi|123', 'T. aculeatus+gi|456'
     alignment2: 'T. aculeatus+gi|999', 'H. sapien+gi|424'
     tree: 'T. aculeatus+gi|456', 'H. sapien'

     If there is no plus, the full sequence identifier will be used for the
     matching (see H. sapien in tree).  The order of sequences in the
     alignments is not important. Also note that we can't split on a colon,
     as would be convenient for pulling sequences from KEGG, because colons
     are special characters in newick.

     A WORD OF WARNING ON SEQUENCE IDENTIFIER CONSTRUCTION:
     A further complication is that in some cases, an organism will have
     multiple copies of proteins involved in a complex, but proteinA from
     locus 1 will not form a functional comples with proteinB from locus 2.
     An example of this is the three T6SSs in P. aeuroginosa. Make sure
     this is handled correctly when building your sequence identifiers!
     Sequence identifiers are used to match the sequences which are
     suspected to form a functional complex, which may not simply mean
     sequences from the same species.

    """
    # Perform general validation step
    coevolve_alignments_validation(
        method,
        alignment1,
        alignment2,
        min_num_seqs,
        max_num_seqs,
        **kwargs,
    )
    # Append alignment 2 to the end of alignment 1 in a new alignment object
    merged_alignment = merge_alignments(alignment1, alignment2)
    validate_alignment(merged_alignment)

    if max_num_seqs and merged_alignment.num_seqs > max_num_seqs:
        try:
            merged_alignment = sequence_filter(merged_alignment, max_num_seqs)
        except TypeError:
            raise ValueError("Too many sequences for covariation analysis.")

    # If the user provided a filepath for the merged alignment, write it to
    # disk. This is sometimes useful for post-processing steps.
    if merged_aln_filepath:
        merged_aln_file = open(merged_aln_filepath, "w")
        merged_aln_file.write(merged_alignment.to_fasta())
        merged_aln_file.close()

    if return_full:
        # If the user requests the full result matrix (inter and intra
        # molecular coevolution data), call coevolve_alignment on the
        # merged alignment. Calling coevolve_alignment ensures that
        # the correct validations are performed, rather than directly
        # calling method.
        result = coevolve_alignment(method, merged_alignment, **kwargs)
        return result

    # Note: we only get here if the above if statement comes back False,
    # i.e., if we only want the intermolecular coevolution and don't care
    # about the intramolecular coevolution.

    # Get the appropriate method (need the pair method,
    # not the alignment method)
    try:
        method = coevolve_alignment_to_coevolve_pair[method]
    except KeyError:
        # may have passed in the coevolve_pair function, so just
        # continue -- will fail (loudly) soon enough if not.
        pass

    # Cache the alignment lengths b/c we use them quite a bit, and build
    # the result object to be filled in.
    len_alignment1 = len(alignment1)
    len_alignment2 = len(alignment2)
    result = array([[DEFAULT_NULL_VALUE] * len_alignment1] * len_alignment2)

    # Some of the methods run much faster if relevant data is computed once,
    # and passed in -- that is done here, but there is a lot of repeated code.
    # I'm interested in suggestions for how to make this block of code more
    # compact (e.g., can I be making better use of kwargs?).
    if method == mi_pair or method == nmi_pair or method == normalized_mi_pair:
        positional_entropies = [
            CategoryCounter(p).entropy for p in merged_alignment.positions
        ]
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(
                    merged_alignment,
                    j + len_alignment1,
                    i,
                    h1=positional_entropies[j + len_alignment1],
                    h2=positional_entropies[i],
                    **kwargs,
                )
    elif method == ancestral_state_pair:
        # Perform method-specific validations so we can safely work
        # directly with method rather than the coevolve_pair wrapper,
        # and thereby avoid validation steps on each call to method.
        ancestral_states_input_validation(merged_alignment, **kwargs)
        ancestral_seqs = get_ancestral_seqs(merged_alignment, kwargs["tree"])
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(
                    aln=merged_alignment,
                    pos1=j + len_alignment1,
                    pos2=i,
                    ancestral_seqs=ancestral_seqs,
                    **kwargs,
                )
    else:
        # Perform method-specific validations so we can safely work
        # directly with method rather than the coevolve_pair wrapper,
        # and thereby avoid validation steps on each call to method.
        if method == sca_pair:
            sca_input_validation(merged_alignment, **kwargs)
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(merged_alignment, j + len_alignment1, i, **kwargs)
    return result


# End intermolecular coevolution analysis

# Start positional coevolution analysis
# coevolve position functions: f(alignment,position,**kwargs) -> 1D array
coevolve_position_functions = {
    "mi": mi_position,
    "nmi": normalized_mi_position,
    "rmi": resampled_mi_position,
    "sca": sca_position,
    "an": ancestral_state_position,
}


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolve_position(method, alignment, position, **kwargs):  # pragma: no cover
    """Apply provided coevolution method to a column in alignment

    method: f(alignment,position,**kwargs) -> array of coevolution scores
    alignment: alignment object for which coevolve scores should be
        calculated (ArrayAlignment)
    position: position of interest for coevolution analysis (int)
    **kwargs: parameters to be passed to method()
    """
    # Perform method-specific validation steps
    if method == sca_position:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_position:
        ancestral_states_input_validation(alignment, **kwargs)
    # Perform general validation steps
    validate_position(alignment, position)
    validate_alignment(alignment)
    # Perform the analysis and return the result vector
    return method(alignment, position=position, **kwargs)


# End positional coevolution analysis


# Start pairwise coevolution analysis
# coevolve pair functions: f(alignment,pos1,pos2,**kwargs) -> float
coevolve_pair_functions = {
    "mi": mi_pair,
    "nmi": normalized_mi_pair,
    "rmi": resampled_mi_pair,
    "sca": sca_pair,
    "an": ancestral_state_pair,
}


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolve_pair(method, alignment, pos1, pos2, **kwargs):  # pragma: no cover
    """Apply provided coevolution method to columns pos1 & pos2 of alignment

    method: f(alignment,pos1,pos2,**kwargs) -> coevolution score
    alignment: alignment object for which coevolve score should be
        calculated (ArrayAlignment)
    pos1, pos2: positions to evaluate coevolution between (int)
    **kwargs: parameters to be passed to method()

    """
    # Perform method-specific validation steps
    if method == sca_pair:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_pair:
        ancestral_states_input_validation(alignment, **kwargs)
    # Perform general validation steps
    validate_position(alignment, pos1)
    validate_position(alignment, pos2)
    validate_alignment(alignment)
    # Perform the analysis and return the result score
    return method(alignment, pos1=pos1, pos2=pos2, **kwargs)


# End pairwise coevolution analysis
# End methods for running coevolutionary analyses on sequence data


# Coevolution matrix filters: the following functions are used as
# post-processing filters for coevolution result matrices.


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def filter_threshold_based_multiple_interdependency(
    aln,
    coevolution_matrix,
    threshold=0.95,
    max_cmp_threshold=1,
    cmp_function=greater_equal,
    intermolecular_data_only=False,
):  # pragma: no cover
    """Filters positions with more than max_cmp_threshold scores >= threshold

    This post-processing filter is based on the idea described in:
     "Using multiple interdependency to separate functional from
      phylogenetic correlations in protein alignments"
      Tillier and Lui, 2003

    The idea is that when a position achieved a high covariation score
     with many other positions, the covariation is more likely to arise
     from the phylogeny than from coevolution. They illustrate that this
     works in their paper, and I plan to test it with my alpha-helix-based
     analysis. Note that you can change cmp_function to change whether
     you're looking for high values to indicate covarying positions
     (cmp_function=greater_equal, used for most coevolution algorithms) or
     low values to indicate covarying positions (cmp_function=less_equal,
     used, e.g., for p-value matrices).

    aln: alignment used to generate the coevolution matrix -- this
     isn't actually used, but is required to maintain the same interface
     as other post-processing filters. Pass None if that's more convenient.
    coevolution_matrix: the 2D numpy array to be filtered. This should
     be a rectangular matrix for intermoelcular coevolution data (in which
     case intermolecular_data_only must be set to True) or a symmetric
     square matrix (when intermolecular_data_only=False)
    threshold: the threshold coevolution score that other scores should be
     compared to
    max_cmp_threshold: the max number of scores that are allowed to be
     True with respect to cmp_function and threshold (e.g., the max number
     of positions that may be greater than the threhsold) before setting
     all values associated that position to gDefaultNullValue (default: 1)
    cmp_function: the function that compares each score in
     coevolution_matrix to threshold (default: ge (greater than)) -
     function should return True if the score is one that your looking
     (e.g. score >= threshold) or False otherwise
    intermolecular_data_only: True if coevolution_matrix is a rectangular
     matrix representing an intermolecular coevolution study, and False
     if the matrix is a symmetric square matrix

    NOTE: IF intermolecular_data_only == True, coevolution_matrix MUST BE
     SYMMETRIC, NOT LOWER TRIANGULAR OR OTHERWISE NON-SYMMETRIC!!
    """
    # Determine which rows need to be filtered (but don't filter them
    # right away or subsequent counts could be off)
    filtered_rows = []
    for row_n in range(coevolution_matrix.shape[0]):
        count_cmp_threshold = 0
        for v in coevolution_matrix[row_n, :]:
            if v != DEFAULT_NULL_VALUE and cmp_function(v, threshold):
                count_cmp_threshold += 1
                if count_cmp_threshold > max_cmp_threshold:
                    filtered_rows.append(row_n)
                    break

    # if the matrix is not symmetric, determine which cols need to be filtered
    if intermolecular_data_only:
        filtered_cols = []
        for col_n in range(coevolution_matrix.shape[1]):
            count_cmp_threshold = 0
            for v in coevolution_matrix[:, col_n]:
                if v != DEFAULT_NULL_VALUE and cmp_function(v, threshold):
                    count_cmp_threshold += 1
                    if count_cmp_threshold > max_cmp_threshold:
                        filtered_cols.append(col_n)
                        break
        # filter the rows and cols in a non-symmetric matrix
        for row_n in filtered_rows:
            coevolution_matrix[row_n, :] = DEFAULT_NULL_VALUE
        for col_n in filtered_cols:
            coevolution_matrix[:, col_n] = DEFAULT_NULL_VALUE
    else:
        # filter the rows and cols in a symmetric matrix
        for row_n in filtered_rows:
            coevolution_matrix[row_n, :] = coevolution_matrix[:, row_n] = (
                DEFAULT_NULL_VALUE
            )

    # return the result
    return coevolution_matrix


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def is_parsimony_informative(
    column_freqs,
    minimum_count=2,
    minimum_differences=2,
    ignored=DEFAULT_EXCLUDES,
    strict=False,
):  # pragma: no cover
    """Return True is aln_position is parsimony informative

    column_freqs: dict of characters at alignmnet position mapped
     to their counts -- this is the output of call alignment.column_freqs()
    minimum_count: the minimum number of times a character must show up
     for it to be acceptable (default: 2)
    minimum_differences: the minimum number of different characters
     that must show up at the alignment position (default: 2)
    ignored: characters that should not be counted toward
     minimum_differences (default are exclude characters)
    strict: if True, requires that all amino acids showing up at least
     once at the alignment position show up at least minimum_counts
     times, rather than only requiring that minimum_differences
     amino acids show up minimum_counts times. (default: False)

    The term parsimony informative comes from Codoner, O'Dea,
     and Fares 2008, Reducing the false positive rate in the non-
     parametric analysis of molecular coevolution. In the paper
     they find that if positions which don't contain at least two
     different amino acids, and where each different amino acid doesnt
     show up at least twice each are ignored (i.e., treated as though
     there is not enough information) the positive predictive value
     (PPV) and sensitivity (SN) increase on simulated alignments. They
     term this quality parsimony informative.
     I implemented this as a filter, but include some generalization.
     To determine if a column in an alignment is parsimony informative
     in the exact manner described in Codoner et al., the following
     parameter settings are required:
      minimum_count = 2 (default)
      minimum_differences = 2 (default)
      strict = True (default is False)
     To generalize this function, minimum_count and minimum_differences
     can be passed in so at least minimum_differences different amino
     acids must show up, and each amino acid must show up at least
     minimum_count times.
     In additional variation, strict=False can be passed requiring
     that only minimum_differences number of amino acids show up at least
     minimum_counts times (opposed to requiring that ALL amino acids show
     up minimum_counts times). This is the default behavior.
     By default, the default exclude characters (- and ?) don't count.

    """
    try:
        column_freqs = column_freqs.to_dict()
    except AttributeError:
        pass
    ignored = None if not ignored else list(set(ignored) & set(column_freqs.keys()))
    if ignored:
        for e_ in ignored:
            try:
                del column_freqs[e_]
            except KeyError:
                pass

    if len(column_freqs) < minimum_differences:
        return False
    count_gte_minimum = 0
    for count in list(column_freqs.values()):
        # if not strict, only minimum_differences of the counts
        # must be greater than or equal to minimum_count, so
        # count those occurrences (this is different than the
        # exact technique presented in Codoner et al.)
        if count >= minimum_count:
            count_gte_minimum += 1
        # if strict, all counts must be greater than minimum_count
        # so return False here if we find one that isn't. This is how
        # the algorithm is described in Codoner et al.
        elif strict:
            return False
    return count_gte_minimum >= minimum_differences


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def filter_non_parsimony_informative(
    aln,
    coevolution_matrix,
    null_value=DEFAULT_NULL_VALUE,
    minimum_count=2,
    minimum_differences=2,
    ignored=DEFAULT_EXCLUDES,
    intermolecular_data_only=False,
    strict=False,
):  # pragma: no cover
    """Replaces scores in coevolution_matrix with null_value for positions
    which are not parsimony informative.

    See is_parsimony_informative doc string for definition of
     parsimony informative.

    aln: the input alignment used to generate the coevolution matrix;
     if the alignment was recoded, this should be the recoded alignment.
    coevolution_matrix: the result matrix
    null_value: the value to place in positions which are not
     parsimony informative
    """
    if intermolecular_data_only:
        len_aln1 = coevolution_matrix.shape[1]
    column_frequencies = aln.counts_per_pos()
    for i in range(len(column_frequencies)):
        if not is_parsimony_informative(
            column_frequencies[i],
            minimum_count,
            minimum_differences,
            ignored,
            strict,
        ):
            if not intermolecular_data_only:
                coevolution_matrix[i, :] = coevolution_matrix[:, i] = null_value
            else:
                try:
                    coevolution_matrix[:, i] = null_value
                except IndexError:
                    coevolution_matrix[i - len_aln1, :] = null_value


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def make_positional_exclude_percentage_function(
    excludes,
    max_exclude_percent,
):  # pragma: no cover
    """return function to identify aln positions with > max_exclude_percent"""
    excludes = {}.fromkeys(excludes)

    def f(col):
        exclude_count = 0
        for c in col:
            if c in excludes:
                exclude_count += 1
        return exclude_count / len(col) > max_exclude_percent

    return f


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def filter_exclude_positions(
    aln,
    coevolution_matrix,
    max_exclude_percent=0.1,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    intermolecular_data_only=False,
):  # pragma: no cover
    """Assign null_value to positions with > max_exclude_percent excludes

    aln: the ArrayAlignment object
    coevolution_matrix: the 2D numpy array -- this will be modified
    max_exclude_percent: the maximimu percent of characters that
     may be exclude characters in any alignment position (column).
     if the percent of exclude characters is greater than this value,
     values in this position will be replaced with null_value
     (default = 0.10)
    null_value: the value to be used as null (default: gDefaultNullValue)
    excludes: the exclude characters (default: gDefaultExcludes)
    intermolecular_data_only: True if the coevolution result
     matrix contains only intermolecular data (default: False)

    """
    # construct the function to be passed to aln.get_position_indices
    f = make_positional_exclude_percentage_function(excludes, max_exclude_percent)
    # identify the positions containing too many exclude characters
    exclude_positions = aln.get_position_indices(f)

    # replace values from exclude_positions with null_value
    if not intermolecular_data_only:
        # if working with intramolecular data (or inter + intra molecular data)
        # this is easy
        for p in exclude_positions:
            coevolution_matrix[p, :] = coevolution_matrix[:, p] = null_value
    else:
        # if working with intermolecular data only, this is more complicated --
        # must convert from alignment positions to matrix positions
        len_aln1 = coevolution_matrix.shape[1]
        for p in exclude_positions:
            try:
                coevolution_matrix[:, p] = null_value
            except IndexError:
                coevolution_matrix[p - len_aln1, :] = null_value


# Functions for archiving/retrieiving coevolve results
# These functions are extremely general -- should they go
# somewhere else, or should I be using pre-existing code?


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def pickle_coevolution_result(
    coevolve_result,
    out_filepath="output.pkl",
):  # pragma: no cover
    """Pickle coevolve_result and store it at output_filepath

    coevolve_result: result from a coevolve_* function (above); this can
     be a float, an array, or a 2D array (most likely it will be one of the
     latter two, as it will usually be fast enough to compute a single
     coevolve value on-the-fly.
    out_filepath: path where the pickled result should be stored
    """
    try:
        infile = open(out_filepath, "wb")
        p = Pickler(infile)
    except OSError:
        err = "Can't access filepath. Do you have write access? " + out_filepath
        raise OSError(err)
    p.dump(coevolve_result)
    infile.close()


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def unpickle_coevolution_result(in_filepath):  # pragma: no cover
    """Read in coevolve_result from a pickled file

    in_filepath: filepath to unpickle
    """
    try:
        infile = open(in_filepath, "rb")
        u = Unpickler(infile)
    except OSError:
        err = (
            "Can't access filepath. Does it exist? Do you have read access? "
            + in_filepath
        )
        raise OSError(err)
    r = u.load()
    infile.close()
    return r


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def coevolution_matrix_to_csv(
    coevolve_matrix,
    out_filepath="output.csv",
):  # pragma: no cover
    """Write coevolve_matrix as csv file at output_filepath

    coevolve_result: result from a coevolve_alignment function (above);
     this should be a 2D numpy array
    out_filepath: path where the csv result should be stored
    """
    try:
        f = open(out_filepath, "w")
    except OSError:
        err = "Can't access filepath. Do you have write access? " + out_filepath
        raise OSError(err)
    f.write("\n".join([",".join([str(v) for v in row]) for row in coevolve_matrix]))
    f.close()


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def csv_to_coevolution_matrix(in_filepath):  # pragma: no cover
    """Read a coevolution matrix from a csv file

    in_filepath: input filepath
    """
    try:
        f = open(in_filepath)
    except OSError:
        err = (
            "Can't access filepath. Does it exist? Do you have read access? "
            + in_filepath
        )
        raise OSError(err)
    result = []
    for line in f:
        values = line.strip().split(",")
        result.append(list(map(float, values)))
    f.close()
    return array(result)


# End functions for archiving/retrieiving coevolve results

# Start functions for analyzing the results of a coevolution run.


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def identify_aln_positions_above_threshold(
    coevolution_matrix,
    threshold,
    aln_position,
    null_value=DEFAULT_NULL_VALUE,
):  # pragma: no cover
    """Returns the list of alignment positions which achieve a
    score >= threshold with aln_position.
    Coevolution matrix should be symmetrical or you
    may get weird results -- scores are pulled from the row describing
    aln_position.
    """
    coevolution_scores = coevolution_matrix[aln_position]
    results = []
    for i in range(len(coevolution_scores)):
        s = coevolution_scores[i]
        if s != null_value and s >= threshold:
            results.append(i)
    return results


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def aln_position_pairs_cmp_threshold(
    coevolution_matrix,
    threshold,
    cmp_function,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):  # pragma: no cover
    """Returns list of position pairs with score >= threshold

    coevolution_matrix: 2D numpy array
    threshold: value to compare matrix positions against
    cmp_function: function which takes a value and theshold
     and returns a boolean (e.g., ge(), le())
    null_value: value representing null scores -- these are
     ignored
    intermolecular_data_only: True if the coevolution result
     matrix contains only intermolecular data (default: False)
    """
    if not intermolecular_data_only:
        assert (
            coevolution_matrix.shape[0] == coevolution_matrix.shape[1]
        ), "Non-square matrices only supported for intermolecular-only data."
    results = []
    # compile the matrix positions with cmp(value,threshold) == True
    for i, row in enumerate(coevolution_matrix):
        for j, value in enumerate(row):
            if value != null_value and cmp_function(value, threshold):
                results.append((i, j))

    # if working with intermolecular data only, need to convert
    # matrix positions to alignment positions
    if intermolecular_data_only:
        # convert matrix positions to alignment positions
        adjustment = coevolution_matrix.shape[1]
        results = [(j, i + adjustment) for i, j in results]
    return results


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def aln_position_pairs_ge_threshold(
    coevolution_matrix,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):  # pragma: no cover
    """wrapper function for aln_position_pairs_cmp_threshold"""
    return aln_position_pairs_cmp_threshold(
        coevolution_matrix,
        threshold,
        greater_equal,
        null_value,
        intermolecular_data_only,
    )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def aln_position_pairs_le_threshold(
    coevolution_matrix,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):  # pragma: no cover
    """wrapper function for aln_position_pairs_cmp_threshold"""
    return aln_position_pairs_cmp_threshold(
        coevolution_matrix,
        threshold,
        less_equal,
        null_value,
        intermolecular_data_only,
    )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def count_cmp_threshold(
    m,
    threshold,
    cmp_function,
    null_value=DEFAULT_NULL_VALUE,
    symmetric=False,
    ignore_diagonal=False,
):  # pragma: no cover
    """Returns a count of the values in m >= threshold, ignoring nulls.

    m: coevolution matrix (numpy array)
    thresold: value to compare against scores in matrix (float)
    cmp_function: function used to compare value to threshold
     (e.g., greater_equal, less_equal)
    """

    total_non_null = 0
    total_hits = 0
    if not symmetric:
        if ignore_diagonal:
            values = [
                m[i, j] for i in range(m.shape[0]) for j in range(m.shape[1]) if i != j
            ]
        else:
            values = m.flat
    elif ignore_diagonal:
        # has to be a better way to do this... tril doesn't work b/c it
        # sets the upper triangle to zero -- if i could get it to set
        # that to null_value, and then apply flat, that'd be fine.
        # values = tril(m,-1)
        values = [m[i, j] for i in range(len(m)) for j in range(i)]
    else:
        # values = tril(m)
        values = [m[i, j] for i in range(len(m)) for j in range(i + 1)]

    if isnan(null_value):

        def is_not_null_value(v):
            return not isnan(v)

    else:

        def is_not_null_value(v):
            return isnan(v) or v != null_value

    for value in values:
        if is_not_null_value(value):
            total_non_null += 1
            if cmp_function(value, threshold):
                total_hits += 1
    return total_hits, total_non_null


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def count_ge_threshold(
    m,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    symmetric=False,
    ignore_diagonal=False,
):  # pragma: no cover
    """wrapper function for count_cmp_threshold"""
    return count_cmp_threshold(
        m,
        threshold,
        greater_equal,
        null_value,
        symmetric,
        ignore_diagonal,
    )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def count_le_threshold(
    m,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    symmetric=False,
    ignore_diagonal=False,
):  # pragma: no cover
    """wrapper function for count_cmp_threshold"""
    return count_cmp_threshold(
        m,
        threshold,
        less_equal,
        null_value,
        symmetric,
        ignore_diagonal,
    )


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def ltm_to_symmetric(m):  # pragma: no cover
    """Copies values from lower triangle to upper triangle"""
    assert (
        m.shape[0] == m.shape[1]
    ), "Making matrices symmetric only supported for square matrices"

    for i in range(len(m)):
        for j in range(i):
            m[j, i] = m[i, j]
    return m


# End functions for analyzing the results of a coevolution run


# Script functionality
@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def build_coevolution_matrix_filepath(
    input_filepath,
    output_dir="./",
    method=None,
    alphabet=None,
    parameter=None,
):  # pragma: no cover
    """ Build filepath from input filename, output dir, and list of suffixes

        input_filepath: filepath to be used for generating the output
            filepath. The path and the final suffix will be stripped to
            get the 'base' filename.
        output_dir: the path to append to the beginning of the base filename
        method: string indicating method that should be appended to filename
        alphabet: string indicating an alphabet recoding which should be
            appended to filename, or None
        parameter: parameter that should be appended to the filename,
            or None (ignored if method doesn't require parameter)

        Examples:
         >>> build_coevolution_matrix_filepath(\
          './p53.fasta','/output/path','mi','charge')
         /output/path/p53.charge.mi
         >>> build_coevolution_matrix_filepath(\
          './p53.new.fasta','/output/path','mi','charge')
         /output/path/p53.new.charge.mi
         >>> build_coevolution_matrix_filepath(\
          './p53.fasta','/output/path','sca','charge',0.75)
         /output/path/p53.charge.sca_75

    """
    if method == "sca":
        try:
            cutoff_str = str(parameter)
            point_index = cutoff_str.rindex(".")
            method = "_".join([method, cutoff_str[point_index + 1 : point_index + 4]])
        except ValueError:
            raise ValueError("Cutoff must be provided when method == 'sca'")

    suffixes = [_f for _f in [alphabet, method] if _f]

    # strip path
    try:
        result = input_filepath[input_filepath.rindex("/") + 1 :]
    except ValueError:
        result = input_filepath
    # strip final suffix
    try:
        result = result[: result.rindex(".")]
    except ValueError:
        pass
    # append output path
    if output_dir.endswith("/"):
        result = "".join([output_dir, result])
    else:
        result = "".join([output_dir, "/", result])
    # append output suffixes
    result = ".".join([_f for _f in [result] + suffixes if _f])

    return result


@c3warn.deprecated_callable("2024.12", reason=_reason_3, is_discontinued=True)
def parse_coevolution_matrix_filepath(filepath):  # pragma: no cover
    """Parses a coevolution matrix filepath into constituent parts.

    Format is very specific. Will only work on filenames such as:
     path/alignment_identifier.alphabet_id.method.pkl
     path/alignment_identifier.alphabet_id.method.csv

    This format is the recommended naming convention for coevolution
     matrices. To ensure filepaths compatible with this function, use
     cogent3.evolve.coevolution.build_coevolution_matrix_filepath to build
     the filepaths for your coevolution matrices.


     Examples:
     parse_coevolution_matrix_filepath('pkls/myosin_995.a1_4.nmi.pkl')
        => ('myosin_995', 'a1_4', 'nmi')
     parse_coevolution_matrix_filepath('p53.orig.mi.csv')
        => ('p53','orig','mi')
    """
    filename = basename(filepath)
    fields = filename.split(".")
    try:
        alignment_id = fields[0]
        alphabet_id = fields[1]
        method_id = fields[2]
        _ = fields[3]  # extension
    except IndexError:
        raise ValueError(
            "output filepath not in parsable format: %s. See doc string for format definition."
            % filepath,
        )

    return (alignment_id, alphabet_id, method_id)

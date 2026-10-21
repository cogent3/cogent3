import itertools
from collections.abc import Callable
from typing import TYPE_CHECKING, Literal

import numba
import numba.types as numba_types
import numpy

from cogent3.core import moltype as c3_moltype
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    get_purine_indices,
    get_pyrimidine_indices,
)

if TYPE_CHECKING:  # pragma: no cover
    # Forward reference to avoid circular import
    from cogent3.core.alignment import Alignment

# turn off code coverage as jit-ted code not accessible to coverage


@numba.jit(cache=True)
def index_to_lower_tri(index):  # pragma: no cover
    """Convert linear index to lower-triangular (i, j) coordinates.

    Maps a flat index over unique pairs 0..n*(n-1)/2-1 to the
    corresponding (row, col) in the lower triangle of an n x n matrix.
    """
    i = int((1 + (1 + 8 * index) ** 0.5) / 2)
    j = index - i * (i - 1) // 2
    return i, j


# fills in a diversity matrix from sequences of integers


@numba.jit
def num_diffs_and_valid(
    seq1: numba_types.uint8[:],
    seq2: numba_types.uint8[:],
    num_states: numba_types.uint8,
) -> tuple[int, int]:  # pragma: no cover
    """returns the number of differences and the number of valid positions"""
    num_valid = 0
    num_diff = 0
    for i in range(len(seq1)):
        if seq1[i] >= num_states or seq2[i] >= num_states:
            continue
        num_valid += 1
        if seq1[i] != seq2[i]:
            num_diff += 1

    return num_diff, num_valid


@numba.jit(cache=True)
def _get_matrix_and_counts(matrix, seq1, seq2):  # pragma: no cover
    """fills the diversity matrix for valid positions.

    Parameters
    ----------
    matrix
        Counts of states
    seq1, seq2
        Sequences to count variations between
    """
    num_valid = 0
    num_diff = 0
    num_states = matrix.shape[0]
    for i in range(len(seq1)):
        if seq1[i] >= num_states or seq2[i] >= num_states:
            continue
        num_valid += 1
        if seq1[i] != seq2[i]:
            num_diff += 1
        matrix[seq1[i], seq2[i]] += 1
    return matrix, num_diff, num_valid


# statistics computed by pairwise_count_matrix
_HAMMING = 0
_PDIST = 1
_JC69 = 2
SHARED_SITES = 3
MATCHING_SITES = 4


@numba.jit
def _count_statistic(method, num_diffs, num_valid, num_states):  # pragma: no cover
    if method == _HAMMING:
        value = float(num_diffs)
    elif method == SHARED_SITES:
        value = float(num_valid)
    elif method == MATCHING_SITES:
        value = float(num_valid - num_diffs)
    elif num_valid == 0:
        value = numpy.nan
    elif method == _PDIST or num_diffs == 0:
        # also gives jc69 a distance of 0.0, not -0.0, for identical pairs
        value = num_diffs / num_valid
    elif num_diffs * num_states >= num_valid * (num_states - 1):
        # integer test, as p * frac rounds below 1 at saturation for some states
        value = numpy.nan
    else:
        p = num_diffs / num_valid
        frac = num_states / (num_states - 1)
        value = -numpy.log(1.0 - frac * p) / frac
    return value


@numba.jit(parallel=True)
def pairwise_count_matrix(
    array_seqs, num_states, method, parallel=True
):  # pragma: no cover
    """returns a matrix of a statistic computed from per pair site counts

    Parameters
    ----------
    array_seqs
        2D array of sequences
    num_states
        Number of states, e.g. 4 for nucleotides. A site is valid for a pair
        when both sequences have a state below this value.
    method
        SHARED_SITES gives the number of valid sites, MATCHING_SITES the
        number of valid sites that are identical. The other codes give the
        hamming, pdist and jc69 distances.
    parallel
        Runs in parallel if True, otherwise reduces to single thread

    Returns
    -------
    numpy.ndarray
        Symmetric float32 matrix with a zero diagonal. pdist and jc69 are
        nan for pairs with no valid sites, jc69 also when the proportion of
        differences is >= (num_states - 1) / num_states.
    """
    if method < _HAMMING or method > MATCHING_SITES:
        msg = "unknown method code"
        raise ValueError(msg)

    num_threads = numba.get_num_threads()
    if not parallel:
        # numba still gets all the threads it can, but it only
        # uses the number specified by set_num_threads
        # we restore the original number of threads at the end
        numba.set_num_threads(1)

    n_seqs = array_seqs.shape[0]
    result = numpy.zeros((n_seqs, n_seqs), dtype=numpy.float32)
    num_pairs = n_seqs * (n_seqs - 1) // 2
    # loop is parallelised
    for index in numba.prange(num_pairs):
        i, j = index_to_lower_tri(index)
        num_diffs, num_valid = num_diffs_and_valid(
            array_seqs[i],
            array_seqs[j],
            num_states,
        )
        value = _count_statistic(method, num_diffs, num_valid, num_states)
        result[i, j] = result[j, i] = value

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return result


def _count_matrix(aln: "Alignment", method: int, parallel: bool) -> numpy.ndarray:
    return pairwise_count_matrix(
        aln.array_seqs, len(aln.moltype.alphabet), method, parallel=parallel
    )


def _to_dist_matrix(
    mat: numpy.ndarray, aln: "Alignment", invalid_raises: bool
) -> DistanceMatrix:
    if invalid_raises and numpy.isnan(mat).any():
        msg = "nan's in matrix"
        raise ArithmeticError(msg)
    return DistanceMatrix.from_array_names(mat, aln.names)


def _check_nucleic(aln: "Alignment", calc: str) -> None:
    if not aln.moltype.is_nucleic:
        msg = (
            f"{calc} distance only works with nucleotide alignments, not {aln.moltype}"
        )
        raise c3_moltype.MolTypeError(msg)


def jc69(
    aln: "Alignment",
    invalid_raises: bool = False,
    parallel: bool = False,
) -> DistanceMatrix:
    """returns JC69 pairwise distances for an alignment

    A distance is numpy.nan when a pair shares no valid sites, or when the
    proportion of differences is at least (n - 1) / n for an alphabet of n
    states, e.g. 0.75 for DNA.

    Parameters
    ----------
    aln
        Alignment instance
    invalid_raises
        If True and nan in result, raises an ArithmeticError
    parallel
        If True, uses parallel processing via numba.
    """
    mat = _count_matrix(aln, _JC69, parallel)
    return _to_dist_matrix(mat, aln, invalid_raises)


@numba.jit
def _calc_tn93_dist(
    matrix,
    pur_coords,
    pyr_coords,
    tv_coords,
    freq_purs,
    freq_pyrs,
    coeff1,
    coeff2,
    coeff3,
    total,
):  # pragma: no cover
    pur_ts_diffs = matrix.take(pur_coords).sum() / total
    pyr_ts_diffs = matrix.take(pyr_coords).sum() / total
    tv_diffs = matrix.flatten().take(tv_coords).sum() / total
    term1 = 1 - pur_ts_diffs / coeff1 - tv_diffs / (2 * freq_purs)
    term2 = 1 - pyr_ts_diffs / coeff2 - tv_diffs / (2 * freq_pyrs)
    term3 = 1 - tv_diffs / (2 * freq_purs * freq_pyrs)

    if term1 <= 0 or term2 <= 0 or term3 <= 0:  # log will fail
        return numpy.float32(numpy.nan)

    return (
        -coeff1 * numpy.log(term1)
        - coeff2 * numpy.log(term2)
        - coeff3 * numpy.log(term3)
    )


@numba.jit(parallel=True)
def tn93_dist_matrix(
    array_seqs,
    num_states,
    pur_freqs,
    pyr_freqs,
    pur_coords,
    pyr_coords,
    tv_coords,
    coeff1,
    coeff2,
    coeff3,
    parallel=True,
):  # pragma: no cover
    num_threads = numba.get_num_threads()
    if not parallel:
        numba.set_num_threads(1)
    n_seqs = array_seqs.shape[0]
    dists = numpy.zeros((n_seqs, n_seqs), dtype=numpy.float32)
    nan = numpy.float32(numpy.nan)
    zero = numpy.float32(0.0)

    # making working matrices for each thread
    matrices = numpy.zeros(
        (numba.get_num_threads(), num_states, num_states),
        dtype=numpy.int32,
    )
    num_pairs = n_seqs * (n_seqs - 1) // 2
    # loop is parallelised
    for index in numba.prange(num_pairs):
        i, j = index_to_lower_tri(index)
        # get a working matrix for this thread
        div_matrix = matrices[numba.get_thread_id()]
        div_matrix.fill(0)
        div_matrix, num_diffs, num_valid = _get_matrix_and_counts(
            div_matrix,
            array_seqs[i],
            array_seqs[j],
        )
        if num_valid == 0:
            d = nan
        elif num_diffs == 0:
            d = zero
        else:
            d = _calc_tn93_dist(
                div_matrix.flatten(),
                pur_coords,
                pyr_coords,
                tv_coords,
                pur_freqs,
                pyr_freqs,
                coeff1,
                coeff2,
                coeff3,
                num_valid,
            )
        dists[i, j] = dists[j, i] = d

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return dists


@numba.jit(parallel=True)
def _count_states(array_seqs, num_states, parallel=True):
    # temp function to remove when the new Alignment.get_motif_probs()
    # is refactored to use numba
    num_threads = numba.get_num_threads()
    if not parallel:
        numba.set_num_threads(1)
    n_seqs = array_seqs.shape[0]
    # making working matrices for each thread
    state_counts = numpy.zeros((numba.get_num_threads(), num_states), dtype=numpy.int32)
    for i in numba.prange(n_seqs):
        states = state_counts[numba.get_thread_id()]
        for j in range(len(array_seqs[i])):
            state = array_seqs[i, j]
            if state < num_states:
                states[state] += 1

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)

    return state_counts.sum(axis=0)


def _get_symmetric_within(indices: numpy.ndarray) -> numpy.ndarray:
    dims = 4, 4
    coords = numpy.array(
        [(i, j) for i, j in itertools.product(indices, indices) if i != j],
    )
    return numpy.ravel_multi_index(coords.T, dims=dims)


def _get_symmetric_between(
    indices1: numpy.ndarray,
    indices2: numpy.ndarray,
) -> numpy.ndarray:
    # coordinates
    dims = 4, 4
    coords = numpy.array(
        list(
            itertools.chain(
                itertools.product(indices1, indices2),
                itertools.product(indices2, indices1),
            ),
        ),
    )
    return numpy.ravel_multi_index(coords.T, dims=dims)


def tn93(
    aln: "Alignment",
    invalid_raises: bool = False,
    parallel: bool = False,
) -> DistanceMatrix:
    """returns TN93 pairwise distances for an alignment

    Parameters
    ----------
    aln
        Alignment instance
    invalid_raises
        If True and nan in result, raises an ArithmeticError
    parallel
        If True, uses parallel processing via numba.
    """
    _check_nucleic(aln, "tn93")
    alpha = aln.moltype.alphabet

    counts = _count_states(aln.array_seqs, len(alpha), parallel=parallel)
    freqs = (counts / counts.sum()).astype(numpy.float32)

    pur_indices = numpy.array(get_purine_indices(aln.moltype))
    pyr_indices = numpy.array(get_pyrimidine_indices(aln.moltype))

    # these constants are used across all the pairwise distance calculations
    pur_freqs = freqs.take(pur_indices).sum()
    pur_prods = freqs.take(pur_indices).prod()
    pyr_freqs = freqs.take(pyr_indices).sum()
    pyr_prods = freqs.take(pyr_indices).prod()
    coeff1 = 2 * pur_prods / pur_freqs
    coeff2 = 2 * pyr_prods / pyr_freqs
    coeff3 = 2 * (
        pur_freqs * pyr_freqs
        - (pur_prods * pyr_freqs / pur_freqs)
        - (pyr_prods * pur_freqs / pyr_freqs)
    )

    # coords are for the flattened matrix
    pur_coords = _get_symmetric_within(pur_indices)
    pyr_coords = _get_symmetric_within(pyr_indices)
    tv_coords = _get_symmetric_between(pur_indices, pyr_indices)

    mat = tn93_dist_matrix(
        aln.array_seqs,
        4,
        pur_freqs,
        pyr_freqs,
        pur_coords,
        pyr_coords,
        tv_coords,
        coeff1,
        coeff2,
        coeff3,
        parallel=parallel,
    )
    return _to_dist_matrix(mat, aln, invalid_raises)


# method codes for the logdet family of distances
_PARALINEAR = 0
_LOGDET_TK = 1
_LOGDET = 2


@numba.jit
def _logdet_family(matrix, num_states, method):  # pragma: no cover
    # we replace the missing diagonal states with a
    # pseudocount of 0.5 then normalise
    frequency = matrix.astype(numpy.float32)
    for i in range(num_states):
        if frequency[i, i] == 0:
            frequency[i, i] = 0.5

    frequency /= frequency.sum()
    det = numpy.linalg.det(frequency)
    if det <= 0:
        return numpy.float32(numpy.nan)

    if method == _LOGDET:
        return numpy.float32(-numpy.log(det) / num_states - numpy.log(num_states))

    freqs_0 = frequency.sum(axis=0)
    freqs_1 = frequency.sum(axis=1)
    log_ratio = numpy.log(det / numpy.sqrt((freqs_0 * freqs_1).prod()))
    if method == _PARALINEAR:
        return numpy.float32(-log_ratio / num_states)

    # Tamura and Kumar 2002 adjustment for unequal state frequencies
    coeff = (((freqs_0 + freqs_1) ** 2).sum() / 4 - 1) / (num_states - 1)
    return numpy.float32(coeff * log_ratio)


@numba.jit(parallel=True)
def logdet_family_distance_matrix(
    array_seqs,
    num_states,
    method,
    parallel=True,
):  # pragma: no cover
    num_threads = numba.get_num_threads()
    if not parallel:
        numba.set_num_threads(1)
    n_seqs = array_seqs.shape[0]
    dists = numpy.zeros((n_seqs, n_seqs), dtype=numpy.float32)
    nan = numpy.float32(numpy.nan)
    zero = numpy.float32(0.0)

    # making working matrices for each thread
    matrices = numpy.zeros(
        (numba.get_num_threads(), num_states, num_states),
        dtype=numpy.int32,
    )
    num_pairs = n_seqs * (n_seqs - 1) // 2
    # loop is parallelised
    for index in numba.prange(num_pairs):
        i, j = index_to_lower_tri(index)
        # get a working matrix for this thread
        div_matrix = matrices[numba.get_thread_id()]
        div_matrix.fill(0)
        div_matrix, num_diffs, num_valid = _get_matrix_and_counts(
            div_matrix,
            array_seqs[i],
            array_seqs[j],
        )
        if num_valid == 0:
            d = nan
        elif num_diffs == 0:
            d = zero
        else:
            d = _logdet_family(div_matrix, num_states, method)
        dists[i, j] = dists[j, i] = d

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return dists


def paralinear(
    aln: "Alignment",
    invalid_raises: bool = False,
    parallel: bool = False,
) -> DistanceMatrix:
    """returns matrix of pairwise paralinear distances

    Parameters
    ----------
    aln
        Alignment instance
    invalid_raises
        If True and nan in result, raises an ArithmeticError
    parallel
        If True, uses parallel processing via numba.

    Notes
    -----
    This is limited to 4-state alphabets for now.
    """
    _check_nucleic(aln, "paralinear")
    mat = logdet_family_distance_matrix(
        aln.array_seqs, len(aln.moltype.alphabet), _PARALINEAR, parallel=parallel
    )
    return _to_dist_matrix(mat, aln, invalid_raises)


def logdet(
    aln: "Alignment",
    invalid_raises: bool = False,
    parallel: bool = False,
    use_tk_adjustment: bool = True,
) -> DistanceMatrix:
    """returns matrix of pairwise LogDet distances

    Parameters
    ----------
    aln
        Alignment instance
    invalid_raises
        If True and nan in result, raises an ArithmeticError
    parallel
        If True, uses parallel processing via numba.
    use_tk_adjustment
        If True, applies the Tamura and Kumar (2002) correction for
        unequal state frequencies.

    Notes
    -----
    This is limited to 4-state alphabets for now.
    """
    _check_nucleic(aln, "logdet")
    method = _LOGDET_TK if use_tk_adjustment else _LOGDET
    mat = logdet_family_distance_matrix(
        aln.array_seqs, len(aln.moltype.alphabet), method, parallel=parallel
    )
    return _to_dist_matrix(mat, aln, invalid_raises)


def hamming(
    aln: "Alignment", parallel: bool = False, invalid_raises: bool = False
) -> DistanceMatrix:
    """returns matrix of hamming distances

    Parameters
    ----------
    aln
        Alignment instance
    parallel
        If True, uses parallel processing via numba.
    invalid_raises
        Ignored, accepted for a common calculator signature. A pair with
        no sites where both sequences have a canonical state is 0.
    """
    mat = _count_matrix(aln, _HAMMING, parallel)
    return DistanceMatrix.from_array_names(mat, aln.names)


def pdist(
    aln: "Alignment", parallel: bool = False, invalid_raises: bool = False
) -> DistanceMatrix:
    """returns matrix of pairwise proportion difference distances

    Parameters
    ----------
    aln
        Alignment instance
    parallel
        If True, uses parallel processing via numba.
    invalid_raises
        If True and nan in result, raises an ArithmeticError. A pair with
        no sites where both sequences have a canonical state is nan.
    """
    mat = _count_matrix(aln, _PDIST, parallel)
    return _to_dist_matrix(mat, aln, invalid_raises)


_calculators: dict[
    Literal["paralinear", "logdet", "jc69", "tn93", "hamming", "pdist"],
    Callable[..., DistanceMatrix],
] = {
    "paralinear": paralinear,
    "logdet": logdet,
    "jc69": jc69,
    "tn93": tn93,
    "hamming": hamming,
    "pdist": pdist,
}


def get_distance_calculator(name: str) -> Callable[..., DistanceMatrix]:
    """returns a pairwise distance calculator

    name is converted to lower case"""
    name = name.lower()
    if name not in _calculators:
        msg = f'Unknown pairwise distance calculator "{name}"'
        raise ValueError(msg)

    return _calculators[name]

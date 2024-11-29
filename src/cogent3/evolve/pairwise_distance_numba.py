import itertools

import numba
import numba.types as numba_types
import numpy

from cogent3.core import new_moltype
from cogent3.evolve.fast_distance import DistanceMatrix

# turn off code coverage as jit-ted code not accessible to coverage


# fills in a diversity matrix from sequences of integers


def _is_nucleic(moltype):
    # todo delete when new_type alignments and moltypes are the norm
    try:
        _ = moltype.complement("A")
    except (ValueError, new_moltype.MolTypeError):
        return False
    return True


@numba.jit
def num_diffs_and_valid(
    seq1: numba_types.uint8[:],
    seq2: numba_types.uint8[:],
    num_states: numba_types.uint8,
) -> numba_types.int32[:]:  # pragma: no cover
    """returns the number of differences and the number of valid positions"""
    num_valid = 0
    num_diff = 0
    for i in range(len(seq1)):
        if seq1[i] >= num_states or seq2[i] >= num_states:
            continue
        num_valid += 1
        if seq1[i] != seq2[i]:
            num_diff += 1

    return numpy.array([num_diff, num_valid], dtype=numpy.int32)


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


@numba.jit(parallel=True)
def jc69_dist_matrix(array_seqs, num_states, parallel=True):  # pragma: no cover
    """Compute JC69 distances between all sequence pairs.

    Parameters
    ----------
    array_seqs
        2D array of sequences
    num_states
        Number of states, e.g. 4 for nucleotides
    parallel
        Runs in parallel if True, otherwise reduces to single thread

    Returns
    -------
    array
        Matrix of pairwise JC69 distances

    Notes
    -----
    If the proportion of differences is >= 0.75, the distance is set
    to numpy.nan.
    """
    num_threads = numba.get_num_threads()
    if not parallel:
        # numba still gets all the threads it can, but it only
        # uses the number specified by set_num_threads
        # we restore the original number of threads at the end
        numba.set_num_threads(1)

    n_seqs = array_seqs.shape[0]
    dists = numpy.zeros((n_seqs, n_seqs), dtype=numpy.float32)
    nan = numpy.float32(numpy.nan)
    zero = numpy.float32(0.0)
    frac = numpy.float32(num_states / (num_states - 1))
    dist_scale = numpy.float32((num_states - 1) / -num_states)

    # outer loop is parallelised
    for i in numba.prange(n_seqs - 1):
        for j in range(i + 1, n_seqs):
            num_diffs, num_valid = num_diffs_and_valid(
                array_seqs[i], array_seqs[j], num_states
            )

            if num_valid == 0:
                d = nan
            elif num_diffs == 0:
                d = zero
            else:
                p = num_diffs / num_valid
                d = nan if p >= 0.75 else dist_scale * numpy.log(1.0 - frac * p)
            dists[i, j] = dists[j, i] = d

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return dists


def jc69(
    aln: "Alignment", invalid_raises: bool = False, parallel: bool = False
) -> DistanceMatrix:
    """returns JC69 pairwise distances for an alignment

    Parameters
    ----------
    aln
        Alignment instance
    invalid_raises
        If True and nan in result, raises an ArithmeticError
    parallel
        If True, uses parallel processing via numba.
    """

    num_states = len(aln.moltype.alphabet)
    mat = jc69_dist_matrix(aln.array_seqs, num_states, parallel=parallel)
    if invalid_raises and numpy.isnan(mat).any():
        raise ArithmeticError("nan's in matrix")
    return DistanceMatrix.from_array_names(mat, aln.names)


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
        (numba.get_num_threads(), num_states, num_states), dtype=numpy.int32
    )
    # outer parallel loop
    for i in numba.prange(n_seqs - 1):
        # get a working matrix for this thread
        div_matrix = matrices[numba.get_thread_id()]
        for j in range(i + 1, n_seqs):
            div_matrix.fill(0)
            div_matrix, num_diffs, num_valid = _get_matrix_and_counts(
                div_matrix, array_seqs[i], array_seqs[j]
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
def _coun_states(array_seqs, num_states, parallel=True):
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
        [(i, j) for i, j in itertools.product(indices, indices) if i != j]
    )
    return numpy.ravel_multi_index(coords.T, dims=dims)


def _get_symmetric_between(
    indices1: numpy.ndarray, indices2: numpy.ndarray
) -> numpy.ndarray:
    # coordinates
    dims = 4, 4
    coords = numpy.array(
        list(
            itertools.chain(
                itertools.product(indices1, indices2),
                itertools.product(indices2, indices1),
            )
        )
    )
    return numpy.ravel_multi_index(coords.T, dims=dims)


def tn93(
    aln: "Alignment", invalid_raises: bool = False, parallel: bool = False
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
    if not _is_nucleic(aln.moltype):
        raise new_moltype.MolTypeError(
            f"tn93 distance only works with nucleotide alignments, not {aln.moltype}"
        )

    alpha = aln.moltype.alphabet

    counts = _coun_states(aln.array_seqs, len(alpha), parallel=parallel)
    freqs = (counts / counts.sum()).astype(numpy.float32)

    pur_indices = alpha.to_indices("AG")
    pyr_indices = alpha.to_indices("CT")

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
    if invalid_raises and numpy.isnan(mat).any():
        raise ArithmeticError("nan's in matrix")
    return DistanceMatrix.from_array_names(mat, aln.names)


@numba.jit
def _paralinear(matrix, num_states):
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

    # the inverse matrix of frequency, every element is squared
    denom = numpy.sqrt((frequency.sum(axis=0) * frequency.sum(axis=1)).prod())

    return -numpy.log(det / denom) / num_states


@numba.jit(parallel=True)
def paralinear_distance_matrix(
    array_seqs,
    num_states,
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
        (numba.get_num_threads(), num_states, num_states), dtype=numpy.int32
    )
    # outer parallel loop
    for i in numba.prange(n_seqs - 1):
        # get a working matrix for this thread
        div_matrix = matrices[numba.get_thread_id()]
        for j in range(i + 1, n_seqs):
            div_matrix.fill(0)
            div_matrix, num_diffs, num_valid = _get_matrix_and_counts(
                div_matrix, array_seqs[i], array_seqs[j]
            )
            if num_valid == 0:
                d = nan
            elif num_diffs == 0:
                d = zero
            else:
                d = _paralinear(
                    div_matrix,
                    num_states,
                )
            dists[i, j] = dists[j, i] = d

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return dists


def paralinear(
    aln: "Alignment", invalid_raises: bool = False, parallel: bool = False
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
    if not _is_nucleic(aln.moltype):
        raise new_moltype.MolTypeError(
            f"paralinear distance only works with nucleotide alignments, not {aln.moltype}"
        )

    num_states = len(aln.moltype.alphabet)

    mat = paralinear_distance_matrix(aln.array_seqs, num_states, parallel=parallel)
    if invalid_raises and numpy.isnan(mat).any():
        raise ArithmeticError("nan's in matrix")
    return DistanceMatrix.from_array_names(mat, aln.names)


@numba.jit(parallel=True)
def simple_distance_matrix(
    array_seqs, num_states, parallel=True, hamming=True
):  # pragma: no cover
    num_threads = numba.get_num_threads()
    if not parallel:
        numba.set_num_threads(1)

    n_seqs = array_seqs.shape[0]
    dists = numpy.zeros((n_seqs, n_seqs), dtype=numpy.float32)
    # making working matrices for each thread
    # outer parallel loop
    for i in numba.prange(n_seqs - 1):
        for j in range(i + 1, n_seqs):
            d, num_valid = num_diffs_and_valid(array_seqs[i], array_seqs[j], num_states)
            if not hamming:
                d /= num_valid
            dists[i, j] = dists[j, i] = d

    if not parallel:
        # restore original number of threads
        numba.set_num_threads(num_threads)
    return dists


def hamming(aln: "Alignment", parallel: bool = False, **kwargs) -> DistanceMatrix:
    """returns matrix of hamming distances

    Parameters
    ----------
    aln
        Alignment instance
    parallel
        If True, uses parallel processing via numba.
    """
    num_states = len(aln.moltype.alphabet)
    mat = simple_distance_matrix(
        aln.array_seqs, num_states, parallel=parallel, hamming=True
    )
    return DistanceMatrix.from_array_names(mat, aln.names)


def pdist(aln: "Alignment", parallel: bool = False, **kwargs) -> DistanceMatrix:
    """returns matrix of pairwise proportion difference distances

    Parameters
    ----------
    aln
        Alignment instance
    parallel
        If True, uses parallel processing via numba.
    """
    num_states = len(aln.moltype.alphabet)
    mat = simple_distance_matrix(
        aln.array_seqs, num_states, parallel=parallel, hamming=False
    )
    return DistanceMatrix.from_array_names(mat, aln.names)


_calculators = {
    "paralinear": paralinear,
    "jc69": jc69,
    "tn93": tn93,
    "hamming": hamming,
    "pdist": pdist,
}


def get_distance_calculator(name):
    """returns a pairwise distance calculator

    name is converted to lower case"""
    name = name.lower()
    if name not in _calculators:
        raise ValueError(f'Unknown pairwise distance calculator "{name}"')

    return _calculators[name]

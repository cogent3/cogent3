import json
import re

import numpy
import pytest

import cogent3

from cogent3._version import __version__
from cogent3.core import (
    new_alphabet,
    new_genetic_code,
    new_moltype,
    new_sequence,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance


@pytest.fixture(scope="function")
def dna_alphabet():
    return new_moltype.DNA.degen_gapped_alphabet


@pytest.fixture(scope="function")
def ascii_alphabet():
    return new_moltype.ASCII.alphabet


@pytest.mark.parametrize("name", ("dna", "rna", "protein", "protein_with_stop", "text"))
def test_moltype_make_seq(name):
    raw = "ACGGA"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


def test_moltype_make_bytes_seq():
    raw = "ACGGA"
    name = "bytes"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


@pytest.fixture
def bytes_alpha():
    return new_moltype.BYTES.most_degen_alphabet()


@pytest.fixture(scope="function")
def integer_seq(bytes_alpha):
    """Used for slicing tests"""
    return new_sequence.SeqView(seq="0123456789", alphabet=bytes_alpha)


@pytest.fixture
def ascii_alpha():
    return new_moltype.ASCII.most_degen_alphabet()


# Tests from test_sequence.py
@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
def test_seqview_initialisation(start, stop, step, bytes_alpha):
    """Initialising a SeqView should work with range of provided values"""
    seq_data = "0123456789"
    got = new_sequence.SeqView(
        seq=seq_data, start=start, stop=stop, step=step, alphabet=bytes_alpha
    )
    expected = seq_data[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("index", (-10, -5, 0, 5, 9))  # -10 and 9 are boundary
def test_seqview_index(index, bytes_alpha):
    """SeqView with default values can be sliced with a single index, when within the length of the sequence"""
    seq_data = "0123456789"
    sv = new_sequence.SeqView(seq=seq_data, alphabet=bytes_alpha)
    got = sv[index]
    expected = seq_data[index]
    assert got.str_value == expected
    assert len(got) == 1


def test_seqview_index_null(ascii_alpha):
    "Indexing a SeqView of length 0 should return an IndexError"
    sv = new_sequence.SeqView(seq="", alphabet=ascii_alpha)
    with pytest.raises(IndexError):
        _ = sv[0]


def test_seqview_step_0(bytes_alpha):
    "Initialising or slicing a SeqView with a step of 0 should return an IndexError"
    sv = new_sequence.SeqView(seq="0123456789", alphabet=bytes_alpha)
    with pytest.raises(ValueError):
        _ = sv[::0]
    with pytest.raises(ValueError):
        _ = new_sequence.SeqView(seq="0123456789", alphabet=bytes_alpha, step=0)


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index(start, bytes_alpha):
    "indexing out of bounds with a forward step should raise an IndexError"
    seq = "0123456789"
    length = abs(start - len(seq))
    pos_boundary_index = length
    neg_boundary_index = -length - 1

    sv = new_sequence.SeqView(seq=seq, start=start, alphabet=bytes_alpha)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index_positive_step_gt_1(start, bytes_alpha):
    "boundary condition for indexing out of bounds with a forward step greater than 1"
    seq = "0123456789"
    step = 2
    length = abs((start - len(seq)) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = new_sequence.SeqView(seq=seq, start=start, step=step, alphabet=bytes_alpha)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -11))
def test_seqview_invalid_index_reverse_step(stop, bytes_alpha):
    "boundary condition for indexing out of bounds with a reverse step"
    seq = "0123456789"
    step = -1
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = new_sequence.SeqView(
        seq=seq, start=start, stop=stop, step=step, alphabet=bytes_alpha
    )
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -6))
def test_seqview_invalid_index_reverse_step_gt_1(stop, bytes_alpha):
    "boundary condition for indexing out of bounds with a reverse step less than -1"
    seq = "0123456789"
    step = -2
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = new_sequence.SeqView(
        seq=seq, start=start, stop=stop, step=step, alphabet=bytes_alpha
    )
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


def test_seqview_slice_null(ascii_alpha):
    sv = new_sequence.SeqView(seq="", alphabet=ascii_alpha)
    assert len(sv) == 0
    got = sv[2:]
    assert len(got) == 0


def test_seqview_start_out_of_bounds(bytes_alpha):
    "boundary condition for start index out of bounds"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 1
    boundary = abs((init_start - init_stop) // init_step)
    sv = new_sequence.SeqView(
        seq=seq, start=init_start, stop=init_stop, step=init_step, alphabet=bytes_alpha
    )
    got = sv[boundary::].str_value
    assert got == ""


def test_seqview_start_out_of_bounds_step_gt_1(bytes_alpha):
    "boundary condition for start index out of bounds with step greater than 1"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 2
    boundary = abs((init_start - init_stop) // init_step)
    sv = new_sequence.SeqView(
        seq=seq, start=init_start, stop=init_stop, step=init_step, alphabet=bytes_alpha
    )
    got = sv[boundary::].str_value
    assert got == ""


def test_seqview_start_out_of_bounds_reverse_step(bytes_alpha):
    "boundary condition for start index out of bounds with reverse step"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, -2
    boundary_pos = abs((init_start - init_stop) // init_step)
    boundary_neg = -abs((init_start - init_stop) // init_step) - 1

    sv = new_sequence.SeqView(
        seq=seq, start=init_start, stop=init_stop, step=init_step, alphabet=bytes_alpha
    )

    assert sv[boundary_pos::].str_value == ""
    assert sv[boundary_neg::].str_value == ""


@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 3, None),
        slice(1, None, None),
        slice(1, 3, None),
        slice(None, None, None),
    ),
)
def test_seqview_defaults(simple_slices, bytes_alpha):
    """SeqView should accept slices with all combinations of default parameters"""
    seq = "0123456789"
    got = new_sequence.SeqView(seq=seq, alphabet=bytes_alpha)[simple_slices]
    expected = seq[simple_slices]
    assert got.str_value == expected


@pytest.mark.parametrize("index", (-8, -5, 0, 5, 8))
@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 10, None),
        slice(1, None, None),
        slice(1, 10, None),
        slice(1, 10, 1),
        slice(None, None, None),
    ),
)
def test_seqview_sliced_index(index, simple_slices, bytes_alpha):
    """SeqView that has been sliced with default parameters, can then be indexed"""
    seq = "0123456789"
    sv = new_sequence.SeqView(seq=seq, alphabet=bytes_alpha)
    got = sv[simple_slices][index]
    expected = seq[simple_slices][index]
    assert got.str_value == expected


@pytest.mark.parametrize("first_step", (1, 2, -1, -2))
@pytest.mark.parametrize("second_step", (1, 2, -1, -2))
def test_seqview_reverse_slice(first_step, second_step, bytes_alpha):
    """subsequent slices may reverse the previous slice"""
    seq = "0123456789"
    sv = new_sequence.SeqView(seq=seq, step=first_step, alphabet=bytes_alpha)
    got = sv[::second_step]
    expected = seq[::first_step][::second_step]
    assert got.str_value == expected


@pytest.mark.parametrize("seq", ("0123456789", "01234567890"))
@pytest.mark.parametrize("index", (-10, -4, 0, 6, 10))
@pytest.mark.parametrize("start", (None, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 9, -10, -11))
@pytest.mark.parametrize("step", (-1, -2))
def test_seqview_rev_sliced_index(index, start, stop, step, seq, bytes_alpha):
    """SeqView that has been reverse sliced, can then be sliced with a single index"""
    seq_data = seq
    try:  # if python slicing raises an index error, we expect SeqView to also throw error
        expected = seq_data[start:stop:step][index]
    except IndexError:
        with pytest.raises(IndexError):
            _ = new_sequence.SeqView(
                seq=seq_data, start=start, stop=stop, step=step, alphabet=bytes_alpha
            )[index].str_value
    else:  # if no index error, SeqView should match python slicing
        got = new_sequence.SeqView(
            seq=seq_data, start=start, stop=stop, step=step, alphabet=bytes_alpha
        )[index].str_value
        assert got == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_init_with_negatives(seq, start, stop, step, bytes_alpha):
    "SeqView initialisation should handle any combination of positive and negative slices"
    got = new_sequence.SeqView(
        seq=seq, start=start, stop=stop, step=step, alphabet=bytes_alpha
    )
    expected = seq[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_slice_with_negatives(seq, start, stop, step, bytes_alpha):
    """SeqView should handle any combination of positive and negative slices"""
    sv = new_sequence.SeqView(seq=seq, alphabet=bytes_alpha)
    got = sv[start:stop:step]
    expected = seq[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("start", (None, 0, 2))
@pytest.mark.parametrize("stop", (None, 5, 7, 10))
@pytest.mark.parametrize("step", (1, 2))
@pytest.mark.parametrize("start_2", (None, 0, 1, 2))
@pytest.mark.parametrize("stop_2", (None, 2, 4, 10))
@pytest.mark.parametrize("step_2", (1, 2))
def test_subsequent_slice_forward(
    start, stop, step, start_2, stop_2, step_2, bytes_alpha
):
    """SeqView should handle subsequent forward slice"""
    seq = "0123456789"
    sv = new_sequence.SeqView(seq=seq, alphabet=bytes_alpha)
    got = sv[start:stop:step][start_2:stop_2:step_2]
    expected = seq[start:stop:step][start_2:stop_2:step_2]
    assert got.str_value == expected
    assert len(got) == len(expected)


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first stop -ve
        (slice(None, -3, None), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -3, None), slice(None, -5, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -5, None), slice(None, -3, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, None)),
        # first stop -ve, second stop +ve, second slice OUTSIDE first
        (slice(None, -6, None), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, None), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, None)),
        # WITH FIRST STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 2), slice(None, -3, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -3, 2), slice(None, -2, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 2), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 2), slice(None, 3, None)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, 2), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, 2), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, None)),
        # WITH SECOND STEP > 1
        # first stop -ve
        (slice(None, -3, None), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -2, None), slice(None, -4, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -4, None), slice(None, -3, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, None), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, None), slice(None, -2, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, 3)),
        # WITH BOTH STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, 3)),
        # NON-ZERO START
        # first stop -ve
        (slice(1, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(1, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(1, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(1, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(1, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(1, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(1, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(1, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(1, 10, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(1, 6, 2), slice(None, -7, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
    ),
)
def test_subsequent_slice_neg_stop(slice_1, slice_2, ascii_alpha):
    """SeqView should handle subsequence slices with >=1 negative stop values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(seq=seq_data, alphabet=ascii_alpha)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first start -ve
        (slice(-6, None, None), slice(None, None, None)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-6, None, None), slice(-4, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-4, None, None), slice(-6, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-8, None, None), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(5, None, None), slice(-6, None, None)),
        # WITH FIRST STEP > 1
        # first start -ve
        (slice(-6, None, 2), slice(None, None, None)),
        # second start -ve
        (slice(None, None, 2), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, 2), slice(-6, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 2), slice(-9, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 2), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 2), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 2), slice(-9, None, None)),
        # WITH SECOND STEP > 1
        # first start -ve
        (slice(-6, None, None), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, None), slice(-6, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, None), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, None), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, None), slice(-9, None, 2)),
        # WITH BOTH STEP > 1
        # first start -ve
        (slice(-6, None, 3), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, 3), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-9, None, 3), slice(-7, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 3), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 3), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 3), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 3), slice(-9, None, 2)),
        (slice(-9, 7, 3), slice(-2, None, None)),
    ),
)
def test_subsequent_slice_neg_start(slice_1, slice_2, ascii_alpha):
    """SeqView should handle subsequence slices with >=1 negative start values,
    subsequent slices may or may not overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(seq=seq_data, alphabet=ascii_alpha)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first step -ve
        (slice(None, None, -1), slice(None, None, None)),
        # second step -ve
        (slice(None, None, None), slice(None, None, -1)),
        # both step -ve, start/stop -ve, second slice WITHIN first
        (slice(-1, -11, -2), slice(-1, -5, -3)),
        # both step -ve, start/stop -ve, second slice OUTSIDE first
        (slice(-1, -11, -2), slice(-1, -11, -3)),
        # both step -ve, start/stop +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(5, 0, -3)),
        # both step -ve, start/stop +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(10, 0, -3)),
        # first step -ve, second step +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(1, 5, 2)),
        # first step -ve, second step +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(0, 10, 2)),
        # first step +ve, second step -ve, second slice WITHIN first
        (slice(0, 10, 2), slice(4, 0, -2)),
        # first step +ve, second step -ve, second slice OUTSIDE first
        (slice(0, 10, 3), slice(10, 0, -2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-8, 11, 2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-19, 0, -2)),
    ),
)
def test_subsequent_slice_neg_step(slice_1, slice_2, ascii_alpha):
    """SeqView should handle subsequence slices with negative step values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "0123456789"
    sv = new_sequence.SeqView(seq=seq_data, alphabet=ascii_alpha)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "sub_slices_triple",
    (
        (slice(None, None, None), slice(None, None, None), slice(None, None, None)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 7, 1)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 9, 1)),
        (slice(1, 9, 1), slice(2, 8, 2), slice(3, 7, -3)),
    ),
)
def test_subslice_3(sub_slices_triple, ascii_alpha):
    """SeqView should handle three subsequent slices"""
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(seq=seq_data, alphabet=ascii_alpha)
    slice_1, slice_2, slice_3 = sub_slices_triple
    assert (
        sv[slice_1][slice_2][slice_3].str_value == seq_data[slice_1][slice_2][slice_3]
    )


@pytest.mark.parametrize("start", (0, 2, -1))
@pytest.mark.parametrize("stop", (7, 10, -11))
@pytest.mark.parametrize("step", (1, -2))
@pytest.mark.parametrize("start_2", (0, 2, -8))
@pytest.mark.parametrize("stop_2", (2, 4))
@pytest.mark.parametrize("step_2", (2, -1))
@pytest.mark.parametrize("start_3", (0, 1, -6))
@pytest.mark.parametrize("stop_3", (4, 10, -10))
@pytest.mark.parametrize("step_3", (2, -2))
def test_triple_slice(
    integer_seq, start, stop, step, start_2, stop_2, step_2, start_3, stop_3, step_3
):
    """SeqView should handle subsequent forward slice"""
    seq = integer_seq.seq
    got = integer_seq[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]
    expected = seq[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]

    assert got.str_value == expected
    assert len(got) == len(expected)


def test_seqview_repr():
    alpha = new_moltype.DNA.most_degen_alphabet()
    # Short sequence, defaults
    seq = "ACGT"
    view = new_sequence.SeqView(seq=seq, alphabet=alpha)
    expected = (
        "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=0, seqid=None, seq_len=4)"
    )
    assert repr(view) == expected

    # Long sequence
    seq = "ACGT" * 10
    view = new_sequence.SeqView(seq=seq, alphabet=alpha)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=0, stop=40, step=1, offset=0, seqid=None, seq_len=40)"
    assert repr(view) == expected

    # Non-zero start, stop, and step values
    seq = "ACGT" * 10
    view = new_sequence.SeqView(seq=seq, start=5, stop=35, step=2, alphabet=alpha)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=5, stop=35, step=2, offset=0, seqid=None, seq_len=40)"
    assert repr(view) == expected

    # offset
    seq = "ACGT"
    view = new_sequence.SeqView(seq=seq, offset=5, alphabet=alpha)
    expected = (
        "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=5, seqid=None, seq_len=4)"
    )
    assert repr(view) == expected

    # seqid
    seq = "ACGT"
    view = new_sequence.SeqView(seq=seq, seqid="seq1", alphabet=alpha)
    expected = "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=0, seqid='seq1', seq_len=4)"
    assert repr(view) == expected


# test_get_kmers_strict() to test_annotate_from_gff()
@pytest.fixture(scope="function")
def one_seq():
    return new_moltype.DNA.make_seq(seq="AACCTGGAACC")


def test_seq_repr(one_seq):
    pat = re.compile("[ACGT]+")
    expect = str(one_seq)
    seq = one_seq

    got = pat.findall(repr(seq))[0]
    assert expect.startswith(got), (expect, got)


def test_seq_repr_rc(one_seq):
    pat = re.compile("[ACGT]+")
    dna = one_seq.moltype
    expect = dna.rc(str(one_seq))
    seq = one_seq.rc()

    got = pat.findall(repr(seq))[0]
    assert expect.startswith(got), (expect, got)


def test_annotation_from_slice_with_stride():
    seq = new_moltype.DNA.make_seq(seq="AAACGCGCGAAAAAAA", name="s1")
    seq.add_feature(biotype="exon", name="ex1", spans=[(3, 9)])
    f = list(seq.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CGCGCG"
    s1 = seq[1::2]
    f = list(s1.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CCC"


def test_absolute_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.absolute_position(5)
    assert got == 5

    # an index outside the range of the sequence should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.absolute_position(20)

    with pytest.raises(IndexError):
        one_seq._seq.absolute_position(-20)


def test_absolute_position_positive(one_seq):
    # with an offset, the abs index should be offset + index
    one_seq.annotation_offset = 2
    got = one_seq._seq.absolute_position(2)
    assert got == 2 + 2

    # with an offset and start, the abs index should be offset + start + index
    view = one_seq[2::]
    view.annotation_offset = 2  # todo: do we want the annotation_offset to be preserved when slicing? I think yes
    got = view._seq.absolute_position(2)
    assert got == 2 + 2 + 2

    # with an offset, start and step, the abs index should be offset + start + index * step
    view = one_seq[2::2]
    view.annotation_offset = 2
    got = view._seq.absolute_position(2)
    assert got == 2 + 2 + 2 * 2


def test_relative_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.relative_position(5)
    assert got == 5

    # a -ve index  should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.relative_position(-5)


def test_relative_position(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView."""

    view = integer_seq[1:9:]
    # view = "12345678"
    got = view.relative_position(0)
    # precedes the view, so should return -1
    assert got == -1
    # exceeds the view, but still returns a value
    got = view.relative_position(10)
    assert got == 9


def test_relative_position_step_GT_one(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView with a step greater than one.
    """

    # precedes the view, with step > 1
    view = integer_seq[2:7:2]
    # view = "246", precedes the view by 1 step
    got = view.relative_position(0)
    assert got == -1
    # precedes the view by 0.5 step, default behaviour is to round up to 0
    got = view.relative_position(1)
    assert got == 0
    # exceeds the view by two steps, len(view) + 2 = 4
    got = view.relative_position(10)
    assert got == 4


@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("rev", (False, True))
def test_seqview_copy(sliced, rev, integer_seq):
    raw_data = integer_seq.seq
    integer_seq = integer_seq[::-1] if rev else integer_seq
    raw_data = raw_data[::-1] if rev else raw_data

    slice_start = 2
    slice_end = 4
    sv = integer_seq[slice_start:slice_end]
    copied = sv.copy(sliced=sliced)

    assert copied.str_value == raw_data[slice_start:slice_end]
    assert copied.is_reversed == integer_seq.is_reversed
    assert sliced and copied.seq is not sv.seq or copied.seq is integer_seq.seq


def test_relative_position_with_remainder(integer_seq):
    """tests relative_position when the index given is excluded from the view as it falls on
    a position that is 'stepped over'"""
    view = integer_seq[1:9:2]
    # view = "1357"
    got = view.relative_position(2)
    # 2 is stepped over in the view, so we return the index of 3 (which is 1)
    assert got == 1

    # setting the arg stop=True will adjust to the largest number, smaller than the given abs value, that is in the view
    got = view.relative_position(8, stop=True)
    # 8 is excluded from the view, so we return the index of 7 (which is 3)
    assert got == 3


@pytest.mark.parametrize("value", (0, 3))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, 1, 2))
@pytest.mark.parametrize("stop", (None, 10, 11))
@pytest.mark.parametrize("step", (None, 1, 2))
def test_absolute_relative_roundtrip(one_seq, value, offset, start, stop, step):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    view = one_seq[start:stop:step]
    view.annotation_offset = offset or 0
    abs_val = view._seq.absolute_position(value)
    rel_val = view._seq.relative_position(abs_val)
    assert rel_val == value


@pytest.mark.parametrize("value", (0, 2))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, -1, -2))
@pytest.mark.parametrize("stop", (None, -10))
@pytest.mark.parametrize("step", (-1, -2))
def test_absolute_relative_roundtrip_reverse(
    integer_seq, value, offset, start, stop, step
):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    view = integer_seq[start:stop:step]
    view.offset = offset or 0
    abs_val = view.absolute_position(value)
    rel_val = view.relative_position(abs_val)
    assert view.offset == (offset or 0)
    assert (view[rel_val]).str_value == view[value].str_value


def test_annotate_gff_nested_features(DATA_DIR):
    """correctly annotate a sequence with nested features"""
    # the synthetic example
    #          1111111111222222222333333333334
    # 1234567890123456789012345678901234567890
    #  **** biological_region
    #                                     ** biological_region
    #                                       * biological_region
    #      *******************************  gene
    #         *********************   mRNA
    #            *********            exon
    #                       *****     exon
    # ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC...
    seq = new_moltype.DNA.make_seq(
        seq="ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22"
    )
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    # TODO: directly assign an annotation_db, annotate_from_gff to be discontinued
    seq.annotate_from_gff(gff3_path)
    # we have 8 records in the gff file
    assert seq.annotation_db.num_matches() == 8

    # get the gene and check it has a single annotation and that
    # its slice is correct
    ann = list(seq.get_features(biotype="gene"))
    assert len(ann) == 1
    ann_seq = ann[0].get_slice()
    assert str(ann_seq) == "GGAAAATTTTTTTTTAAGGGGGAAAAAAAAA"
    # the gene has 1 transcript
    gene = ann[0]
    mrna = list(gene.get_children(biotype="mRNA"))
    assert len(mrna) == 1
    mrna = mrna[0]
    ann_seq = mrna.get_slice()
    assert str(ann_seq) == "AAAATTTTTTTTTAAGGGGGAAA"

    # the transcript has 2 exons, from the parent feature
    exons = list(mrna.get_children(biotype="exon"))
    assert len(exons) == 2
    # or the sequence
    ann = list(seq.get_features(biotype="exon"))
    assert len(ann) == 2
    exon_seqs = ("TTTTTTTTT", "GGGGG")
    assert tuple(str(ex.get_slice()) for ex in exons) == exon_seqs


@pytest.mark.xfail(
    reason="AttributeError: 'MolType' object has no attribute 'coerce_str'"
)
def test_to_moltype_dna():
    """to_moltype("dna") ensures conversion from T to U"""
    seq = new_moltype.DNA.make_seq(seq="AAAAGGGGTTT", name="seq1")
    rna = seq.to_moltype("rna")

    assert "T" not in rna


@pytest.mark.xfail(
    reason="AttributeError: 'MolType' object has no attribute 'coerce_str'"
)
def test_to_moltype_rna():
    """to_moltype("rna") ensures conversion from U to T"""
    seq = new_moltype.RNA.make_seq(seq="AAAAGGGGUUU", name="seq1")
    rna = seq.to_moltype("dna")

    assert "U" not in rna


def test_to_rich_dict():
    """Sequence to_dict works"""
    dna = new_moltype.DNA
    r = dna.make_seq(seq="AAGGCC", name="seq1")
    got = r.to_rich_dict()
    seq = new_sequence.SeqView(
        seq="AAGGCC", seqid="seq1", alphabet=dna.most_degen_alphabet()
    ).to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": r.moltype.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


def test_to_json():
    """to_json roundtrip recreates to_dict"""
    dna = new_moltype.DNA
    r = dna.make_seq(seq="AAGGCC", name="seq1")
    got = json.loads(r.to_json())
    seq = new_sequence.SeqView(
        seq="AAGGCC", seqid="seq1", alphabet=dna.most_degen_alphabet()
    ).to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": dna.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


@pytest.mark.xfail(
    reason="NotImplementedError: deserialising 'cogent3.core.new_sequence.DnaSequence' from json"
)
def test_offset_with_multiple_slices(DATA_DIR):
    seq = new_moltype.DNA.make_seq(
        seq="ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22"
    )
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    # TODO: directly assign an annotation_db, annotate_from_gff to be discontinued
    seq.annotate_from_gff(gff3_path)
    rd = seq[2:].to_rich_dict()
    s1 = deserialise_object(rd)
    assert s1.annotation_offset == 2
    rd = s1[3:].to_rich_dict()
    s2 = deserialise_object(rd)
    assert s2.annotation_offset == 5
    expect = {(f.seqid, f.biotype, f.name) for f in seq.get_features(start=5)}
    got = {(f.seqid, f.biotype, f.name) for f in s2.get_features()}
    assert got == expect


@pytest.mark.parametrize("coord", ("start", "stop"))
def test_seqview_to_rich_dict(coord, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)
    plus = sv.to_rich_dict()
    minus = sv[::-1].to_rich_dict()
    plus = plus.pop("init_args")
    minus = minus.pop("init_args")
    assert plus.pop("seq") == minus.pop("seq")
    assert plus["step"] == -minus["step"]
    assert coord not in plus
    assert coord not in minus


@pytest.mark.parametrize("reverse", (False, True))
def test_seqview_round_trip(reverse, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)
    sv = sv[::-1] if reverse else sv

    rd = sv.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_sequence.SeqView)
    assert got.to_rich_dict() == sv.to_rich_dict()


@pytest.mark.parametrize("reverse", (False, True))
def test_sliced_seqview_rich_dict(reverse, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sl = slice(2, 13)
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)[sl]
    sv = sv[::-1] if reverse else sv
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seq"] == parent[sl]
    assert rd["init_args"]["offset"] == 2


@pytest.mark.parametrize(
    "sl",
    (
        slice(2, 5, 1),  # positive indices, positive step
        slice(-8, -5, 1),  # negative indices, positive step
        slice(4, 1, -1),  # positive indices, negative step
        slice(-6, -9, -1),  # negative indices, negative step
    ),
)
@pytest.mark.parametrize("offset", (4, 0))
def test_parent_start_stop(sl, offset, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = "234" if sl.step > 0 else "432"
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv.offset = offset
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (2 + offset, 5 + offset)


@pytest.mark.parametrize(
    "sl",
    (
        slice(None, None, 1),  # slice whole sequence plus strand
        slice(None, None, -1),  # slice whole sequence minus strand
    ),
)
def test_parent_start_stop_limits(sl, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 10)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_start_stop_empty(rev, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = ""
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[0 : 0 : -1 if rev else 1]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 0)


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("index", range(9))
def test_parent_start_stop_singletons(index, rev, ascii_alphabet):
    data = "0123456789"
    start, stop = (-(10 - index), -(10 - index + 1)) if rev else (index, index + 1)
    sl = slice(start, stop, -1 if rev else 1)
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (index, index + 1)


def test_get_drawable(DATA_DIR):
    seq = cogent3.load_seq(DATA_DIR / "annotated_seq.gb")
    seq = seq[2000:4000]
    biotypes = "CDS", "gene", "mRNA"
    for feat in seq.get_features(biotype=biotypes, allow_partial=True):
        draw = feat.get_drawable()
        assert "(incomplete)" in draw.text

    full = seq.get_drawable(biotype=biotypes)
    # should only include elements that overlap the segment
    assert len(full.traces) == len(biotypes)
    # and their names should indicate they're incomplete
    for trace in full.traces:
        assert "(incomplete)" in trace.text


@pytest.mark.parametrize("gc,seq", ((1, "TCCTGA"), (1, "ACGTAA---"), (2, "TCCAGG")))
def test_has_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert seq.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_has_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert not seq.has_terminal_stop(gc=gc)


def test_has_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.has_terminal_stop(gc=gc, strict=True)


@pytest.mark.parametrize(
    "gc,seq",
    (
        (2, "TCCAGG"),
        (1, "TAATGA"),
        (1, "ACGTGA---"),
        (1, "--AT-CTGA"),
    ),
)
def test_trim_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    expect = re.sub("(TGA|AGG)(?=[-]*$)", "---" if "-" in seq else "", seq)

    seq = new_moltype.DNA.make_seq(seq=seq)
    got = str(seq.trim_stop_codon(gc=gc))
    assert got == expect


@pytest.mark.parametrize("gc,seq", ((1, "T?CTGC"), (2, "TCCAAG")))
def test_trim_terminal_stop_nostop(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    got = seq.trim_stop_codon(gc=gc)
    assert str(got) == str(seq)
    # since there's no stop, we just return the same object
    assert got is seq


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_trim_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert str(seq.trim_stop_codon(gc=gc)) == str(seq)


def test_trim_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.trim_stop_codon(gc=gc, strict=True)


@pytest.mark.parametrize("cast", (int, numpy.int32, numpy.int64, numpy.uint8))
def test_index_a_seq(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    got = seq[cast(1)]
    assert isinstance(got, new_sequence.Sequence)


@pytest.mark.parametrize("cast", (float, numpy.float32))
def test_index_a_seq_float_fail(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    index = cast(1)
    with pytest.raises(TypeError):
        seq[index]  # pylint: disable=W0104


@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_same_moltype(moltype):
    moltype = new_moltype.get_moltype(moltype)
    seq = moltype.make_seq(seq="TCCAG")
    got = seq.to_moltype(moltype)
    assert got is seq


def test_gapped_by_map_segment_iter():
    moltype = new_moltype.DNA
    m, seq = moltype.make_seq(seq="-TCC--AG").parse_out_gaps()
    g = list(seq.gapped_by_map_segment_iter(m, allow_gaps=True, recode_gaps=False))
    assert g == ["-", "TCC", "--", "AG"]


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("start_stop", ((None, None), (3, 7)))
def test_copied_parent_coordinates(sliced, rev, start_stop):
    orig_name = "orig"
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name=orig_name)
    start, stop = start_stop
    start = start or 0
    stop = stop or len(seq)
    sl = slice(start, stop)
    seq = seq[sl]
    sliced_name = "sliced"
    seq.name = sliced_name
    assert seq.name == sliced_name
    seq = seq.rc() if rev else seq
    copied = seq.copy(sliced=sliced)
    assert copied.name == sliced_name
    # matches original
    assert copied.parent_coordinates() == seq.parent_coordinates()
    # and expected -- the coordinate name always reflects the underlying sequence
    assert copied.parent_coordinates() == (orig_name, start, stop, -1 if rev else 1)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_coordinates(one_seq, rev):
    seq = one_seq[1:1]
    seq = seq.rc() if rev else seq
    seq.name = "sliced"  # this assignment does not affect the
    # note that when a sequence has zero length, the parent seqid is None
    assert seq.parent_coordinates() == (None, 0, 0, 1)


@pytest.mark.parametrize("cls", (str, bytes))
def test_coerce_to_seqview_str_bytes(cls, dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    s = bytes(seq, "utf8") if cls == bytes else seq
    got = new_sequence._coerce_to_seqview(s, seqid, alphabet=dna_alphabet)
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_sequence(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_moltype.DNA.make_seq(seq=seq), seqid, alphabet=dna_alphabet
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_already_seqview(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_sequence.SeqView(seq=seq, alphabet=dna_alphabet),
        seqid,
        alphabet=dna_alphabet,
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_seqview_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", alphabet=dna_alphabet)
    assert sv.seqid is None

    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    assert sv.seqid == "seq1"


def test_seqview_rich_dict_round_trip_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == "seq1"

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid == "seq1"

    sv = new_sequence.SeqView(seq="ACGGTGGGAC", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] is None

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid is None


def test_seqview_slice_propagates_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    sliced_sv = sv[1:8:2]
    assert sliced_sv.seqid == "seq1"

    copied_sv = sliced_sv.copy(sliced=False)
    assert copied_sv.seqid == "seq1"

    copied_sliced_sv = sliced_sv.copy(sliced=True)
    assert copied_sliced_sv.seqid == "seq1"


def test_sequences_propogates_seqid():
    # creating a name Sequence propagates the seqid to the SeqView.
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name="seq1")
    assert seq._seq.seqid == "seq1"

    # renaming the Sequence doesnt change the seqid of the SeqView.
    seq.name = "seq2"
    assert seq.name == "seq2"
    assert seq._seq.seqid == "seq1"


@pytest.mark.xfail(reason="no SeqView dispatch for new_alphabet.to_indices")
def test_sequences_propogates_seqid_seqview():
    # creating a Sequence with a seqview does not change the seqid of the SeqView.
    seq = new_moltype.DNA.make_seq(
        seq=new_sequence.SeqView(
            seq="ACGGTGGGAC", seqid="parent_name", alphabet=dna_alphabet
        ),
        name="seq_name",
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid == "parent_name"

    # creating a Sequence with an unnamed seqview does not name the SeqView.
    seq = new_moltype.DNA.make_seq(
        new_sequence.SeqView(seq="ACGGTGGGAC"), name="seq_name"
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid is None


def test_make_seq_assigns_to_seqview():
    seq = new_moltype.DNA.make_seq(seq="ACGT", name="s1")
    assert seq.name == seq._seq.seqid == "s1"


def test_empty_seqview_translate_position(dna_alphabet):
    sv = new_sequence.SeqView(seq="", alphabet=dna_alphabet)
    assert sv.absolute_position(0) == 0
    assert sv.relative_position(0) == 0


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
@pytest.mark.parametrize("length", (1, 8, 999))
def test_seqview_seq_len_init(start, stop, step, length, dna_alphabet):
    # seq_len is length of seq when None
    seq_data = "A" * length
    sv = new_sequence.SeqView(
        seq=seq_data, start=start, stop=stop, step=step, alphabet=dna_alphabet
    )
    expect = len(seq_data)
    # Check property and slot
    assert sv.seq_len == expect
    assert sv._seq_len == expect


@pytest.mark.parametrize("seq, seq_len", [("A", 0), ("", 1), ("A", 2)])
def test_seqview_seq_len_mismatch(seq, seq_len, dna_alphabet):
    # If provided, seq_len must match len(seq)
    with pytest.raises(AssertionError):
        new_sequence.SeqView(seq=seq, seq_len=seq_len, alphabet=dna_alphabet)


def test_seqview_copy_propagates_seq_len(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq, alphabet=dna_alphabet)
    copied = sv.copy()
    assert copied.seq_len == len(seq)


def test_seqview_seq_len_modified_seq(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq, alphabet=dna_alphabet)

    sv.seq = "ATGC"  # this should not modify seq_len
    assert sv.seq_len == len(seq)


def test_sequence_str_bytes_array():
    data = "ACGGTGGGAC"
    seq = new_moltype.DNA.make_seq(seq=data)
    print(type(seq))
    assert str(seq) == data
    assert bytes(seq) == data.encode("utf8")
    assert numpy.array_equal(
        numpy.array(seq), new_moltype.DNA.alphabet.to_indices(data)
    )


@pytest.mark.xfail(reason="remove when Sequence.__str__ complements sequence")
@pytest.mark.parametrize("seq,rc", (("ATGTTT", False), ("AAACAT", True)))
def test_translation(seq, rc):
    seq = new_moltype.DNA.make_seq(seq=seq)
    if rc:
        seq = seq.rc()
    get_str = str(seq)
    assert get_str == "ATGTTT"
    aa = seq.get_translation()
    assert str(aa) == "MF"


def test_get_translation_include_stop():
    s = new_moltype.DNA.make_seq(seq="ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True)
    assert str(aa) == "I*L"


def test_get_translation_trim_stop():
    s = new_moltype.DNA.make_seq(seq="ATTTCCTGA", name="s1")
    aa = s.get_translation(trim_stop=True)
    assert str(aa) == "IS"
    # no effect on internal stops
    s = new_moltype.DNA.make_seq(seq="ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True, trim_stop=True)
    assert str(aa) == "I*L"

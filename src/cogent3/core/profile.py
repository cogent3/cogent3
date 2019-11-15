import numpy

from numpy import digitize
from numpy.random import random

from cogent3.maths.util import safe_log, safe_p_log_p
from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.util.misc import extend_docstring_from


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def validate_freqs_array(data, axis=None):
    if (data < 0).any():
        raise ValueError("negative frequency not allowed")

    # we explicitly ignore nan
    result = data.sum(axis=axis)
    if not numpy.allclose(result[numpy.isnan(result) == False], 1):
        raise ValueError("invalid frequencies, sum(axis=1) is not equal to 1")


class _MotifNumberArray(DictArray):
    def __init__(self, data, motifs, row_indices=None, dtype=None):
        """
        data
            series of numbers, can be numpy array, CategoryCounter, dict instances
        row_indices
            row_indices correspond to original indexes, defaults to length of
            motif
        """
        # todo validate that motifs are strings and row_indices are ints or
        # strings
        # todo change row_indices argument name to row_keys
        if isinstance(data, numpy.ndarray):
            some_data = data.any()
        else:
            some_data = any(data)

        if not some_data or len(data) == 0:
            raise ValueError("Must provide data")

        try:
            len(data[0])
        except TypeError:
            ndim = 1
        else:
            ndim = 2
        num_elements = len(data) if ndim == 1 else len(data[0])
        if num_elements != len(motifs):
            raise ValueError(f"number of data elements {len(data[0])} != {len(motifs)}")
        motifs = tuple(motifs)

        # create template
        if row_indices is None and ndim == 2:
            row_indices = len(data)

        template = DictArrayTemplate(row_indices, motifs)
        darr = template.wrap(data)
        try:
            darr.array.astype(dtype, casting="safe")
        except TypeError as err:
            raise ValueError(err)
        self.__dict__.update(darr.__dict__)
        self.motifs = motifs
        self.motif_length = len(motifs[0])

    def __getitem__(self, names):
        (index, remaining) = self.template.interpret_index(names)
        result = self.array[index]
        row_indices = None
        motifs = self.motifs
        if remaining is not None:
            if type(names) != tuple:
                # we're indexing a row, motifs unchanged
                row_indices = None
                if result.ndim == 2:
                    row_indices = remaining.names[0]
            elif type(names[0]) == type(names[1]) == slice:
                # slicing rows and, motifs
                row_indices, motifs = remaining.names
            elif type(names[0]) == slice:
                # slicing rows, indexing a motif
                row_indices = remaining.names[0]
                motifs = [names[1]]
                result = result.reshape(result.shape[0], 1)
            elif type(names[1]) == slice:
                # slicing motifs, indexing rows
                row_indices = None
            else:
                raise RuntimeError(result)
            result = self.__class__(result, motifs, row_indices=row_indices)
        return result

    def take(self, indices, negate=False, axis=1):
        """returns new array

        Parameters
        ----------
        indices
            ints or keys corresponding to row (axis=0) col (axis=1)
        negate
            excludes the indicated indices from the result
        axis
            indicates row/column
        """
        assert 0 <= axis <= 1, "invalid axis"
        try:
            indices[0]
        except TypeError:
            raise ValueError("must provide indexable series to take")

        one_dim = self.array.ndim == 1
        if one_dim:
            current = self.template.names[0]
            axis = None
        else:
            current = self.template.names[1] if axis else self.template.names[0]

        # are indices a subset of of indicated axis
        if not set(indices) <= set(current):
            if (
                isinstance(indices[0], int)
                and 0 <= min(indices)
                and max(indices) < len(current)
            ):
                current = list(range(len(current)))
            elif isinstance(indices[0], int):
                raise IndexError(f"{indices} out of bounds")
            else:
                raise ValueError("unexpected indices")

        if negate:
            indices = tuple(v for v in current if v not in indices)

        indices = [current.index(v) for v in indices]

        result = self.array.take(indices, axis=axis)
        if one_dim:
            motifs = numpy.take(self.template.names[0], indices)
            row_order = None
        else:
            motifs = self.template.names[1]
            row_order = self.template.names[0]
            if axis:
                motifs = numpy.take(motifs, indices)
            else:
                row_order = numpy.take(row_order, indices)

        return self.__class__(result, motifs=motifs, row_indices=row_order)


def _get_ordered_motifs_from_tabular(data, index=1):
    """backend motif extraction function for motif_counts, motif_freqs and pssm
       assumed index 1 are motif strings; motif returned in order of occurrence"""

    chars = []
    for entry in data:
        if not entry[index] in chars:
            chars.append(entry[index])
    return chars


def _get_data_from_tabular(tab_data, motifs, dtype):
    """backend data extraction function for motif_counts, motif_freqs and pssm"""
    num_motifs = len(motifs)
    num_pos = len(tab_data) // num_motifs
    result = numpy.zeros((num_pos, num_motifs), dtype=dtype)
    for pos, motif, value in tab_data:
        motif_index = motifs.index(motif)
        result[pos, motif_index] = value
    return result


def make_motif_counts_from_tabular(tab_data):
    """converts data to indicated motif number class

    Parameters
    ----------
    tab_data : numpy array
       tab_data is numpy array, with tab_data.shape must be (n, 3)
   """
    motif = _get_ordered_motifs_from_tabular(tab_data)
    data = _get_data_from_tabular(tab_data, motif, "int")
    return MotifCountsArray(data, motif)


@extend_docstring_from(make_motif_counts_from_tabular)
def make_motif_freqs_from_tabular(tab_data):
    motif = _get_ordered_motifs_from_tabular(tab_data)
    data = _get_data_from_tabular(tab_data, motif, "float")
    return MotifFreqsArray(data, motif)


@extend_docstring_from(make_motif_counts_from_tabular)
def make_pssm_from_tabular(tab_data):
    motif = _get_ordered_motifs_from_tabular(tab_data)
    data = _get_data_from_tabular(tab_data, motif, "float")
    return PSSM(data, motif)


class MotifCountsArray(_MotifNumberArray):
    def __init__(self, counts, motifs, row_indices=None):
        super(MotifCountsArray, self).__init__(counts, motifs, row_indices, dtype=int)

    def _to_freqs(self, pseudocount=0):
        data = self.array
        if pseudocount:
            data = data + pseudocount
        row_sum = data.sum(axis=1)
        freqs = data / numpy.vstack(row_sum)
        return freqs

    def to_freq_array(self, pseudocount=0):
        """
        Parameters
        ----------
        pseudocount
            added to every element prior to normalising

        Returns
        -------
        a MotifFreqsArray
        """
        freqs = self._to_freqs(pseudocount=pseudocount)
        return MotifFreqsArray(
            freqs, self.template.names[1], row_indices=self.template.names[0]
        )

    def to_pssm(self, background=None, pseudocount=0):
        """returns a PSSM array

        Parameters
        ----------
        background
            array of numbers representing the background frequency distribution
        pseudocount
            added to every element prior to normalising
        """
        # make freqs, then pssm
        freqs = self._to_freqs(pseudocount=pseudocount)

        return PSSM(
            freqs,
            self.motifs,
            row_indices=self.template.names[0],
            background=background,
        )

    def motif_totals(self):
        """returns totalled motifs"""
        col_sums = self.array.sum(axis=0)
        return MotifCountsArray(col_sums, self.motifs)

    def row_totals(self):
        """returns totalled row values"""
        row_sums = self.array.sum(axis=1)
        template = DictArray(1, row_sums.shape[0])
        return template.wrap(row_sums)


class MotifFreqsArray(_MotifNumberArray):
    def __init__(self, data, motifs, row_indices=None):
        super(MotifFreqsArray, self).__init__(data, motifs, row_indices, dtype=float)
        axis = 0 if self.array.ndim == 1 else 1
        validate_freqs_array(self.array, axis=axis)

    def entropy(self):
        """Shannon entropy per position using log2"""
        entropies = safe_p_log_p(self.array)
        return entropies.sum(axis=1)

    def information(self):
        """returns information as -max_entropy - entropy"""
        n = self.shape[1]
        max_val = -numpy.log2(1 / n)
        return max_val - self.entropy()

    def simulate_seq(self):
        """returns a simulated sequence as a string"""
        cumsum = self.array.cumsum(axis=1)
        series = [
            digitize(random(), cumsum[i], right=False) for i in range(self.shape[0])
        ]
        seq = "".join([self.motifs[i] for i in series])
        return seq

    def to_pssm(self, background=None):
        """returns a PSSM array

        Parameters
        ----------
        background
            array of numbers representing the background frequency distribution
        """
        return PSSM(
            self.array,
            self.motifs,
            row_indices=self.template.names[0],
            background=background,
        )


class PSSM(_MotifNumberArray):
    """position specific scoring matrix

    A log-odds matrix"""

    def __init__(self, data, motifs, row_indices=None, background=None):
        data = numpy.array(data)
        row_sum = data.sum(axis=1)

        # are we dealing with counts data?
        if 0 <= data.min() and 1 < data.max():
            # convert to freqs data
            data = data / numpy.vstack(row_sum)
            row_sum = data.sum(axis=1)

        # are we dealing with freqs data?
        if (data >= 0).all() and numpy.allclose(
            row_sum[numpy.isnan(row_sum) == False], 1
        ):
            # standard PSSM object creation
            if background is None:
                background = numpy.ones(len(motifs), dtype=float) / len(motifs)
            self._background = numpy.array(background)
            assert len(background) == len(
                motifs
            ), "Mismatch between number of motifs and the background"
            validate_freqs_array(self._background)
            pssm = safe_log(data) - safe_log(self._background)
            super(PSSM, self).__init__(
                pssm, motifs, row_indices=row_indices, dtype=float
            )
            self._indices = numpy.arange(self.shape[0])  # used for scoring
            return

        if not (data.min() < 0 < data.max()):
            raise ValueError("PSSM has been supplied invalid data")

        # we dealing with pssm data
        super(PSSM, self).__init__(data, motifs, row_indices=row_indices, dtype=float)
        self._indices = numpy.arange(self.shape[0])  # used for scoring

    def get_indexed_seq(self, seq):
        """converts seq to numpy array of int
        characters in seq not present in motifs are assigned out-of-range index
        """
        get_index = {c: i for i, c in enumerate(self.motifs)}.get
        num_motifs = len(self.motifs)
        if self.motif_length == 1:
            indexed = [get_index(c, num_motifs) for c in seq]
        else:
            indexed = []
            for i in range(0, self.shape[0] - self.motif_length + 1, self.motif_length):
                indexed.append(get_index(seq[i : i + self.motif_length], num_motifs))
        indexed = numpy.array(indexed)
        return indexed

    def score_seq(self, seq):
        """return score for a sequence"""
        indexed = self.get_indexed_seq(seq)
        return self.score_indexed_seq(indexed)

    def score_indexed_seq(self, indexed):
        """return score for a sequence already converted to integer indices"""
        if len(indexed) < self.shape[1]:
            msg = f"sequence length {len(indexed)} shorter than PSSM {self.shape[1]}"
            raise ValueError(msg)
        indexed = numpy.array(indexed)
        num_motifs = len(self.motifs)
        scores = []

        for i in range(0, indexed.shape[0] - self.shape[0] + 1):
            segment = indexed[i : i + self.shape[0]]
            un_ambig = segment < num_motifs
            if not un_ambig.all():
                pssm = self.array[un_ambig, :]
                segment = segment[un_ambig]
            else:
                pssm = self.array
            score = pssm[self._indices[: pssm.shape[0]], segment].sum()
            scores.append(score)
        return scores

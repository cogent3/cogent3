from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.maths.util import safe_p_log_p, safe_log
import numpy
from numpy.random import random
from numpy import digitize

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def validate_freqs_array(data, axis=None):
    if (data < 0).any():
        raise ValueError('negative frequency not allowed')

    if not numpy.allclose(data.sum(axis=axis), 1):
        raise ValueError(
            'invalid frequencies, sum(axis=1) is not equal to 1')


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
            darr.array.astype(dtype, casting='safe')
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
        assert 0 <= axis <= 1, 'invalid axis'
        try:
            indices[0]
        except TypeError:
            raise ValueError('must provide indexable series to take')

        one_dim = self.array.ndim == 1
        if one_dim:
            current = self.template.names[0]
            axis = None
        else:
            current = self.template.names[1] if axis else self.template.names[0]

        # are indices a subset of of indicated axis
        if not set(indices) <= set(current):
            if (isinstance(indices[0], int) and 0 <= min(indices) and
                    max(indices) < len (current)):
                current = list(range(len(current)))
            elif isinstance(indices[0], int):
                raise IndexError(f'{indices} out of bounds')
            else:
                raise ValueError('unexpected indices')

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

class MotifCountsArray(_MotifNumberArray):
    def __init__(self, counts, motifs, row_indices=None):
        super(MotifCountsArray, self).__init__(counts, motifs, row_indices,
                                               dtype=int)

    def _to_freqs(self):
        row_sum = self.array.sum(axis=1)
        freqs = self.array / numpy.vstack(row_sum)
        return freqs

    def to_freq_array(self):
        """returns a MotifFreqsArray"""
        freqs = self._to_freqs()
        return MotifFreqsArray(freqs, self.template.names[1],
                               row_indices=self.template.names[0])

    def to_pssm(self, background=None):
        """returns a PSSM array

        Parameters
        ----------
        background
            array of numbers representing the background frequency distribution
        """
        # make freqs, then pssm
        freqs = self._to_freqs()

        return PSSM(freqs, self.motifs, row_indices=self.template.names[0],
                    background=background)

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
        super(MotifFreqsArray, self).__init__(data, motifs, row_indices,
                                              dtype=float)
        validate_freqs_array(self.array, axis=1)

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
        series = [digitize(random(), cumsum[i], right=False)
                  for i in range(self.shape[0])]
        seq = ''.join([self.motifs[i] for i in series])
        return seq


class PSSM(_MotifNumberArray):
    """position specific scoring matrix

    A log-odds matrix"""

    def __init__(self, data, motifs, row_indices=None, background=None):
        freqs = MotifFreqsArray(data, motifs, row_indices=row_indices)
        if background is None:
            background = numpy.ones(len(motifs), dtype=float) / len(motifs)
        self._background = numpy.array(background)
        assert len(background) == len(motifs), \
            "Mismatch between number of motifs and the background"
        validate_freqs_array(self._background)
        pssm = safe_log(freqs.array) - safe_log(self._background)
        super(PSSM, self).__init__(pssm, motifs, row_indices,
                                   dtype=float)
        self._indices = numpy.arange(self.shape[0])  # used for scoring

    def score_seq(self, seq):
        """return score for a sequence"""
        get_index = self.motifs.index
        if self.motif_length == 1:
            indexed = list(map(get_index, seq))
        else:
            indexed = []
            for i in range(0, self.shape[0] - self.motif_length, self.motif_length):
                indexed.append(get_index(seq[i: i + self.motif_length]))
        indexed = numpy.array(indexed)
        return self.score_indexed_seq(indexed)

    def score_indexed_seq(self, indexed):
        """return score for a sequence already converted to integer indices"""
        indexed = numpy.array(indexed)
        scores = []
        for i in range(0, indexed.shape[0] - self.shape[0] + 1):
            segment = indexed[i: i + self.shape[0]]
            score = self.array[self._indices, segment].sum()
            scores.append(score)
        return scores


import typing
from collections import defaultdict, namedtuple
from numbers import Number

import numba
import numpy
from numpy import array, diag, dot, eye, float64, int32, log, sqrt, zeros
from numpy.linalg import det, inv

import cogent3
from cogent3._version import __version__
from cogent3.core import moltype as c3_moltype
from cogent3.core.table import Table
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.dict_array import DictArray
from cogent3.util.misc import get_object_provenance
from cogent3.util.progress_display import display_wrap

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.tree import PhyloNode

PySeq = typing.Sequence
PySeqStr = PySeq[str]


@numba.jit(cache=True)
def fill_diversity_matrix(matrix, seq1, seq2) -> None:  # pragma: no cover
    """fills the diversity matrix for valid positions.

    Assumes the provided sequences have been converted to indices with
    invalid characters being negative numbers (use get_moltype_index_array
    plus seq_to_indices)."""

    for i in range(len(seq1)):
        if seq1[i] < 0 or seq2[i] < 0:
            continue
        matrix[seq1[i], seq2[i]] += 1.0


def _same_moltype(ref, query):
    """if ref and query have the same states"""
    return set(ref) == set(query)


def get_pyrimidine_indices(moltype):
    """returns pyrimidine indices for the moltype"""
    alpha = moltype.alphabet
    if len(alpha) != 4:
        msg = "Non-nucleic acid MolType"
        raise RuntimeError(msg)
    pyrimidines = "CT" if moltype.label == "dna" else "CU"
    return alpha.to_indices(pyrimidines).tolist()


def get_purine_indices(moltype):
    """returns purine indices for the moltype"""
    alpha = moltype.alphabet
    if len(alpha) != 4:
        msg = "Non-nucleic acid MolType"
        raise RuntimeError(msg)
    return alpha.to_indices("AG").tolist()


def get_matrix_diff_coords(indices):
    """returns coordinates for off diagonal elements"""
    return [(i, j) for i in indices for j in indices if i != j]


def get_moltype_index_array(moltype, invalid=-9):
    """returns the index array for a molecular type"""
    canonical_chars = list(moltype)
    # maximum ordinal for an allowed character, this defines the length of
    # the required numpy array
    # refactor: simplify
    # added for compatibility with both new and old moltypes - should be removed
    # when old moltype is removed
    from cogent3.core.moltype import MolType as c3_moltype

    if isinstance(moltype, c3_moltype):
        max_ord = max(list(map(ord, list(moltype.most_degen_alphabet()))))
    else:
        max_ord = max(list(map(ord, list(moltype.All.keys()))))
    char_to_index = zeros(max_ord + 1, int32)
    # all non canonical_chars are ``invalid''
    char_to_index.fill(invalid)

    for i in range(len(canonical_chars)):
        c = canonical_chars[i]
        o = ord(c)
        char_to_index[o] = i

    return char_to_index


def seq_to_indices(seq, char_to_index):
    """returns an array with sequence characters replaced by their index"""
    ords = list(map(ord, seq))
    return char_to_index.take(ords)


def _hamming(matrix):
    """computes the edit distance
    Parameters
    ----------
    matrix : array
        2D numpy array of counts
    Returns
    -------
    total of the matrix, the proportion of changes, hamming distance, variance
    (the variance calculation is not yet implemented)
    """
    # TODO implement the estimate of the variance
    invalid = None, None, None, None
    total = matrix.sum()
    dist = total - diag(matrix).sum()
    if total == 0:
        return invalid

    p = dist / total

    return total, p, dist, None


def _jc69_from_matrix(matrix):
    """computes JC69 stats from a diversity matrix"""
    invalid = None, None, None, None
    total = matrix.sum()
    diffs = total - diag(matrix).sum()
    if total == 0:
        return invalid

    p = diffs / total
    if p >= 0.75:  # cannot take log
        return invalid

    factor = 1 - (4 / 3) * p

    dist = -3.0 * log(factor) / 4
    var = p * (1 - p) / (factor * factor * total)
    return total, p, dist, var


def _tn93_from_matrix(
    matrix,
    freqs,
    pur_indices,
    pyr_indices,
    pur_coords,
    pyr_coords,
    tv_coords,
):
    invalid = None, None, None, None

    total = matrix.sum()
    freqs = matrix.sum(axis=0) + matrix.sum(axis=1)
    freqs /= 2 * total

    if total == 0:
        return invalid

    p = matrix.take(pur_coords + pyr_coords + tv_coords).sum() / total

    freq_purs = freqs.take(pur_indices).sum()
    prod_purs = freqs.take(pur_indices).prod()
    freq_pyrs = freqs.take(pyr_indices).sum()
    prod_pyrs = freqs.take(pyr_indices).prod()

    # purine transition diffs
    pur_ts_diffs = matrix.take(pur_coords).sum()
    pur_ts_diffs /= total
    # pyr transition  diffs
    pyr_ts_diffs = matrix.take(pyr_coords).sum()
    pyr_ts_diffs /= total
    # transversions
    tv_diffs = matrix.take(tv_coords).sum() / total

    coeff1 = 2 * prod_purs / freq_purs
    coeff2 = 2 * prod_pyrs / freq_pyrs
    coeff3 = 2 * (
        freq_purs * freq_pyrs
        - (prod_purs * freq_pyrs / freq_purs)
        - (prod_pyrs * freq_purs / freq_pyrs)
    )

    term1 = 1 - pur_ts_diffs / coeff1 - tv_diffs / (2 * freq_purs)
    term2 = 1 - pyr_ts_diffs / coeff2 - tv_diffs / (2 * freq_pyrs)
    term3 = 1 - tv_diffs / (2 * freq_purs * freq_pyrs)

    if term1 <= 0 or term2 <= 0 or term3 <= 0:  # log will fail
        return invalid

    dist = -coeff1 * log(term1) - coeff2 * log(term2) - coeff3 * log(term3)
    v1 = 1 / term1
    v2 = 1 / term2
    v3 = 1 / term3
    v4 = (
        (coeff1 * v1 / (2 * freq_purs))
        + (coeff2 * v2 / (2 * freq_pyrs))
        + (coeff3 * v3 / (2 * freq_purs * freq_pyrs))
    )
    var = (
        v1**2 * pur_ts_diffs
        + v2**2 * pyr_ts_diffs
        + v4**2 * tv_diffs
        - (v1 * pur_ts_diffs + v2 * pyr_ts_diffs + v4 * tv_diffs) ** 2
    )
    var /= total

    return total, p, dist, var


def _logdetcommon(matrix):
    invalid = (None,) * 5

    total = matrix.sum()
    diffs = total - matrix.diagonal().sum()
    if total == 0:
        return invalid

    p = diffs / total

    if diffs == 0:  # seqs indentical
        return invalid

    # we replace the missing diagonal states with a frequency of 0.5,
    # then normalise
    frequency = matrix.copy()
    frequency[(frequency == 0) * eye(*matrix.shape, dtype=bool)] = 0.5
    frequency /= frequency.sum()

    if det(frequency) <= 0:  # if the result is nan
        return invalid

    # the inverse matrix of frequency, every element is squared
    M_matrix = inv(frequency) ** 2
    freqs = [frequency.sum(axis=axis) for axis in (0, 1)]
    var_term = dot(M_matrix, frequency).diagonal().sum()

    return total, p, frequency, freqs, var_term


def _paralinear(matrix):
    """the paralinear distance from a diversity matrix"""

    invalid = (None,) * 4

    total, p, frequency, freqs, var_term = _logdetcommon(matrix)
    if frequency is None:
        return invalid

    r = matrix.shape[0]
    d_xy = -log(det(frequency) / sqrt((freqs[0] * freqs[1]).prod())) / r
    var = (var_term - (1 / sqrt(freqs[0] * freqs[1])).sum()) / (r**2 * total)

    return total, p, d_xy, var


def _logdet(matrix, use_tk_adjustment=True):
    """returns the LogDet from a diversity matrix

    Parameters
    ----------
    use_tk_adjustment
        when True, unequal state frequencies are allowed

    """

    invalid = (None,) * 4

    total, p, frequency, freqs, var_term = _logdetcommon(matrix)
    if frequency is None:
        return invalid

    r = matrix.shape[0]
    if use_tk_adjustment:
        coeff = (sum(sum(freqs) ** 2) / 4 - 1) / (r - 1)
        d_xy = coeff * log(det(frequency) / sqrt((freqs[0] * freqs[1]).prod()))
        var = None
    else:
        d_xy = -log(det(frequency)) / r - log(r)
        var = (var_term / r**2 - 1) / total

    return total, p, d_xy, var


def _number_formatter(template):
    """flexible number formatter"""

    def call(val):
        try:
            result = template % val
        except TypeError:
            result = val
        return result

    return call


Stats = namedtuple("Stats", ["length", "fraction_variable", "dist", "variance"])


def _make_stat_table(stats, names, **kwargs) -> Table:
    from cogent3.core.table import Table

    header = ["Seq1 \\ Seq2", *names]
    rows = zeros((len(names), len(names)), dtype="O")
    for i in range(len(names) - 1):
        n1 = names[i]
        for j in range(i + 1, len(names)):
            n2 = names[j]
            val = stats[(n1, n2)]
            rows[i, j] = val
            rows[j, i] = val
    rows = rows.tolist()
    for i in range(len(names)):
        rows[i].insert(0, names[i])

    return Table(
        header=header,
        data=rows,
        index_name=r"Seq1 \ Seq2",
        missing_data="*",
        **kwargs,
    )


class _PairwiseDistance:
    """base class for computing pairwise distances"""

    valid_moltypes = ()

    def __init__(
        self,
        moltype,
        invalid=-9,
        alignment=None,
        invalid_raises=False,
    ) -> None:
        super().__init__()
        moltype = cogent3.get_moltype(moltype)
        if moltype.label not in self.valid_moltypes:
            name = self.__class__.__name__
            msg = (
                f"Invalid moltype for {name}: '{moltype.label}' not "
                f"in {self.valid_moltypes}"
            )
            raise ValueError(msg)

        self.moltype = moltype
        self.char_to_indices = get_moltype_index_array(moltype, invalid=invalid)
        self._dim = len(list(moltype))
        self._dists = None
        self._dupes = None
        self._duped = None
        self._invalid_raises = invalid_raises

        self.names = None
        self.indexed_seqs = None

        if alignment is not None:
            self._convert_seqs_to_indices(alignment)

        self._func_args = []

    def _convert_seqs_to_indices(self, alignment) -> None:
        assert isinstance(
            alignment.moltype,
            type(self.moltype) | c3_moltype.MolType,
        ), "Alignment does not have correct MolType"

        self._dists = {}
        self.names = alignment.names[:]
        indexed_seqs = []
        for name in self.names:
            seq = alignment.get_gapped_seq(name)
            indexed = seq_to_indices(str(seq), self.char_to_indices)
            indexed_seqs.append(indexed)

        self.indexed_seqs = array(indexed_seqs)

    @property
    def duplicated(self):
        """returns mapping of IDs to duplicates as {id:[dupe1, ..], },
        or None"""
        return self._duped

    @staticmethod
    def func() -> None:
        pass  # over ride in subclasses

    @display_wrap
    def run(self, alignment=None, ui=None) -> None:
        """computes the pairwise distances"""
        self._dupes = None
        self._duped = None

        dupes = set()
        duped = defaultdict(list)

        if alignment is not None:
            self._convert_seqs_to_indices(alignment)

        names = self.names[:]
        matrix = zeros((self._dim, self._dim), float64)
        off_diag = [
            (i, j) for i in range(self._dim) for j in range(self._dim) if i != j
        ]
        off_diag = tuple(tuple(a) for a in zip(*off_diag, strict=False))

        done = 0.0
        to_do = (len(names) ** 2 - 1) / 2
        for i in range(len(names) - 1):
            if i in dupes:
                continue

            name_1 = names[i]
            s1 = self.indexed_seqs[i]
            for j in range(i + 1, len(names)):
                if j in dupes:
                    continue

                name_2 = names[j]
                ui.display(f"{name_1} vs {name_2}", done / to_do)
                done += 1
                matrix.fill(0)
                s2 = self.indexed_seqs[j]
                fill_diversity_matrix(matrix, s1, s2)
                if not (matrix[off_diag] > 0).any():
                    # j is a duplicate of i
                    dupes.update([j])
                    duped[i].append(j)
                    continue

                total, p, dist, var = self.func(matrix, *self._func_args)
                if self._invalid_raises and not isinstance(dist, Number):
                    msg = f"distance could not be calculated for {name_1} - {name_2}"
                    raise ArithmeticError(msg)
                result = Stats(total, p, dist, var)
                self._dists[(name_1, name_2)] = result
                self._dists[(name_2, name_1)] = result

        self._dupes = [names[i] for i in dupes] or None
        if duped:
            self._duped = {}
            for k, v in duped.items():
                key = names[k]
                vals = [names[i] for i in v]
                self._duped[key] = vals

            # clean the distances so only unique seqs included
            remove = set(self._dupes)
            keys = list(self._dists.keys())
            for key in keys:
                if set(key) & remove:
                    del self._dists[key]

    __call__ = run

    def get_pairwise_distances(self, include_duplicates=True):
        """returns a matrix of pairwise distances.

        Parameters
        ----------
        include_duplicates : bool
            all seqs included in the distances, otherwise only unique sequences
            are included.
        """
        if self._dists is None:
            return None

        dists = {k: self._dists[k].dist for k in self._dists}
        if include_duplicates:
            dists = self._expand(dists)

        return DistanceMatrix(dists)

    def _expand(self, pwise):
        """returns a pwise statistic dict that includes duplicates"""
        if not self.duplicated:
            # no duplicates, nothing to do
            return pwise

        redundants = {}
        for k in self.duplicated:
            for r in self.duplicated[k]:
                redundants[r] = k

        names = self.names[:]
        for add, alias in redundants.items():
            for name in names:
                if name == add:
                    continue
                val = 0 if name == alias else pwise.get((alias, name), None)
                pwise[(add, name)] = pwise[(name, add)] = val

        return pwise

    @property
    def dists(self):
        if self._dists is None:
            return None

        return self.get_pairwise_distances(include_duplicates=True)

    @property
    def stderr(self):
        if self._dists is None:
            return None

        stats = {k: sqrt(self._dists[k].variance) for k in self._dists}
        stats = self._expand(stats)
        kwargs = {"title": "Standard Error of Pairwise Distances", "digits": 4}
        return _make_stat_table(stats, self.names, **kwargs)

    @property
    def variances(self):
        if self._dists is None:
            return None

        stats = {k: self._dists[k].variance for k in self._dists}
        stats = self._expand(stats)
        kwargs = {"title": "Variances of Pairwise Distances", "digits": 4}
        t = _make_stat_table(stats, self.names, **kwargs)
        var_formatter = _number_formatter("%.2e")
        for name in self.names:
            t.format_column(name, var_formatter)
        return t

    @property
    def proportions(self):
        if self._dists is None:
            return None

        stats = {k: self._dists[k].fraction_variable for k in self._dists}
        stats = self._expand(stats)
        kwargs = {"title": "Proportion variable sites", "digits": 4}
        return _make_stat_table(stats, self.names, **kwargs)

    @property
    def lengths(self):
        if self._dists is None:
            return None

        stats = {k: self._dists[k].length for k in self._dists}
        stats = self._expand(stats)
        kwargs = {"title": "Pairwise Aligned Lengths", "digits": 0}
        return _make_stat_table(stats, self.names, **kwargs)


class HammingPair(_PairwiseDistance):
    """Hamming distance calculator for pairwise alignments"""

    valid_moltypes = ("dna", "rna", "protein", "text", "bytes")

    def __init__(self, moltype="text", *args, **kwargs) -> None:
        """states: the valid sequence states"""
        super().__init__(moltype, *args, **kwargs)
        self.func = _hamming


class ProportionIdenticalPair(_PairwiseDistance):
    """Proportion sites identical distance calculator for pairwise alignments"""

    valid_moltypes = ("dna", "rna", "protein", "text", "bytes")

    def __init__(self, moltype="text", *args, **kwargs) -> None:
        """states: the valid sequence states"""
        super().__init__(moltype, *args, **kwargs)
        self.func = _hamming

    def get_pairwise_distances(self, include_duplicates=True):
        """returns a matrix of pairwise distances.

        Parameters
        ----------
        include_duplicates : bool
            all seqs included in the distances, otherwise only unique sequences
            are included.
        """
        if self._dists is None:
            return None

        dists = {k: self._dists[k].fraction_variable for k in self._dists}
        if include_duplicates:
            dists = self._expand(dists)

        return DistanceMatrix(dists)


class _NucleicSeqPair(_PairwiseDistance):
    """base class pairwise distance calculator for nucleic acid seqs"""

    valid_moltypes = ("dna", "rna")

    def __init__(self, moltype="dna", *args, **kwargs) -> None:
        super().__init__(moltype, *args, **kwargs)
        if not _same_moltype(
            cogent3.get_moltype("dna"),
            self.moltype,
        ) and not _same_moltype(
            cogent3.get_moltype("rna"),
            self.moltype,
        ):
            msg = "Invalid MolType for this metric"
            raise RuntimeError(msg)


class JC69Pair(_NucleicSeqPair):
    """JC69 distance calculator for pairwise alignments"""

    def __init__(self, moltype="dna", *args, **kwargs) -> None:
        """states: the valid sequence states"""
        super().__init__(moltype, *args, **kwargs)
        self.func = _jc69_from_matrix


class TN93Pair(_NucleicSeqPair):
    """TN93 calculator for pairwise alignments"""

    def __init__(self, moltype="dna", *args, **kwargs) -> None:
        """states: the valid sequence states"""
        super().__init__(moltype, *args, **kwargs)
        self._freqs = zeros(self._dim, float64)
        self.pur_indices = get_purine_indices(self.moltype)
        self.pyr_indices = get_pyrimidine_indices(self.moltype)

        # matrix coordinates
        self.pyr_coords = get_matrix_diff_coords(self.pyr_indices)
        self.pur_coords = get_matrix_diff_coords(self.pur_indices)

        self.tv_coords = get_matrix_diff_coords(list(range(self._dim)))
        for coord in self.pur_coords + self.pyr_coords:
            self.tv_coords.remove(coord)

        # flattened
        self.pyr_coords = [i * 4 + j for i, j in self.pyr_coords]
        self.pur_coords = [i * 4 + j for i, j in self.pur_coords]
        self.tv_coords = [i * 4 + j for i, j in self.tv_coords]

        self.func = _tn93_from_matrix
        self._func_args = [
            self._freqs,
            self.pur_indices,
            self.pyr_indices,
            self.pur_coords,
            self.pyr_coords,
            self.tv_coords,
        ]


class LogDetPair(_PairwiseDistance):
    """computes logdet distance between sequence pairs"""

    valid_moltypes = ("dna", "rna", "protein")

    def __init__(self, moltype="dna", use_tk_adjustment=True, *args, **kwargs) -> None:
        """Arguments:
        - moltype: string or moltype instance (must be dna or rna)
        - use_tk_adjustment: use the correction of Tamura and Kumar 2002
        """
        super().__init__(moltype, *args, **kwargs)
        self.func = _logdet
        self._func_args = [use_tk_adjustment]

    def run(self, use_tk_adjustment=None, *args, **kwargs) -> None:
        if use_tk_adjustment is not None:
            self._func_args = [use_tk_adjustment]

        super().run(*args, **kwargs)


class ParalinearPair(_PairwiseDistance):
    """computes the paralinear distance (Lake 1994) between sequence pairs"""

    valid_moltypes = ("dna", "rna", "protein")

    def __init__(self, moltype="dna", *args, **kwargs) -> None:
        super().__init__(moltype, *args, **kwargs)
        self.func = _paralinear


_calculators = {
    "paralinear": ParalinearPair,
    "logdet": LogDetPair,
    "jc69": JC69Pair,
    "tn93": TN93Pair,
    "hamming": HammingPair,
    "pdist": ProportionIdenticalPair,
}


def get_distance_calculator(name, *args, **kwargs):
    """returns a pairwise distance calculator

    name is converted to lower case"""
    name = name.lower()
    if "moltype" in kwargs and kwargs.get("moltype") is None:
        kwargs.pop("moltype")
    if name not in _calculators:
        msg = f'Unknown pairwise distance calculator "{name}"'
        raise ValueError(msg)

    calc = _calculators[name]
    return calc(*args, **kwargs)


def available_distances() -> Table:
    """returns Table listing available fast pairwise genetic distance calculator

    Notes
    -----
    For more complicated genetic distance methods, see the evolve.models module.
    """
    from cogent3.core.table import Table

    rows = []
    for n, c in _calculators.items():
        rows.append([n, ", ".join(c.valid_moltypes)])

    return Table(
        header=["Abbreviation", "Suitable for moltype"],
        data=rows,
        title=(
            "Specify a pairwise genetic distance calculator "
            "using 'Abbreviation' (case insensitive)."
        ),
        index_name="Abbreviation",
    )


class DistanceMatrix(DictArray):
    """pairwise distance matrix"""

    def __init__(
        self, dists, invalid: bool | None = None, source: str | None = None
    ) -> None:
        super().__init__(dists, dtype=float)

        self._invalid = invalid
        self.source = source

    @classmethod
    def from_array_names(
        cls, matrix: numpy.ndarray, names: PySeqStr, invalid=None
    ) -> "DistanceMatrix":
        """construct a distance matrix from numpy array and names"""
        darr = DictArray.from_array_names(matrix, names, names)
        return cls(darr, invalid=invalid)

    def __setitem__(self, names, value) -> None:
        index, _ = self.template.interpret_index(names)
        self.array[index] = value

    def __getitem__(self, names):
        index, remaining = self.template.interpret_index(names)
        result = self.array[index]
        if remaining is not None:
            result = self.__class__(result, remaining)
            self.template.names = array(self.template.names)[index]
            result.template = self.template
        return result

    @property
    def names(self) -> list[str]:
        return self.template.names[0]

    def to_table(self) -> Table:
        """converted to a Table"""
        from cogent3.core.table import Table

        data = {"names": self.names}
        for i, name in enumerate(self.names):
            column = self.array[:, i]
            data[name] = column
        header = ["names", *list(self.names)]
        return Table(header=header, data=data, index_name="names")

    def to_dict(self, **kwargs) -> dict[tuple[str, str], float]:
        """Returns a flattened dict with diagonal elements removed"""
        result = super().to_dict(flatten=True)
        for n1 in self.names:
            del result[(n1, n1)]
        return result

    def to_rich_dict(self) -> dict:
        # because dicts with tuples as keys cannot be json'ed, we convert to
        # a list of tuples
        dists = self.to_dict()
        json_safe = [(k[0], k[1], dists[k]) for k in dists]
        return {
            "dists": json_safe,
            "invalid": self._invalid,
            "type": get_object_provenance(self),
            "version": __version__,
        }

    def take_dists(
        self, names: list[str] | str, negate: bool = False
    ) -> "DistanceMatrix":
        """
        Parameters
        ----------
        names
            series of names
        negate : bool
            if True, elements in names will be excluded
        Returns
        -------
        DistanceMatrix for names x names
        """
        if isinstance(names, str):
            names = [names]

        if list(self.names) == list(names):
            return self

        current_names = array(self.names)
        if negate:
            keep = [self.names.index(n) for n in current_names if n not in names]
        else:
            keep = [self.names.index(n) for n in names if n in current_names]

        if len(keep) <= 1:
            return None

        data = self.array.take(keep, axis=0)
        data = data.take(keep, axis=1)
        names = current_names.take(keep)
        return self.from_array_names(data, names)

    def drop_invalid(self) -> "DistanceMatrix":
        """drops all rows / columns with an invalid entry"""
        if (
            self.shape[0] != self.shape[1]
            or self.template.names[0] != self.template.names[1]
        ):
            msg = "Must be a square matrix"
            raise RuntimeError(msg)
        names = array(self.names)
        # NaN is an invalid value
        cols = numpy.isnan(self.array).sum(axis=0)
        exclude = names[cols != 0].tolist()
        rows = numpy.isnan(self.array).sum(axis=1)
        exclude += names[rows != 0].tolist()
        exclude = set(exclude)
        keep = set(names) ^ exclude
        return self.take_dists(keep)

    def quick_tree(self, use_hook: str | None = None) -> "PhyloNode":
        """Returns a phylogenetic tree

        Parameters
        ----------
        use_hook
            name of a third-party package that implements the quick_tree
            hook. If not specified, defaults to the first available hook or
            the cogent3 quick_tree() app. To force default, set
            use_hook="cogent3".

        Notes
        -----
        Invalid distances are dropped prior to building the tree.
        Defaults to the Neighbour Joining Tree algorithm.
        If only two sequences are present, a simple tree is returned
        and hooks are not used.
        """
        qtree = cogent3._plugin.get_quick_tree_hook(name=use_hook)  # noqa: SLF001
        dists = self.drop_invalid()
        if dists.shape[0] == 2:
            # NJ not needed for 2 seqs
            dist = dists[0, 1]
            n1, n2 = dists.names
            dist = dists[n1, n2] / 2
            return cogent3.make_tree(f"({n1}:{dist},{n2}:{dist})")

        # we directly use the app main method, as this means any errors will be
        # raised as exceptions
        return qtree.main(dists)

    def max_pair(self) -> tuple[str, str]:
        """returns the pair of names with the maximum distance

        Returns
        -------
        tuple of the two names corresponding to the maximum distance

        Notes
        -----
        In case of multiple occurrences of the maximum values, the names
        corresponding to the first occurrence are returned.
        """
        max_index_flat = numpy.argmax(self)
        max_index_1, max_index_2 = numpy.unravel_index(max_index_flat, self.shape)
        return self.names[max_index_1], self.names[max_index_2]

    def min_pair(self) -> tuple[str, str]:
        """returns the pair of names with the minimum distance (excluding the diagonal)

        Returns
        -------
        tuple of the two names corresponding to the minimum distance

        Notes
        -----
        In case of multiple occurrences of the minimum values, the names
        corresponding to the first occurrence are returned.
        """
        dmat_copy = self.array.copy()
        numpy.fill_diagonal(
            dmat_copy,
            numpy.inf,
        )  # Exclude diagonal by setting diagonal elements to infinity
        min_index_flat = numpy.argmin(dmat_copy)
        min_index_1, min_index_2 = numpy.unravel_index(min_index_flat, self.shape)
        return self.names[min_index_1], self.names[min_index_2]


@register_deserialiser(get_object_provenance(DistanceMatrix))
def deserialise_tabular(data: dict) -> DistanceMatrix:
    """deserialising DistanceMatrix from a rich dict"""
    data.pop("version", None)
    _ = data.pop("type", None)
    # dists is a list of simple dists from which we reconstruct a 1D dict
    dists = {}
    for element in data["dists"]:
        key = tuple(element[:2])
        value = element[2]
        dists[key] = value
    data["dists"] = dists
    return DistanceMatrix(**data)

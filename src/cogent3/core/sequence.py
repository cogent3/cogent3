"""Contains classes that represent biological sequence data.

Notes
-----
These should be created via MolType.make_seq()
"""

from __future__ import annotations

import contextlib
import json
import re
import warnings
from collections import defaultdict
from functools import total_ordering
from operator import eq, ne
from random import shuffle
from typing import TYPE_CHECKING, Any, Self, SupportsIndex, cast

import numba
import numpy
import numpy.typing as npt

from cogent3._version import __version__
from cogent3.core import alphabet as c3_alphabet
from cogent3.core import genetic_code as c3_genetic_code
from cogent3.core import moltype as c3_moltype
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    AnnotatableMixin,
    FeatureDataType,
    SqliteAnnotationDbMixin,
    SupportsFeatures,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import (
    FeatureMap,
    IndelMap,
    LostSpan,
    Span,
    Strand,
    _LostSpan,
)
from cogent3.core.seqview import SeqView, SeqViewABC
from cogent3.core.slice_record import SliceRecord
from cogent3.draw.drawable import Shape
from cogent3.format.fasta import seqs_to_fasta
from cogent3.maths.stats.contingency import CategoryCounts, TestResult
from cogent3.maths.stats.number import CategoryCounter
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.dict_array import DictArray
from cogent3.util.misc import (
    DistanceFromMatrix,
    get_object_provenance,
    get_setting_from_environ,
    is_float,
    is_int,
)
from cogent3.util.transform import for_seq, per_shortest

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Mapping

    from cogent3.core.alignment import Aligned
    from cogent3.draw.drawable import Drawable, Shape

NumpyIntArrayType = npt.NDArray[numpy.integer]


# standard distance functions: left  because generally useful
frac_same = for_seq(f=eq, aggregator=sum, normalizer=per_shortest)
frac_diff = for_seq(f=ne, aggregator=sum, normalizer=per_shortest)


def _moltype_seq_from_rich_dict(
    data: dict[str, str | bytes | NumpyIntArrayType | dict[str, str]],
) -> tuple[c3_moltype.MolType[Any], str | bytes | NumpyIntArrayType]:
    """returns moltype and seq and mutates data so it can serve as kwargs to Sequence constructor"""
    data.pop("type")
    data.pop("version")
    data.pop("annotation_db", None)
    data.pop("annotation_offset", 0)

    moltype: c3_moltype.MolTypeLiteral | c3_moltype.MolType[Any] = cast(
        "c3_moltype.MolTypeLiteral", data.pop("moltype")
    )
    moltype = c3_moltype.get_moltype(moltype)

    seq = cast("str | bytes | NumpyIntArrayType", data.pop("seq"))
    return moltype, seq


@numba.jit(cache=True, nogil=True)
def count_kmers(
    seq: npt.NDArray[numpy.integer],
    num_states: int,
    k: int,
    dtype: npt.DTypeLike = numpy.uint64,
) -> NumpyIntArrayType:  # pragma: no cover
    """return freqs of valid k-mers using 1D indices

    Parameters
    ----------
    seq
        numpy array of uint8, assumed that canonical characters have
        indexes which are all < num_states
    num_states
        defines range of possible ints at a position
    k
        k-mer size
    dtype
        numpy dtype for the returned array
    """
    coeffs = c3_alphabet.coord_conversion_coeffs(num_states, k)
    kfreqs = numpy.zeros(num_states**k, dtype=dtype)
    if len(seq) < k:
        return kfreqs

    skip_until = 0
    for i in range(k):
        if seq[i] >= num_states:
            skip_until = i + 1

    idx = -1  # -1 means an invalid index
    biggest_coeff = coeffs[0]

    for i in range(len(seq) - k + 1):
        gained_char = seq[i + k - 1]
        if gained_char >= num_states:
            # we reset the kmer index to invalid
            # until we get a kmer with only canonical chars
            idx = -1
            skip_until = i + k

        if i < skip_until:
            continue

        if idx < 0:
            idx = (seq[i : i + k] * coeffs).sum()
            kfreqs[idx] += 1
            continue

        dropped_char = seq[i - 1]
        idx = (idx - dropped_char * biggest_coeff) * num_states + gained_char
        kfreqs[idx] += 1

    return kfreqs


# This is a design note about the annotation_offset attribute on sequences
#  and the related offset attribute on the sequence view.

# The SeqView object is responsible for performing slices and keeping track of
#  those slice operations relative to a parent. This is done via the start, stop
#  and step attributes.

# The annotation_offset attribute is intended to facilitate relating sequences
#  to annotations stored in coordinates of an annotated parent sequence. Consider
#  for example a gene that lies on a chromosome. If a user has an instance of that
#  gene's sequence in a sequence object, then the annotation_offset would be the
#  start position of that gene on the plus strand of the original chromosome.

# It's the case that an annotation_offset can only be specified by the user as
#  an argument to the sequence constructor. The design decision was made to
#  have the responsibility for providing the coordinates on the parent lie
#  with the SeqView class. These values are provided by the parent_start,
#  parent_stop and absolute_position methods.

# The state flow is in the constructor for the sequence class. There is
#  an assignment of the annotation_offset to the provided sequence_view.

# The only time the offset value is not at the SeqView level is the
#  sequence object is being serialised.


@total_ordering
class Sequence(AnnotatableMixin):
    """Holds the standard Sequence object. Immutable.

    Notes
    -----
    Sequences should be constructed by a MolType instance.
    """

    __slots__ = (
        "_annotation_db",
        "_repr_policy",
        "_seq",
        "info",
        "moltype",
        "name",
    )

    def __init__(
        self,
        moltype: c3_moltype.MolType[Any],
        seq: str | bytes | NumpyIntArrayType | SeqViewABC,
        *,
        name: str | None = None,
        info: dict[str, Any] | InfoClass | None = None,
        annotation_offset: int = 0,
        annotation_db: SupportsFeatures | None = None,
    ) -> None:
        """Initialize a sequence.

        Parameters
        ----------
        moltype
            MolType instance
        seq
            the raw sequence string, default is ''
        name
            the sequence name
        info
            Info object or dict
        annotation_offset
            integer indicating start position relative to annotations
        annotation_db
            optional annotation database

        Notes
        -----
        If a user attempts to add a feature and the annotation_db is None,
        a default BasicAnnotationDb instance will be created and used.
        """
        self.moltype = moltype
        self.name = name
        self._seq = _coerce_to_seqview(
            seq,
            name,
            self.moltype.most_degen_alphabet(),
            annotation_offset,
        )

        self.info = InfoClass(**(info or {}))
        self._repr_policy = {"num_pos": 60}
        self._annotation_db: list[SupportsFeatures] = self._init_annot_db_value(
            annotation_db
        )

    def __str__(self) -> str:
        result = numpy.array(self)
        return self.moltype.most_degen_alphabet().from_indices(result)

    def __bytes__(self) -> bytes:
        return str(self).encode("utf8")

    def __array__(
        self,
        dtype: numpy.dtype[numpy.integer] | None = None,
        copy: bool | None = None,
    ) -> NumpyIntArrayType:
        # using the array_value attribute means we can have
        # a aligned or seq data view here and the outcome will be the
        # same -- just the sequence is returned
        result = self._seq.array_value
        if self._seq.slice_record.is_reversed:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def to_array(self, apply_transforms: bool = True) -> NumpyIntArrayType:
        """returns the numpy array

        Parameters
        ----------
        apply_transforms
            if True, applies any reverse complement operation

        Notes
        -----
        Use this method with apply_transforms=False if you are
        creating data for storage in a SeqData instance.
        """
        if apply_transforms:
            return numpy.array(self)

        arr = self._seq.array_value
        if self._seq.is_reversed:
            # the reversal will have been applied in the SeqView
            # array_value method, so we undo that here.
            arr = arr[::-1]

        return arr

    def to_fasta(
        self,
        make_seqlabel: Callable[[Sequence], str] | None = None,
        block_size: int = 60,
    ) -> str:
        """Return string of self in FASTA format, no trailing newline

        Parameters
        ----------
        make_seqlabel
            callback function that takes the seq object and
            returns a label str

        """

        label = "0"

        if make_seqlabel is not None:
            label = make_seqlabel(self)
        elif hasattr(self, "label") and self.label:
            label = self.label
        elif hasattr(self, "name") and self.name:
            label = self.name
        return seqs_to_fasta({label: str(self)}, block_size=block_size)

    def to_phylip(self, name_len: int = 28, label_len: int = 30) -> str:
        """Return string of self in one line for PHYLIP, no newline.

        Default: max name length is 28, label length is 30.
        """
        return str(self.name)[:name_len].ljust(label_len) + str(self)

    def to_rich_dict(
        self,
        exclude_annotations: bool = True,
    ) -> dict[str, str | dict[str, str]]:
        """returns {'name': name, 'seq': sequence, 'moltype': moltype.label}

        Notes
        -----
        Deserialisation of the sequence object will not include the annotation_db
        even if exclude_annotations=False.
        """
        info: InfoClass | dict[str, Any] = {} if self.info is None else self.info
        if info.get("Refs", None) is not None and "Refs" in info:
            info.pop("Refs")

        data: dict[str, Any] = {
            "name": self.name,
            "seq": str(self),
            "moltype": self.moltype.label,
            "info": info or None,
            "type": get_object_provenance(self),
            "version": __version__,
        }
        if hasattr(self, "annotation_offset"):
            offset = int(self._seq.slice_record.parent_start)
            data |= {"annotation_offset": offset}

        if (
            hasattr(self, "annotation_db")
            and self.annotation_db
            and not exclude_annotations
        ):
            data["annotation_db"] = self.annotation_db.to_rich_dict()

        return data

    @classmethod
    def from_rich_dict(cls, data: dict[str, Any]) -> Sequence:
        """create a Sequence object from a rich dict"""
        if isinstance(data["seq"], dict):
            # this is a rich dict from the old type sequences
            data["seq"] = data.pop("seq")["init_args"]["seq"]

        moltype, seq = _moltype_seq_from_rich_dict(data)

        return cls(moltype=moltype, seq=seq, **data)

    def to_json(self) -> str:
        """returns a json formatted string"""
        return json.dumps(self.to_rich_dict())

    def count(self, item: str) -> int:
        """count() delegates to self._seq."""
        return str(self).count(item)

    def counts(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        warn: bool = False,
    ) -> CategoryCounter[str | bytes]:
        """returns dict of counts of motifs

        only non-overlapping motifs are counted.

        Parameters
        ----------
        motif_length
            number of elements per character.
        include_ambiguity
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gaps
            if True, motifs containing a gap character are included.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """
        data = numpy.array(self)
        if not len(data):
            return CategoryCounter()

        if motif_length != 1:
            if warn and len(data) % motif_length != 0:
                warnings.warn(
                    f"{self.name} length not divisible by {motif_length}, truncating",
                    stacklevel=2,
                )
            limit = (len(data) // motif_length) * motif_length
            data = data[:limit]
            data = data.reshape(-1, motif_length)

        unique_values, counts = numpy.unique(data, axis=0, return_counts=True)
        if motif_length == 1:
            unique_values = unique_values[:, None]

        indices = numpy.zeros(unique_values.shape[0], dtype=bool)
        if not allow_gap:
            indices |= numpy.apply_along_axis(
                self.moltype.is_gapped, axis=1, arr=unique_values
            )

        if not include_ambiguity:
            indices |= numpy.apply_along_axis(
                self.moltype.is_degenerate, axis=1, arr=unique_values
            )

        unique_values = unique_values[~indices]
        counts = counts[~indices]

        result: dict[str | bytes, int] = {}
        alpha = self.moltype.most_degen_alphabet()
        # The following approach is used because of bytes moltype.
        # We can't use the from_indices method on the bytes moltype
        # as that creates a string, but individual elements are
        # bytes, not all of which can be decoded. Thus we get an
        # empty instance and use join.
        monomer_type: str | bytes = type(alpha[0])()
        for motif, count in zip(unique_values, counts, strict=True):
            key = monomer_type.join([alpha[i] for i in motif])
            result[key] = int(count)
        return CategoryCounter(result)

    def count_ambiguous(self) -> int:
        """Returns the number of ambiguous characters in the sequence."""
        data = numpy.array(self)
        gap_index = cast("int", self.moltype.most_degen_alphabet().gap_index)
        return int(numpy.sum(data > gap_index))

    def count_kmers(
        self,
        k: int = 1,
        use_hook: str | None = None,
        **kwargs,  # noqa: ANN003
    ) -> NumpyIntArrayType:
        """return array of counts of all possible kmers of length k

        Parameters
        ----------
        k
            length of kmers to count
        use_hook
            name of a third-party package that implements the quick_tree
            hook. If not specified, defaults to the first available hook or
            the cogent3 quick_tree() app. To force default, set
            use_hook="cogent3".
        **kwargs
            additional arguments to pass to the hook app constructor

        Notes
        -----
        Only states in moltype.alphabet are allowed in a kmer (the canonical
        states). If using cogent3, to get the order of kmers as strings, use
        self.moltype.alphabet.get_kmer_alphabet(k). See the documentation
        for details of any third-party hooks.
        """
        from cogent3._plugin import get_count_kmers_hook

        kwargs = {"k": k, **kwargs}

        app = get_count_kmers_hook(name=use_hook, **kwargs)
        if app is not None:
            result = app([self])
            # make sure result from a sequence is 1D
            return result.flatten()

        seqarray = numpy.array(self)
        return count_kmers(seqarray, len(self.moltype.alphabet), k)

    def __lt__(self, other: Sequence) -> bool:
        """compares based on the sequence string."""
        return str(self) < str(other)

    def __eq__(self, other: object) -> bool:
        """compares based on the sequence string."""
        return str(self) == str(other)

    def __ne__(self, other: object) -> bool:
        """compares based on the sequence string."""
        return str(self) != str(other)

    def __hash__(self) -> int:
        """__hash__ behaves like the sequence string for dict lookup."""
        return hash(str(self))

    def __contains__(self, other: object) -> bool:
        """__contains__ checks whether other is in the sequence string."""
        if not isinstance(other, str):
            msg = f"Must use type str for __contains__, got {type(other)}."
            raise TypeError(msg)
        return other in str(self)

    def shuffle(self) -> Self:
        """returns a randomized copy of the Sequence object"""
        randomized_copy_list = list(self)
        shuffle(randomized_copy_list)
        return self.__class__(
            moltype=self.moltype,
            seq="".join(randomized_copy_list),
            info=self.info,
        )

    def strip_degenerate(self) -> Self:
        """Removes degenerate bases by stripping them out of the sequence."""
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.strip_degenerate(bytes(self)),
            info=self.info,
        )

    def strip_bad(self) -> Self:
        """Removes any symbols not in the alphabet."""
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.strip_bad(str(self)),
            info=self.info,
        )

    def strip_bad_and_gaps(self) -> Self:
        """Removes any symbols not in the alphabet, and any gaps. As the missing
        character could be a gap, this method will remove it as well."""
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.strip_bad_and_gaps(str(self)),
            info=self.info,
        )

    def is_gapped(self) -> bool:
        """Returns True if sequence contains gaps."""
        return self.moltype.is_gapped(str(self))

    def is_degenerate(self) -> bool:
        """Returns True if sequence contains degenerate characters."""
        return self.moltype.is_degenerate(str(self))

    def is_valid(self) -> bool:
        """Returns True if sequence contains no items absent from alphabet."""
        return self.moltype.is_valid(numpy.array(self))

    def is_strict(self) -> bool:
        """Returns True if sequence contains only monomers."""
        return self.moltype.alphabet.is_valid(numpy.array(self))

    def disambiguate(self, method: str = "strip") -> Self:
        """Returns a non-degenerate sequence from a degenerate one.

        Parameters
        ----------
        seq
            the sequence to be disambiguated
        method
            how to disambiguate the sequence, one of "strip", "random"
            strip: deletes any characters not in monomers or gaps
            random: assigns the possibilities at random, using equal frequencies
        """
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.disambiguate(str(self), method),
            info=self.info,
        )

    def degap(self) -> Self:
        """Deletes all gap characters from sequence."""
        result = self.__class__(
            moltype=self.moltype,
            seq=self.moltype.degap(bytes(self)),
            name=self.name,
            info=self.info,
        )
        result.annotation_db = self.annotation_db
        return result

    def gap_indices(self) -> NumpyIntArrayType:
        """Returns array of the indices of all gaps in the sequence"""
        return numpy.where(self.gap_vector())[0]

    def gap_vector(self) -> list[bool]:
        """Returns vector of True or False according to which pos are gaps or missing."""
        degen_gapped_alphabet = cast(
            "c3_alphabet.CharAlphabet[Any]", self.moltype.degen_gapped_alphabet
        )
        return (
            (numpy.array(self) == degen_gapped_alphabet.gap_index)
            | (numpy.array(self) == degen_gapped_alphabet.missing_index)
        ).tolist()

    def count_gaps(self) -> int:
        """Counts the gaps in the specified sequence."""
        return self.moltype.count_gaps(bytes(self))

    def count_degenerate(self) -> int:
        """Counts the degenerate bases in the specified sequence.

        Notes
        -----
        gap and missing characters are counted as degenerate.
        """
        # design: refactor
        # should gap and missing characters be counted as degenerate?
        return self.moltype.count_degenerate(bytes(self))

    def count_variants(self) -> int:
        """Counts number of possible sequences matching the sequence, given
        any ambiguous characters in the sequence.

        Notes
        -----
        Uses self.ambiguitues to decide how many possibilities there are at
        each position in the sequence and calculates the permutations.
        """
        return self.moltype.count_variants(str(self))

    def mw(self, method: str = "random", delta: float | None = None) -> float:
        """Returns the molecular weight of (one strand of) the sequence.

        Parameters
        ----------
        method
            If the sequence is ambiguous, uses method (random or strip) to
            disambiguate the sequence.
        delta
            If delta is passed in, adds delta per strand. Default is None, which
            uses the alphabet default. Typically, this adds 18 Da for terminal
            water. However, note that the default nucleic acid weight assumes
            5' monophosphate and 3' OH: pass in delta=18.0 if you want 5' OH as
            well.

        Notes
        -----
        this method only calculates the MW of the coding strand. If you want
        the MW of the reverse strand, add self.rc().mw(). DO NOT just multiply
        the MW by 2: the results may not be accurate due to strand bias, e.g.
        in mitochondrial genomes.
        """
        return self.moltype.mw(str(self), method, delta)

    def can_match(self, other: Self) -> bool:
        """Returns True if every pos in self could match same pos in other.

        Truncates at length of shorter sequence.
        gaps are only allowed to match other gaps.
        """
        return self.moltype.can_match(str(self), str(other))

    def diff(self, other: Self) -> float:
        """Returns number of differences between self and other.

        Notes
        -----
        Truncates at the length of the shorter sequence.
        """
        return self.distance(other)

    def distance(
        self,
        other: Self,
        function: Callable[[str, str], float] | None = None,
    ) -> float:
        """Returns distance between self and other using function(i,j).

        Parameters
        ----------
        other
            a sequence to compare to self
        function
            takes two seq residues and returns a number. To turn a 2D matrix into
            a function, use cogent3.util.miscs.DistanceFromMatrix(matrix).

        Notes
        -----
        Truncates at the length of the shorter sequence.

        The function acts on two _elements_ of the sequences, not the two
        sequences themselves (i.e. the behavior will be the same for every
        position in the sequences, such as identity scoring or a function
        derived from a distance matrix as suggested above). One limitation of
        this approach is that the distance function cannot use properties of
        the sequences themselves: for example, it cannot use the lengths of the
        sequences to normalize the scores as percent similarities or percent
        differences.

        If you want functions that act on the two sequences themselves, there
        is no particular advantage in making these functions methods of the
        first sequences by passing them in as parameters like the function
        in this method. It makes more sense to use them as standalone functions.
        The factory function cogent3.util.transform.for_seq is useful for
        converting per-element functions into per-sequence functions, since it
        takes as parameters a per-element scoring function, a score aggregation
        function, and a normalization function (which itself takes the two
        sequences as parameters), returning a single function that combines
        these functions and that acts on two complete sequences.
        """
        if function is None:
            # use identity scoring function
            def function(a: str, b: str) -> bool:
                return a != b

        distance: float = 0
        for first, second in zip(self, other, strict=False):
            distance += function(first, second)
        return distance

    def matrix_distance(
        self,
        other: Self,
        matrix: Mapping[str, Mapping[str, float]],
    ) -> float:
        """Returns distance between self and other using a score matrix.

        Warnings
        --------
        The matrix must explicitly contain scores for the case where
        a position is the same in self and other (e.g. for a distance matrix,
        an identity between U and U might have a score of 0). The reason the
        scores for the 'diagonals' need to be passed explicitly is that for
        some kinds of distance matrices, e.g. log-odds matrices, the 'diagonal'
        scores differ from each other. If these elements are missing, this
        function will raise a KeyError at the first position that the two
        sequences are identical.
        """
        return self.distance(other, DistanceFromMatrix(matrix))

    def frac_same(self, other: Self) -> float:
        """Returns fraction of positions where self and other are the same.

        Notes
        -----
        Truncates at length of shorter sequence. Will return  0 if one sequence
        is empty.
        """
        return frac_same(self, other)

    def frac_diff(self, other: Self) -> float:
        """Returns fraction of positions where self and other differ.

        Notes
        -----
        Truncates at length of shorter sequence. Will return  0 if one sequence
        is empty.
        """
        return frac_diff(self, other)

    def frac_same_gaps(self, other: Self) -> float:
        """Returns fraction of positions where self and other share gap states.

        In other words, if self and other are both all gaps, or both all
        non-gaps, or both have gaps in the same places, frac_same_gaps will
        return 1.0. If self is all gaps and other has no gaps, frac_same_gaps
        will return 0.0. Returns 0 if one sequence is empty.

        Uses self's gap characters for both sequences.
        """
        if not self or not other:
            return 0.0

        is_gap = self.moltype.gaps.__contains__
        return sum(
            [is_gap(i) == is_gap(j) for i, j in zip(self, other, strict=False)],
        ) / min(
            len(self),
            len(other),
        )

    def frac_diff_gaps(self, other: Self) -> float:
        """Returns frac. of positions where self and other's gap states differ.

        In other words, if self and other are both all gaps, or both all
        non-gaps, or both have gaps in the same places, frac_diff_gaps will
        return 0.0. If self is all gaps and other has no gaps, frac_diff_gaps
        will return 1.0.

        Returns 0 if one sequence is empty.

        Uses self's gap characters for both sequences.
        """
        if not self or not other:
            return 0.0
        return 1.0 - self.frac_same_gaps(other)

    def frac_same_non_gaps(self, other: Self) -> float:
        """Returns fraction of non-gap positions where self matches other.

        Doesn't count any position where self or other has a gap.
        Truncates at the length of the shorter sequence.

        Returns 0 if one sequence is empty.
        """
        # refactor: simplify
        # refactor: array - make use of self._seq.array_value
        if not self or not other:
            return 0.0

        is_gap = self.moltype.gaps.__contains__
        count = 0
        identities = 0
        for i, j in zip(self, other, strict=False):
            if is_gap(i) or is_gap(j):
                continue
            count += 1
            if i == j:
                identities += 1

        if count:
            return identities / count
        # there were no positions that weren't gaps
        return 0

    def frac_diff_non_gaps(self, other: Self) -> float:
        """Returns fraction of non-gap positions where self differs from other.

        Doesn't count any position where self or other has a gap.
        Truncates at the length of the shorter sequence.

        Returns 0 if one sequence is empty. Note that this means that
        frac_diff_non_gaps is _not_ the same as 1 - frac_same_non_gaps, since both
        return 0 if one sequence is empty.
        """
        if not self or not other:
            return 0.0

        is_gap = self.moltype.gaps.__contains__
        count = 0
        diffs = 0
        for i, j in zip(self, other, strict=False):
            if is_gap(i) or is_gap(j):
                continue
            count += 1
            if i != j:
                diffs += 1

        if count:
            return diffs / count
        # there were no positions that weren't gaps
        return 0

    def frac_similar(
        self,
        other: Self,
        similar_pairs: dict[tuple[str, str], Any],
    ) -> float:
        """Returns fraction of positions where self[i] is similar to other[i].

        similar_pairs must be a dict such that d[(i,j)] exists if i and j are
        to be counted as similar. Use PairsFromGroups in cogent3.util.misc to
        construct such a dict from a list of lists of similar residues.

        Truncates at the length of the shorter sequence.

        Note: current implementation re-creates the distance function each
        time, so may be expensive compared to creating the distance function
        using for_seq separately.

        Returns 0 if one sequence is empty.
        """
        if not self or not other:
            return 0.0

        def is_in_similar_pairs(x: str, y: str) -> bool:
            return (x, y) in similar_pairs

        return for_seq(f=is_in_similar_pairs, normalizer=per_shortest)(
            self,
            other,
        )

    def with_termini_unknown(self) -> Self:
        """Returns copy of sequence with terminal gaps remapped as missing."""
        gaps = self.gap_vector()
        first_nongap = last_nongap = None
        for i, state in enumerate(gaps):
            if not state:
                if first_nongap is None:
                    first_nongap = i
                last_nongap = i
        missing = cast("str", self.moltype.missing)
        if first_nongap is None:
            return self.__class__(
                moltype=self.moltype,
                seq=missing * len(self),
                info=self.info,
            )
        last_nongap = cast("int", last_nongap)

        prefix = missing * first_nongap
        mid = str(self)[first_nongap : last_nongap + 1]
        suffix = missing * (len(self) - last_nongap - 1)
        return self.__class__(
            moltype=self.moltype,
            seq=prefix + mid + suffix,
            info=self.info,
        )

    def _repr_html_(self) -> str:
        settings = self._repr_policy.copy()
        env_vals = get_setting_from_environ(
            "COGENT3_ALIGNMENT_REPR_POLICY",
            {"num_pos": int},
        )
        settings.update(env_vals)
        return self.to_html(limit=settings["num_pos"])

    def to_html(
        self,
        wrap: int = 60,
        limit: int | None = None,
        colors: Mapping[str, str] | None = None,
        font_size: int = 12,
        font_family: str = "Lucida Console",
    ) -> str:
        """returns html with embedded styles for sequence colouring

        Parameters
        ----------
        wrap
            maximum number of printed bases, defaults to
            alignment length
        limit
            truncate alignment to this length
        colors
            dict of {char: color} to use for coloring
        font_size
            in points. Affects labels and sequence and line spacing
            (proportional to value)
        font_family
            string denoting font family

        To display in jupyter notebook:

            >>> from IPython.core.display import HTML
            >>> HTML(aln.to_html())
        """
        css, styles = self.moltype.get_css_style(
            colors=colors,
            font_size=font_size,
            font_family=font_family,
        )

        seq = str(self)
        seq = seq if limit is None else seq[:limit]
        seqlen = len(seq)
        if gaps := self.moltype.gaps:
            non_hyphen = "".join(gaps - {"-"})
            # make sure hyphen at end of negated character group
            # so it's not interpreted as a character range
            chars = f"{non_hyphen}-" if "-" in gaps else non_hyphen
            # looking for gaps at the the start of the seq
            start_gap = re.search(f"^[{chars}]+", "".join(seq))
            # looking for gaps at the end of the seq
            end_gap = re.search(f"[{gaps}]+$", "".join(seq))

            start = 0 if start_gap is None else start_gap.end()
            end = len(seq) if end_gap is None else end_gap.start()
        else:
            start = 0
            end = len(seq)

        seq_style: list[str] = []
        template = '<span class="%s">%%s</span>'
        styled_seq: list[str] = []
        for i in range(seqlen):
            char = seq[i]
            if i < start or i >= end:
                style = f"terminal_ambig_{self.moltype.label}"
            else:
                style = styles[char]

            seq_style.append(template % style)
            styled_seq.append(seq_style[-1] % char)

        # make a html table
        seq_array = numpy.array(styled_seq, dtype="O")
        table = ["<table>"]
        seq_ = "<td>%s</td>"
        label_ = '<td class="label">%s</td>'
        num_row_ = '<tr class="num_row"><td></td><td><b>{:,d}</b></td></tr>'
        for i in range(0, seqlen, wrap):
            table.append(num_row_.format(i))
            seqblock = seq_array[i : i + wrap].tolist()
            seqblock = "".join(seqblock)
            row = "".join([label_ % self.name, seq_ % seqblock])
            table.append(f"<tr>{row}</tr>")
        table.append("</table>")
        class_name = self.__class__.__name__
        if limit and limit < len(self):
            summary = f"{class_name}, length={len(self):,} (truncated to {limit if limit else len(self)})"
        else:
            summary = f"{class_name}, length={len(self):,}"

        text = [
            "<style>",
            ".c3seq table {margin: 10px 0;}",
            ".c3seq td { border: none !important; text-align: left !important; }",
            ".c3seq tr:not(.num_row) td span {margin: 0 2px;}",
            ".c3seq tr:nth-child(even) {background: #f7f7f7;}",
            ".c3seq .num_row {background-color:rgba(161, 195, 209, 0.5) !important; border-top: solid 1px black; }",
            f".c3seq .label {{ font-size: {font_size}pt ; text-align: right !important; "
            "color: black !important; padding: 0 4px; }}",
            "\n".join([f".c3seq {style}" for style in css]),
            "</style>",
            '<div class="c3seq">',
            "\n".join(table),
            f"<p><i>{summary}</i></p>",
            "</div>",
        ]
        return "\n".join(text)

    def __add__(self, other: Self) -> Self:
        """Adds two sequences (other can be a string as well)."""
        if hasattr(other, "moltype") and self.moltype != other.moltype:
            msg = f"MolTypes don't match: ({self.moltype},{other.moltype})"
            raise ValueError(
                msg,
            )
        other_seq = str(other)

        if not self.moltype.most_degen_alphabet().is_valid(str(other)):
            msg = (
                f"Invalid sequence characters in other for moltype={self.moltype.label}"
            )
            raise c3_alphabet.AlphabetError(
                msg,
            )

        # If two sequences with the same name are being added together the name should not be None
        if type(other) is type(self):
            name = self.name if self.name == other.name else None
        else:
            name = None

        return self.__class__(
            moltype=self.moltype,
            seq=str(self) + other_seq,
            name=name,
        )

    @property
    def annotation_offset(self) -> int:
        """
        The offset between annotation coordinates and sequence coordinates.

        The offset can be used to adjust annotation coordinates to match the position
        of the given Sequence within a larger genomic context. For example, if the
        Annotations are with respect to a chromosome, and the sequence represents
        a gene that is 100 bp from the start of a chromosome, the offset can be set to
        100 to ensure that the gene's annotations are aligned with the appropriate
        genomic positions.


        Returns:
            int: The offset between annotation coordinates and sequence coordinates.
        """
        parent_offset = self._seq.parent_offset
        slice_start = self._seq.slice_record.parent_start
        return parent_offset + slice_start

    def get_features(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        allow_partial: bool = False,
        **kwargs: Any,
    ) -> Iterator[Feature[Sequence]]:
        """yields Feature instances

        Parameters
        ----------
        biotype
            biotype of the feature
        name
            name of the feature
        start, stop
            start, stop positions to search between, relative to offset
            of this sequence. If not provided, entire span of sequence is used.
        kwargs
            keyword arguments passed to annotation_db.get_features_matching()

        Notes
        -----
        When dealing with a nucleic acid moltype, the returned features will
        yield a sequence segment that is consistently oriented irrespective
        of strand of the current instance.
        """

        if self._annotation_db is None:
            return None

        start = start or 0
        stop = stop or len(self)

        # handle negative start / stop
        start = start + len(self) if start < 0 else start
        stop = stop + len(self) if stop < 0 else stop

        start, stop = (start, stop) if start < stop else (stop, start)
        stop = min(
            len(self), stop
        )  # otherwise index error from absolute_position method

        # we set include_boundary=False because start is inclusive indexing,
        # i,e., the start cannot be equal to the length of the view
        (
            # parent_id can differ from self.name if there was a
            # rename operation
            parent_id,
            *_,
            _,
        ) = self.parent_coordinates()
        parent_offset = self._seq.parent_offset
        sr = self._seq.slice_record
        query_start = (
            sr.absolute_position(
                start,
                include_boundary=False,
            )
            + parent_offset
        )
        # we set include_boundary=True because stop is exclusive indexing,
        # i,e., the stop can be equal to the length of the view
        query_stop = (
            sr.absolute_position(
                stop,
                include_boundary=True,
            )
            + parent_offset
        )

        query_start, query_stop = (
            min(query_stop, query_start),
            max(query_stop, query_start),
        )
        query_start = max(query_start, 0)
        # in the underlying db, features are always plus strand oriented
        # (need to check that for user defined features)
        # a '-' strand feature in the db may be [(2, 4), (6, 8)]
        # so we would take 7-6, 3-2 and complement it
        # if the current sequence is reverse complemented
        # then a '+' feature from the db needs to be nucleic reversed
        # and a '-' feature from the db needs to be nucleic reversed too

        # making this more complicated is self.make_feature() assumes
        # all coordinates are with respect to the current orientation
        # so self.make_feature(data(spans=[(2,4)])) has different meaning if
        # self is rc'ed, it would correspond to len(self)-2, etc...
        # To piggy-back on that method we need to convert our feature spans
        # into the current orientation. HOWEVER, we also have the reversed
        # flag which comes back from the db
        kwargs |= {"allow_partial": allow_partial}
        for feature in self.annotation_db.get_features_matching(
            seqid=parent_id,
            name=name,
            biotype=biotype,
            start=query_start,
            stop=query_stop,
            **kwargs,
        ):
            # spans need to be converted from absolute to relative positions
            # DO NOT do adjustment in make_feature since that's user facing,
            # and we expect them to make a feature manually wrt to their
            # current view
            spans = numpy.array(feature["spans"], dtype=int)
            for i, v in enumerate(spans.ravel()):
                rel_pos = sr.relative_position(v) - parent_offset
                spans.ravel()[i] = rel_pos

            if sr.is_reversed:
                # see above comment
                spans = len(self) - spans

            feature["spans"] = spans.tolist()
            yield self.make_feature(feature)

    def make_feature(self, feature: FeatureDataType, *args: Any) -> Feature[Sequence]:
        """
        return an Feature instance from feature data

        Parameters
        ----------
        feature
            dict of key data to make an Feature instance

        Notes
        -----
        Unlike add_feature(), this method does not add the feature to the
        database.
        We assume that spans represent the coordinates for this instance!
        """
        feature_dict: dict[str, Any] = dict(feature)
        seq_rced = self._seq.is_reversed
        spans = feature_dict.pop("spans", None)
        revd = Strand.from_value(feature_dict.pop("strand", None)) is Strand.MINUS
        feature_dict["strand"] = (
            Strand.PLUS.value if revd == seq_rced else Strand.MINUS.value
        )

        vals = numpy.array(spans)
        pre = abs(vals.min()) if vals.min() < 0 else 0
        post = abs(vals.max() - len(self)) if vals.max() > len(self) else 0

        # we find the spans > len(self)
        new_spans: list[tuple[int, int]] = []
        for coord in vals:
            new = coord[:]
            if coord.min() < 0 < coord.max():
                new = coord[:]
                new[new < 0] = 0
            elif coord.min() < len(self) < coord.max():
                new[new > len(self)] = len(self)
            elif coord[0] == coord[1] or coord.min() > len(self) or coord.max() < 0:
                # would really to adjust pre and post to be up to the next span
                continue
            new_spans.append(new.tolist())

        fmap = FeatureMap.from_locations(locations=new_spans, parent_length=len(self))
        if pre or post:
            # create a lost span to represent the segment missing from
            # the instance
            spans = list(fmap.iter_spans())
            if pre:
                spans.insert(0, LostSpan(pre))
            if post:
                spans.append(LostSpan(post))
            fmap = FeatureMap(spans=spans, parent_length=len(self))

        if seq_rced:
            fmap = fmap.nucleic_reversed()

        feature_dict.pop("on_alignment", None)
        feature_dict.pop("seqid", None)
        return Feature(
            parent=self, seqid=cast("str", self.name), map=fmap, **feature_dict
        )

    def add_feature(
        self,
        *,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]],
        parent_id: str | None = None,
        strand: str | None = None,
        on_alignment: bool = False,
        seqid: str | None = None,
    ) -> Feature[Sequence]:
        """
        add a feature to annotation_db

        Parameters
        ----------
        biotype
            biological type
        name
            name of the feature
        spans
            coordinates for this sequence
        parent_id
            name of the feature parent
        strand
            '+' or '-', defaults to '+'
        on_alignment
            whether the feature spans are alignment coordinates
        seqid
            ignored since the feature is added to this sequence

        Returns
        -------
        Feature instance
        """
        local_vars = locals()
        local_vars["strand"] = Strand.from_value(strand).value

        feature_data = FeatureDataType(
            seqid=cast("str", self.name),
            **{n: v for n, v in local_vars.items() if n not in ("self", "seqid")},  # type: ignore[typeddict-item]
        )

        self.annotation_db.add_feature(**cast("dict[str, Any]", feature_data))
        for discard in ("on_alignment", "parent_id"):
            feature_data.pop(discard)
        return self.make_feature(feature_data)

    def to_moltype(
        self, moltype: c3_moltype.MolTypeLiteral | c3_moltype.MolType[Any]
    ) -> Sequence:
        """returns copy of self with moltype seq

        Parameters
        ----------
        moltype
            molecular type

        Notes
        -----
        This method cannot convert between nucleic acids and proteins. Use
        get_translation() for that.

        When applied to a sequence in a SequenceCollection, the resulting
        sequence will no longer be part of the collection.
        """

        if not moltype:
            msg = f"unknown moltype '{moltype}'"
            raise ValueError(msg)

        moltype = c3_moltype.get_moltype(moltype)

        if moltype is self.moltype:
            return self

        seq = numpy.array(self)
        self_alpha = self.moltype.most_degen_alphabet()
        other_alpha = moltype.most_degen_alphabet()
        if len(self.moltype.alphabet) != len(moltype.alphabet):
            # use alphabet converter
            seq = self_alpha.convert_seq_array_to(
                seq=seq,
                alphabet=other_alpha,
                check_valid=False,
            )

        if not other_alpha.is_valid(seq):
            msg = (
                f"Changing from old moltype={self.moltype.label!r} to new "
                f"moltype={moltype.label!r} is not valid for this data"
            )
            raise c3_moltype.MolTypeError(
                msg,
            )
        sv = SeqView(
            parent=seq,
            parent_len=len(seq),
            alphabet=moltype.most_degen_alphabet(),
        )
        return moltype.make_sequence(
            seq=sv,
            name=self.name,
            check_seq=False,
            info=self.info,
            annotation_db=self.annotation_db,
        )

    def copy_annotations(self, seq_db: SqliteAnnotationDbMixin) -> None:
        """copy annotations into attached annotation db

        Parameters
        ----------
        seq_db
            compatible annotation db

        Notes
        -----
        Only copies annotations for records with seqid equal to self.name
        """
        if not isinstance(seq_db, SupportsFeatures):
            msg = f"type {type(seq_db)} does not match SupportsFeatures interface"
            raise TypeError(
                msg,
            )

        if not seq_db.num_matches(seqid=self.name):
            return

        if self._annotation_db and not self.annotation_db.compatible(seq_db):
            msg = f"type {type(seq_db)} != {type(self.annotation_db)}"
            raise TypeError(msg)

        self.annotation_db.update(seq_db, seqids=self.name)

    def copy(
        self,
        exclude_annotations: bool = False,
        sliced: bool = True,
    ) -> Self:
        """returns a copy of self

        Parameters
        ----------
        sliced
            Slices underlying sequence with start/end of self coordinates. The offset
            is retained.
        exclude_annotations
            drops annotation_db when True
        """
        # when slicing a seqview with sliced=True, seqview discards the
        # offset attribute. Sequence then needs to be provided to its
        # constructor from its current state.
        offset = self.annotation_offset if sliced else 0
        data = self._seq.copy(sliced=sliced)
        return self.__class__(
            moltype=self.moltype,
            seq=data,
            name=self.name,
            info=self.info,
            annotation_offset=offset,
            annotation_db=None if exclude_annotations else self.annotation_db,
        )

    def with_masked_annotations(
        self,
        biotypes: str | Iterable[str],
        mask_char: str | None = None,
        shadow: bool = False,
    ) -> Self:
        """returns a sequence with annot_types regions replaced by mask_char
        if shadow is False, otherwise all other regions are masked.

        Parameters
        ----------
        annot_types
            annotation type(s)
        mask_char
            must be a character valid for the seq MolType. The
            default value is the most ambiguous character, eg. '?' for DNA
        shadow
            whether to mask the annotated regions, or everything but
            the annotated regions
        """
        if mask_char is None:
            assert self.moltype.ambiguities is not None
            mask_char = (
                self.moltype.missing
                or max(self.moltype.ambiguities.items(), key=lambda x: len(x[1]))[0]
            )
        assert mask_char in self.moltype.most_degen_alphabet(), (
            f"Invalid mask_char {mask_char}"
        )

        annotations: list[Feature[Sequence]] = []
        biotypes = [biotypes] if isinstance(biotypes, str) else biotypes
        for annot_type in biotypes:
            annotations += list(
                self.get_features(biotype=annot_type, allow_partial=True),
            )
        if not annotations:
            return self

        region = annotations[0].union(annotations[1:])

        if shadow:
            region = region.shadow()

        i = 0
        segments: list[str] = []
        coords = region.map.get_coordinates()
        for b, e in coords:
            segments.extend((str(self[i:b]), mask_char * (e - b)))
            i = e
        segments.append(str(self[i:]))

        new = self.__class__(
            moltype=self.moltype,
            seq="".join(segments),
            name=self.name,
            info=self.info,
        )
        new.annotation_db = self.annotation_db
        return new

    def gapped_by_map_segment_iter(
        self,
        segment_map: IndelMap | FeatureMap,
        allow_gaps: bool = True,
        recode_gaps: bool = False,
    ) -> Iterator[str]:
        if not allow_gaps and not segment_map.complete:
            msg = f"gap(s) in map {segment_map}"
            raise ValueError(msg)

        for span in segment_map.iter_spans():
            if span.lost:
                span = cast("_LostSpan", span)
                unknown = "?" if span.terminal or recode_gaps else "-"
                seg = unknown * span.length
            else:
                span = cast("Span", span)
                seg = str(self[span.start : span.end])

            yield seg

    def gapped_by_map_motif_iter(
        self,
        segment_map: IndelMap,
    ) -> Iterator[str]:
        for segment in self.gapped_by_map_segment_iter(segment_map):
            yield from segment

    def gapped_by_map(
        self,
        segment_map: IndelMap | FeatureMap,
        recode_gaps: bool = False,
    ) -> Self:
        segments = self.gapped_by_map_segment_iter(segment_map, True, recode_gaps)
        return self.__class__(
            moltype=self.moltype,
            seq="".join(segments),
            name=self.name,
            info=self.info,
        )

    def _mapped(self, segment_map: FeatureMap) -> Self:
        # Called by generic __getitem__
        seq: SeqViewABC | str
        if segment_map.num_spans == 1:
            seq = self._seq[segment_map.start : segment_map.end]
        else:
            segments = self.gapped_by_map_segment_iter(segment_map, allow_gaps=False)
            seq = "".join(segments)

        return self.__class__(
            moltype=self.moltype,
            seq=seq,
            name=self.name,
            info=self.info,
        )

    def __repr__(self) -> str:
        myclass = f"{self.__class__.__name__}"
        myclass = myclass.split(".")[-1]
        seq = f"{str(self)[:7]}... {len(self):,}" if len(self) > 10 else str(self)
        return f"{myclass}({seq})"

    def __getitem__(
        self, index: Feature[Sequence] | FeatureMap | slice | SupportsIndex
    ) -> Self:
        preserve_offset = False
        if isinstance(index, Feature):
            if index.parent is not self:
                msg = "cannot slice Feature not bound to self"
                raise ValueError(msg)
            return cast("Self", index.get_slice())

        if isinstance(index, FeatureMap):
            new = self._mapped(index)
            # annotations have no meaning if disjoint slicing segments
            preserve_offset = index.num_spans == 1
        elif isinstance(index, slice) or is_int(index):
            new = self.__class__(
                moltype=self.moltype,
                seq=self._seq[index],
                name=self.name,
                info=self.info,
            )
            stride = getattr(index, "step", 1) or 1
            preserve_offset = stride > 0
        else:
            msg = f"Cannot slice using a {type(index)}"
            raise TypeError(msg)

        if isinstance(index, list | tuple):
            msg = "cannot slice using list or tuple"
            raise TypeError(msg)

        if self._annotation_db and preserve_offset:
            new.replace_annotation_db(self.annotation_db, check=False)

        if is_float(index):
            msg = "cannot slice using float"
            raise TypeError(msg)

        if hasattr(self, "_repr_policy"):
            new._repr_policy.update(self._repr_policy)

        return new

    def __iter__(self) -> Iterator[str]:
        yield from iter(str(self))

    def get_name(self) -> str | None:
        """Return the sequence name -- should just use name instead."""
        return self.name

    def __len__(self) -> int:
        return len(self._seq)

    def get_type(self) -> str:
        """Return the sequence type as moltype label."""
        return self.moltype.label

    def resolved_ambiguities(self) -> list[set[str]]:
        """Returns a list of sets of strings."""
        ambigs = cast("dict[str, frozenset[str]]", self.moltype.ambiguities)
        return [set(ambigs.get(motif, motif)) for motif in str(self)]

    def iter_kmers(self, k: int, strict: bool = True) -> Iterator[str]:
        """generates all overlapping k-mers.
        When strict is True, the characters in the k-mer must be
        a subset of the canonical characters for the moltype"""
        if k <= 0:
            msg = f"k must be an int > 0, not {k}"
            raise ValueError(msg)

        if not isinstance(k, int):
            msg = f"k must be an int, not {k}"
            raise ValueError(msg)

        md = self.moltype.most_degen_alphabet()
        gap_index = md.gap_index
        arr_to_str = md.from_indices
        if gap_index is None:
            gap_index = len(self.moltype)
        seq = numpy.array(self)
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if not strict or kmer.max() < gap_index:
                yield arr_to_str(kmer)

    def get_kmers(self, k: int, strict: bool = True) -> list[str]:
        """return all overlapping k-mers"""
        return list(self.iter_kmers(k, strict))

    def sliding_windows(
        self,
        window: int,
        step: int,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Self]:
        """Generator function that yield new sequence objects
        of a given length at a given interval.

        Parameters
        ----------
        window
            The length of the returned sequence
        step
            The interval between the start of the returned
            sequence objects
        start
            first window start position
        end
            last window start position

        """
        if start is None:
            start = 0
        if end is None:
            end = len(self) - window + 1
        end = min(len(self) - window + 1, end)
        if start < end and len(self) - end >= window - 1:
            for pos in range(start, end, step):
                yield self[pos : pos + window]

    def get_in_motif_size(
        self, motif_length: int = 1, warn: bool = False
    ) -> list[str] | str:
        """returns sequence as list of non-overlapping motifs

        Parameters
        ----------
        motif_length
            length of the motifs
        warn
            whether to notify of an incomplete terminal motif
        """
        seq: SeqViewABC | str = self._seq
        if isinstance(seq, SeqViewABC):
            seq = str(self)
        if motif_length == 1:
            return seq

        length = len(seq)
        remainder = length % motif_length
        if remainder and warn:
            warnings.warn(
                f'Dropped remainder "{seq[-remainder:]}" from end of sequence',
                stacklevel=2,
            )
        return [
            seq[i : i + motif_length]
            for i in range(0, length - remainder, motif_length)
        ]

    def parse_out_gaps(self) -> tuple[IndelMap, Self]:
        """returns Map corresponding to gap locations and ungapped Sequence"""
        gap = re.compile(f"[{re.escape(cast('str', self.moltype.gap))}]+")
        seq_str = str(self)
        gap_pos_list: list[int] = []
        cum_lengths_list: list[int] = []
        for match in gap.finditer(seq_str):
            pos = match.start()
            gap_pos_list.append(pos)
            cum_lengths_list.append(match.end() - pos)

        gap_pos = numpy.array(gap_pos_list)
        cum_lengths = numpy.array(cum_lengths_list).cumsum()
        gap_pos[1:] = gap_pos[1:] - cum_lengths[:-1]

        seq = self.__class__(
            moltype=self.moltype,
            seq=gap.sub("", seq_str),
            name=self.get_name(),
            info=self.info,
        )
        indel_map = IndelMap(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_lengths,
            parent_length=len(seq),
        )
        seq.annotation_db = self.annotation_db
        return indel_map, seq

    def is_annotated(
        self,
        biotype: str | tuple[str] | None = None,
    ) -> bool:
        """returns True if sequence parent name has any annotations

        Parameters
        ----------
        biotype
            amend condition to return True only if the sequence is
            annotated with one of provided biotypes.
        """
        if not self._annotation_db:
            return False
        with contextlib.suppress(AttributeError):
            return (
                self.annotation_db.num_matches(seqid=self._seq.seqid, biotype=biotype)
                != 0
            )
        return False

    def annotate_matches_to(
        self,
        pattern: str,
        biotype: str,
        name: str,
        allow_multiple: bool = False,
    ) -> list[Feature[Sequence]]:
        """Adds an annotation at sequence positions matching pattern.

        Parameters
        ----------
        pattern
            The search string for which annotations are made. IUPAC ambiguities
            are converted to regex on sequences with the appropriate MolType.
        biotype
            The type of the annotation (e.g. "domain").
        name
            The name of the annotation.
        allow_multiple
            If True, allows multiple occurrences of the input pattern. Otherwise,
            only the first match is used.

        Returns
        -------
        Returns a list of Feature instances.
        """
        with contextlib.suppress(ValueError):
            # assume ValueError is due to already being a regex
            pattern = self.moltype.to_regex(seq=pattern)

        pos = [m.span() for m in re.finditer(pattern, str(self))]
        if not pos:
            return []

        num_match = len(pos) if allow_multiple else 1
        return [
            self.add_feature(
                biotype=biotype,
                name=f"{name}:{i}" if allow_multiple else name,
                spans=[pos[i]],
            )
            for i in range(num_match)
        ]

    def get_drawables(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
    ) -> dict[str, list[Shape]]:
        """returns a dict of drawables, keyed by type

        Parameters
        ----------
        biotype
            passed to get_features(biotype). Can be a single biotype or
            series. Only features matching this will be included.
        """
        # make sure the drawables are unique by adding to a set
        features = set(self.get_features(biotype=biotype, allow_partial=True))
        result: dict[str, list[Shape]] = defaultdict(list)
        for f in features:
            result[f.biotype].append(f.get_drawable())
        return result

    def get_drawable(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        width: float = 600,
        vertical: bool = False,
    ) -> Drawable | None:
        """make a figure from sequence features

        Parameters
        ----------
        biotype
            passed to get_features(biotype). Can be a single biotype or
            series. Only features matching this will be included.
        width
            width in pixels
        vertical
            rotates the drawable

        Returns
        -------
        a Drawable instance

        Notes
        -----
        If provided, the biotype is used for plot order.
        """
        from cogent3.draw.drawable import Drawable

        drawables = self.get_drawables(biotype=biotype)
        if not drawables:
            return None

        biotype = list(drawables) if biotype is None else biotype
        biotypes = (biotype,) if isinstance(biotype, str) else biotype

        # we order by tracks
        top: float = 0
        space = 0.25
        annotes: list[Shape] = []
        annott = None
        for feature_type in biotypes:
            new_bottom = top + space
            for i, annott in enumerate(drawables[feature_type]):
                annott.shift(y=new_bottom - annott.bottom)
                if i > 0:
                    annott._showlegend = False
                annotes.append(annott)

            top = cast("Shape", annott).top

        top += space
        height = max((top / len(self)) * width, 300)
        xaxis: dict[str, Any] = {
            "range": [0, len(self)],
            "zeroline": False,
            "showline": True,
        }
        yaxis: dict[str, Any] = {
            "range": [0, top],
            "visible": False,
            "zeroline": True,
            "showline": True,
        }

        if vertical:
            all_traces = [t.T.as_trace() for t in annotes]
            width, height = height, width
            xaxis, yaxis = yaxis, xaxis
        else:
            all_traces = [t.as_trace() for t in annotes]

        drawer = Drawable(
            title=self.name,
            traces=all_traces,
            width=width,
            height=height,
        )
        drawer.layout.update(xaxis=xaxis, yaxis=yaxis)
        return drawer

    def parent_coordinates(
        self, apply_offset: bool = False, **kwargs: Any
    ) -> tuple[str, int, int, int]:
        """returns seqid, start, stop, strand of this sequence on its parent

        Parameters
        ----------
        apply_offset
            if True, adds annotation offset from parent

        Notes
        -----
        seqid is the identifier of the parent. Returned coordinates are with
        respect to the plus strand, irrespective of whether the sequence has
        been reversed complemented or not.

        Returns
        -------
        seqid, start, end, strand of this sequence on the parent. strand is either
        -1 or 1.
        """
        strand = Strand.MINUS.value if self._seq.is_reversed else Strand.PLUS.value
        seqid, start, stop, _ = self._seq.parent_coords(apply_offset=apply_offset)
        return seqid, start, stop, strand

    def sample(
        self,
        *,
        n: int | None = None,
        with_replacement: bool = False,
        motif_length: int = 1,
        randint: Callable[
            [int, int | None, int | None], NumpyIntArrayType
        ] = numpy.random.randint,
        permutation: Callable[[int], NumpyIntArrayType] = numpy.random.permutation,
    ) -> Self:
        """Returns random sample of positions from self, e.g. to bootstrap.

        Parameters
        ----------
        n
            number of positions to sample. If None, all positions are sampled.
        with_replacement
            if True, samples with replacement.
        motif_length
            number of positions to sample as a single motif. Starting point
            of each sampled motif is modulo motif_length in the original sequence.
        randint
            random number generator, default is numpy.randint
        permutation
            function to generate a random permutation of positions, default is
            numpy.permutation

        Notes
        -----
        By default (resampling all positions without replacement), generates
        a permutation of the positions of the alignment.
        """
        population_size = len(self) // motif_length
        if not with_replacement and n and n > population_size:
            msg = f"cannot sample without replacement when {n=} > {population_size=}"
            raise ValueError(msg)

        n = n or population_size

        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            locations = permutation(population_size)[:n]

        if motif_length == 1:
            positions = locations
        else:
            positions = numpy.empty(n * motif_length, dtype=int)
            for i, loc in enumerate(locations):
                positions[i * motif_length : (i + 1) * motif_length] = range(
                    loc * motif_length,
                    (loc + 1) * motif_length,
                )

        sampled = numpy.array(self).take(positions)
        return cast(
            "Self",
            self.moltype.make_sequence(
                seq=sampled, name=f"{self.name}-randomised", check_seq=False
            ),
        )


class ProteinSequence(Sequence):
    """Holds the standard Protein sequence."""

    # constructed by PROTEIN moltype


class ByteSequence(Sequence):
    """Holds the bytes sequence."""

    # constructed by BYTES moltype


class ProteinWithStopSequence(Sequence):
    """Holds the standard Protein sequence, allows for stop codon."""

    # constructed by PROTEIN_WITH_STOP moltype


class NucleicAcidSequenceMixin(Sequence):
    """Mixin class for DNA and RNA sequences."""

    def can_pair(self, other: Self) -> bool:
        """Returns True if self and other could pair.

        other
            sequence, must be able to be cast to string.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.

        Gaps are only allowed to pair with other gaps.

        Weak pairs (like GU) are considered as possible pairs.
        """
        return self.moltype.can_pair(str(self), str(other))

    def can_mispair(self, other: Self) -> bool:
        """Returns True if any position in self could mispair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.

        Gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        return self.moltype.can_mispair(str(self), str(other))

    def must_pair(self, other: Self) -> bool:
        """Returns True if all positions in self must pair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.
        """
        return not self.moltype.can_mispair(str(self), str(other))

    def complement(self) -> Self:
        """Returns complement of self, using data from MolType.

        Always tries to return same type as item: if item looks like a dict,
        will return list of keys.
        """
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.complement(bytes(self)),
            info=self.info,
        )

    def reverse_complement(self) -> Self:
        """Converts a nucleic acid sequence to its reverse complement.
        Synonymn for rc."""
        return self.rc()

    def rc(self) -> Self:
        """Converts a nucleic acid sequence to its reverse complement."""
        return self.__class__(
            moltype=self.moltype,
            seq=self._seq[::-1],
            name=self.name,
            info=self.info,
            annotation_db=self.annotation_db,
        )

    def has_terminal_stop(
        self,
        gc: c3_genetic_code.GeneticCode | int = 1,
        strict: bool = False,
    ) -> bool:
        """Return True if the sequence has a terminal stop codon.

        Parameters
        ----------
        gc
            valid input to c3_genetic_code.get_code(), a genetic code object, number
            or name
        strict
            If True, raises an exception if length not divisible by 3
        """
        gc = c3_genetic_code.get_code(gc)
        _, s = self.parse_out_gaps()

        divisible_by_3 = len(s) % 3 == 0
        if divisible_by_3:
            end3 = str(s[-3:])
            return gc.is_stop(end3)

        if strict:
            msg = f"{self.name!r} length not divisible by 3"
            raise c3_alphabet.AlphabetError(msg)

        return False

    def trim_stop_codon(
        self,
        gc: c3_genetic_code.GeneticCode | int = 1,
        strict: bool = False,
    ) -> Self:
        """Removes a terminal stop codon from the sequence

        Parameters
        ----------
        gc
            valid input to c3_genetic_code.get_code(), a genetic code object, number
            or name
        strict
            If True, raises an exception if length not divisible by 3

        Notes
        -----
        If sequence contains gap characters, the result preserves the sequence
        length by adding gap characters at the end.
        """
        if not self.has_terminal_stop(gc=gc, strict=strict):
            return self

        gc = c3_genetic_code.get_code(gc)
        m, s = self.parse_out_gaps()

        divisible_by_3 = len(s) % 3 == 0
        if not divisible_by_3:
            return self

        end = str(s[-3:])

        if not gc.is_stop(end):
            return self

        if not m.num_gaps:
            # has zero length if no gaps
            return self[:-3]

        # determine terminal gap needed to fill in the sequence
        str_s = str(self)
        gaps = "".join(self.moltype.gaps)
        pattern = f"({'|'.join(gc['*'])})[{gaps}]*$"
        terminal_stop = re.compile(pattern)
        if match := terminal_stop.search(str_s):
            diff = len(str_s) - match.start()
            str_s = terminal_stop.sub("-" * diff, str_s)

        result = self.__class__(
            moltype=self.moltype,
            seq=str_s,
            name=self.name,
            info=self.info,
        )
        result.annotation_db = self.annotation_db
        return result

    def get_translation(
        self,
        gc: c3_genetic_code.GeneticCode | int = 1,
        incomplete_ok: bool = False,
        include_stop: bool = False,
        trim_stop: bool = True,
    ) -> ProteinSequence:
        """translate to amino acid sequence

        Parameters
        ----------
        gc
            valid input to cogent3.get_code(), a genetic code object, number
            or name
        incomplete_ok
            codons that are mixes of nucleotide and gaps converted to '-'.
            codons containing ambiguous nucleotides are translated as 'X'.
            raises a AlphabetError if False
        include_stop
            allows stop codons in translation
        trim_stop
            trims a terminal stop codon if it exists

        Returns
        -------
        sequence of PROTEIN moltype

        Raises
        ------
        AlphabetError if include_stop is False and a stop codon occurs
        """
        if not self.moltype.is_nucleic:
            msg = f"moltype must be a DNA/RNA, not {self.moltype.name!r}"
            raise c3_moltype.MolTypeError(
                msg,
            )

        protein = c3_moltype.get_moltype(
            "protein_with_stop" if include_stop else "protein",
        )
        gc = c3_genetic_code.get_code(gc)

        if trim_stop:
            seq = self.trim_stop_codon(gc=gc, strict=not incomplete_ok)
        else:
            seq = self

        # since we are realising the view, reverse complementing will be
        # dealt with, so rc=False
        pep = gc.translate(numpy.array(seq), rc=False, incomplete_ok=incomplete_ok)

        if not include_stop and "*" in pep:
            msg = f"{self.name!r} has a stop codon in the translation"
            raise c3_alphabet.AlphabetError(msg)

        if not incomplete_ok and "X" in pep:
            msg = (
                f"{self.name!r} has an incomplete codon or contains an ambiguity, set incomplete_ok=True to "
                "allow translation"
            )
            raise c3_alphabet.AlphabetError(
                msg,
            )
        return cast("ProteinSequence", protein.make_sequence(seq=pep, name=self.name))

    def to_rna(self) -> Sequence:
        """Returns copy of self as RNA."""
        return self.to_moltype("rna")

    def to_dna(self) -> Sequence:
        """Returns copy of self as DNA."""
        return self.to_moltype("dna")

    def strand_symmetry(self, motif_length: int = 1) -> TestResult:
        """returns G-test for strand symmetry"""
        counts = self.counts(motif_length=motif_length)
        ssym_pairs = self.moltype.strand_symmetric_motifs(motif_length=motif_length)

        obs: list[NumpyIntArrayType] = []
        motifs: list[str] = []
        for plus, minus in sorted(ssym_pairs):
            row = numpy.array([counts[plus], counts[minus]], dtype=int)
            if row.max() == 0:  # we ignore motifs missing on both strands
                continue

            obs.append(row)
            motifs.append(plus)

        d_array = DictArray.from_array_names(numpy.array(obs), motifs, ["+", "-"])
        cat = CategoryCounts(d_array)
        return cat.G_fit()


class DnaSequence(NucleicAcidSequenceMixin):
    """Holds the standard DNA sequence."""

    # constructed by DNA moltype


class RnaSequence(NucleicAcidSequenceMixin):
    """Holds the standard RNA sequence."""

    # constructed by RNA moltype


def _coerce_to_seqview(
    data: Aligned
    | SeqViewABC
    | Sequence
    | str
    | bytes
    | NumpyIntArrayType
    | tuple[str, ...]
    | list[str],
    seqid: str | None,
    alphabet: c3_alphabet.CharAlphabet[Any],
    offset: int,
) -> SeqViewABC:
    from cogent3.core.alignment import Aligned

    if isinstance(data, Sequence):
        data = data._seq

    if isinstance(data, SeqViewABC):
        # we require the indexes of shared states in alphabets to be the same
        # SeqView has an alphabet but SeqViewABC does NOT because that is
        # more general and covers the case where the SeqsData collection has the
        # alphabet
        if hasattr(data, "alphabet"):
            n = min(len(data.alphabet), len(alphabet))
            if data.alphabet[:n] != alphabet[:n]:
                msg = f"element order {data.alphabet=} != to that in {alphabet=} for {data=!r}"
                raise c3_alphabet.AlphabetError(
                    msg,
                )

        if offset and data.offset:
            msg = f"cannot set {offset=} on a SeqView with an offset {data.offset=}"
            raise ValueError(
                msg,
            )
        if offset:
            return data.with_offset(offset)
        return data

    if isinstance(data, (tuple, list)):
        data = "".join(data)

    if isinstance(data, Aligned):
        data = str(data)

    if isinstance(data, bytes):
        data = data.decode("utf8")

    if isinstance(data, str):
        return SeqView(
            parent=data,
            parent_len=len(data),
            seqid=seqid,
            alphabet=alphabet,
            offset=offset,
        )

    if isinstance(data, numpy.ndarray):
        return SeqView(
            parent=data.astype(alphabet.dtype),
            parent_len=len(data),
            seqid=seqid,
            alphabet=alphabet,
            offset=offset,
        )

    msg = f"{type(data)}"
    raise NotImplementedError(msg)


cls_map = {
    get_object_provenance(cls): cls
    for cls in (
        Sequence,
        DnaSequence,
        RnaSequence,
        ProteinSequence,
        ProteinWithStopSequence,
        ByteSequence,
    )
} | {
    "cogent3.core.sequence.Sequence": Sequence,
    "cogent3.core.sequence.DnaSequence": DnaSequence,
    "cogent3.core.sequence.RnaSequence": RnaSequence,
    "cogent3.core.sequence.ProteinSequence": ProteinSequence,
    "cogent3.core.sequence.ProteinWithStopSequence": ProteinWithStopSequence,
    "cogent3.core.sequence.ByteSequence": ByteSequence,
}


@register_deserialiser(*cls_map.keys())
def deserialise_sequence(data: dict[str, Any]) -> Sequence:
    cls = cls_map[data["type"]]
    return cls.from_rich_dict(data)

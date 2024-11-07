"""Contains classes that represent biological sequence data.

Notes
-----
These should be created via MolType.make_seq()
"""

from __future__ import annotations

import contextlib
import copy
import json
import os
import re
import typing
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from functools import singledispatch, total_ordering
from operator import eq, ne
from random import shuffle

import numpy
from numpy import array, floating, integer, issubdtype

from cogent3._version import __version__
from cogent3.core import new_alphabet, new_moltype
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    BasicAnnotationDb,
    FeatureDataType,
    GenbankAnnotationDb,
    SupportsFeatures,
    load_annotations,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import (
    FeatureMap,
    IndelMap,
    LostSpan,
    _input_vals_neg_step,
    _input_vals_pos_step,
)
from cogent3.format.fasta import seqs_to_fasta
from cogent3.maths.stats.contingency import CategoryCounts
from cogent3.maths.stats.number import CategoryCounter
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.util.misc import (
    DistanceFromMatrix,
    get_object_provenance,
    get_setting_from_environ,
)
from cogent3.util.transform import for_seq, per_shortest

OptGeneticCodeType = typing.Optional[typing.Union[int, str, "GeneticCode"]]
OptStr = typing.Optional[str]
OptInt = typing.Optional[int]
OptFloat = typing.Optional[float]
IntORFloat = typing.Union[int, float]
StrORIterableStr = typing.Union[str, typing.Iterable[str]]
StrORBytesORArray = typing.Union[str, bytes, numpy.ndarray]
ARRAY_TYPE = type(array(1))
DEFAULT_ANNOTATION_DB = BasicAnnotationDb

# standard distance functions: left  because generally useful
frac_same = for_seq(f=eq, aggregator=sum, normalizer=per_shortest)
frac_diff = for_seq(f=ne, aggregator=sum, normalizer=per_shortest)


def _is_int(val) -> bool:
    """whether val is builtin, or numpy, integer"""
    return issubdtype(val.__class__, integer) or isinstance(val, int)


def _is_float(val) -> bool:
    """whether val is builtin, or numpy, integer"""
    return issubdtype(val.__class__, floating) or isinstance(val, float)


def _moltype_seq_from_rich_dict(data):
    from cogent3.core import new_moltype

    data.pop("type")
    data.pop("version")
    data.pop("annotation_db", None)
    offset = data.pop("annotation_offset", 0)

    moltype = data.pop("moltype")
    moltype = new_moltype.get_moltype(moltype)

    seqview_data = data.pop("seq")
    seq = _coerce_to_seqview(
        seqview_data["init_args"]["parent"],
        data["name"],
        moltype.most_degen_alphabet(),
        offset,
    )
    seq = seq[:: seqview_data["init_args"]["slice_record"]["init_args"]["step"]]

    return moltype, seq


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
class Sequence:
    """Holds the standard Sequence object. Immutable.

    Notes
    -----
    Sequences should be constructed by a MolType instance.
    """

    def __init__(
        self,
        moltype: "MolType",
        seq: typing.Union[StrORBytesORArray, SeqViewABC],
        *,
        name: OptStr = None,
        info: typing.Optional[typing.Union[dict, InfoClass]] = None,
        annotation_offset: int = 0,
    ):
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
        """
        self.moltype = moltype
        self.name = name
        self._seq = _coerce_to_seqview(
            seq, name, self.moltype.most_degen_alphabet(), annotation_offset
        )

        info = info or {}
        self.info = InfoClass(**info)
        self._repr_policy = dict(num_pos=60)
        self._annotation_db = DEFAULT_ANNOTATION_DB()

    def __str__(self):
        result = str(self._seq)
        if self._seq.is_reversed:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def __bytes__(self):
        result = bytes(self._seq)
        if self._seq.is_reversed:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def __array__(self, dtype=None, copy=None) -> numpy.ndarray[int]:
        result = array(self._seq, dtype=dtype)
        if self._seq.is_reversed:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def to_array(self, apply_transforms: bool = True) -> numpy.ndarray[int]:
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

    def to_fasta(self, make_seqlabel=None, block_size=60) -> str:
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
        self, exclude_annotations: bool = True
    ) -> dict[str, str | dict[str, str]]:
        """returns {'name': name, 'seq': sequence, 'moltype': moltype.label}

        Notes
        -----
        Deserialisation of the sequence object will not include the annotation_db
        even if exclude_annotations=False.
        """
        info = {} if self.info is None else self.info
        if info.get("Refs", None) is not None and "Refs" in info:
            info.pop("Refs")

        info = info or None
        seq = self._seq.to_rich_dict() if hasattr(self, "_seq") else str(self)
        data = dict(
            name=self.name,
            seq=seq,
            moltype=self.moltype.label,
            info=info,
            type=get_object_provenance(self),
            version=__version__,
        )
        if hasattr(self, "annotation_offset"):
            offset = int(self._seq.slice_record.parent_start)
            data |= dict(annotation_offset=offset)

        if (
            hasattr(self, "annotation_db")
            and self.annotation_db
            and not exclude_annotations
        ):
            data["annotation_db"] = self.annotation_db.to_rich_dict()

        return data

    @classmethod
    def from_rich_dict(cls, data: dict):
        moltype, seq = _moltype_seq_from_rich_dict(data)

        return cls(moltype=moltype, seq=seq, **data)

    def to_json(self) -> str:
        """returns a json formatted string"""
        return json.dumps(self.to_rich_dict())

    def count(self, item: str):
        """count() delegates to self._seq."""
        return str(self).count(item)

    def counts(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
        warn: bool = False,
    ) -> CategoryCounter:
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
        exclude_unobserved
            if True, unobserved motif combinations are excluded.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """
        # refactor: array

        data = str(self)

        if motif_length == 1:
            counts = CategoryCounter(data)
        else:
            if warn and len(data) % motif_length != 0:
                warnings.warn(
                    f"{self.name} length not divisible by {motif_length}, truncating"
                )
            limit = (len(data) // motif_length) * motif_length
            data = data[:limit]

            counts = CategoryCounter(
                data[i : i + motif_length] for i in range(0, limit, motif_length)
            )

        exclude = []
        if not include_ambiguity or not allow_gap:
            is_degen = self.moltype.is_degenerate
            is_gap = self.moltype.is_gapped
            for motif in counts:
                if not include_ambiguity and is_degen(motif):
                    exclude.append(motif)
                elif not allow_gap and is_gap(motif):
                    exclude.append(motif)

        for motif in exclude:
            del counts[motif]

        return counts

    def __lt__(self, other):
        """compares based on the sequence string."""
        return str(self) < str(other)

    def __eq__(self, other):
        """compares based on the sequence string."""
        return str(self) == str(other)

    def __ne__(self, other):
        """compares based on the sequence string."""
        return str(self) != str(other)

    def __hash__(self):
        """__hash__ behaves like the sequence string for dict lookup."""
        return hash(str(self))

    def __contains__(self, other):
        """__contains__ checks whether other is in the sequence string."""
        return other in str(self)

    def shuffle(self):
        """returns a randomized copy of the Sequence object"""
        randomized_copy_list = list(self)
        shuffle(randomized_copy_list)
        return self.__class__(
            moltype=self.moltype, seq="".join(randomized_copy_list), info=self.info
        )

    def strip_degenerate(self):
        """Removes degenerate bases by stripping them out of the sequence."""
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.strip_degenerate(bytes(self)),
            info=self.info,
        )

    def strip_bad(self):
        """Removes any symbols not in the alphabet."""
        return self.__class__(
            moltype=self.moltype, seq=self.moltype.strip_bad(str(self)), info=self.info
        )

    def strip_bad_and_gaps(self):
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
        return self.moltype.is_valid(self)

    def is_strict(self) -> bool:
        """Returns True if sequence contains only monomers."""
        return self.moltype.alphabet.is_valid(array(self))

    def disambiguate(self, method: str = "strip"):
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

    def degap(self):
        """Deletes all gap characters from sequence."""
        result = self.__class__(
            moltype=self.moltype,
            seq=self.moltype.degap(bytes(self)),
            name=self.name,
            info=self.info,
        )
        result.annotation_db = self.annotation_db
        return result

    def gap_indices(self) -> numpy.ndarray:
        """Returns array of the indices of all gaps in the sequence"""
        return numpy.where(self.gap_vector())[0]

    def gap_vector(self) -> list[bool]:
        """Returns vector of True or False according to which pos are gaps or missing."""
        return (
            (array(self) == self.moltype.degen_gapped_alphabet.gap_index)
            | (array(self) == self.moltype.degen_gapped_alphabet.missing_index)
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

    def count_variants(self):
        """Counts number of possible sequences matching the sequence, given
        any ambiguous characters in the sequence.

        Notes
        -----
        Uses self.ambiguitues to decide how many possibilities there are at
        each position in the sequence and calculates the permutations.
        """
        return self.moltype.count_variants(str(self))

    def mw(self, method: str = "random", delta: OptFloat = None) -> float:
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
        return self.moltype.mw(self, method, delta)

    def can_match(self, other) -> bool:
        """Returns True if every pos in self could match same pos in other.

        Truncates at length of shorter sequence.
        gaps are only allowed to match other gaps.
        """
        return self.moltype.can_match(str(self), str(other))

    def diff(self, other) -> int:
        """Returns number of differences between self and other.

        Notes
        -----
        Truncates at the length of the shorter sequence.
        """
        return self.distance(other)

    def distance(
        self, other, function: typing.Callable[[str, str], IntORFloat] = None
    ) -> IntORFloat:
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
            function = lambda a, b: a != b
        distance = 0
        for first, second in zip(self, other):
            distance += function(first, second)
        return distance

    def matrix_distance(self, other, matrix) -> IntORFloat:
        """Returns distance between self and other using a score matrix.

        Warning
        -------
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

    def frac_same(self, other) -> float:
        """Returns fraction of positions where self and other are the same.

        Notes
        -----
        Truncates at length of shorter sequence. Will return  0 if one sequence
        is empty.
        """
        return frac_same(self, other)

    def frac_diff(self, other) -> float:
        """Returns fraction of positions where self and other differ.

        Notes
        -----
        Truncates at length of shorter sequence. Will return  0 if one sequence
        is empty.
        """
        return frac_diff(self, other)

    def frac_same_gaps(self, other):
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
        return sum([is_gap(i) == is_gap(j) for i, j in zip(self, other)]) / min(
            len(self), len(other)
        )

    def frac_diff_gaps(self, other):
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

    def frac_same_non_gaps(self, other):
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
        for i, j in zip(self, other):
            if is_gap(i) or is_gap(j):
                continue
            count += 1
            if i == j:
                identities += 1

        if count:
            return identities / count
        else:  # there were no positions that weren't gaps
            return 0

    def frac_diff_non_gaps(self, other):
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
        for i, j in zip(self, other):
            if is_gap(i) or is_gap(j):
                continue
            count += 1
            if i != j:
                diffs += 1

        if count:
            return diffs / count
        else:  # there were no positions that weren't gaps
            return 0

    def frac_similar(self, other, similar_pairs: dict[(str, str), typing.Any]):
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

        return for_seq(f=lambda x, y: (x, y) in similar_pairs, normalizer=per_shortest)(
            self, other
        )

    def with_termini_unknown(self):
        """Returns copy of sequence with terminal gaps remapped as missing."""
        gaps = self.gap_vector()
        first_nongap = last_nongap = None
        for i, state in enumerate(gaps):
            if not state:
                if first_nongap is None:
                    first_nongap = i
                last_nongap = i
        missing = self.moltype.missing
        if first_nongap is None:
            return self.__class__(
                moltype=self.moltype, seq=missing * len(self), info=self.info
            )
        prefix = missing * first_nongap
        mid = str(self)[first_nongap : last_nongap + 1]
        suffix = missing * (len(self) - last_nongap - 1)
        return self.__class__(
            moltype=self.moltype, seq=prefix + mid + suffix, info=self.info
        )

    def _repr_html_(self):
        settings = self._repr_policy.copy()
        env_vals = get_setting_from_environ(
            "COGENT3_ALIGNMENT_REPR_POLICY",
            dict(num_pos=int),
        )
        settings.update(env_vals)
        return self.to_html(limit=settings["num_pos"])

    def to_html(
        self,
        wrap: int = 60,
        limit: OptInt = None,
        colors: typing.Optional[typing.Mapping[str, str]] = None,
        font_size: int = 12,
        font_family: str = "Lucida Console",
    ):
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
            colors=colors, font_size=font_size, font_family=font_family
        )

        seq = str(self)
        seq = seq if limit is None else seq[:limit]
        gaps = "".join(self.moltype.gaps)
        seqlen = len(seq)
        start_gap = re.search(f"^[{gaps}]+", "".join(seq))
        end_gap = re.search(f"[{gaps}]+$", "".join(seq))

        start = 0 if start_gap is None else start_gap.end()
        end = len(seq) if end_gap is None else end_gap.start()
        seq_style = []
        template = '<span class="%s">%%s</span>'
        styled_seq = []
        for i in range(seqlen):
            char = seq[i]
            if i < start or i >= end:
                style = f"terminal_ambig_{self.moltype.label}"
            else:
                style = styles[char]

            seq_style.append(template % style)
            styled_seq.append(seq_style[-1] % char)

        # make a html table
        seq = array(styled_seq, dtype="O")
        table = ["<table>"]
        seq_ = "<td>%s</td>"
        label_ = '<td class="label">%s</td>'
        num_row_ = '<tr class="num_row"><td></td><td><b>{:,d}</b></td></tr>'
        for i in range(0, seqlen, wrap):
            table.append(num_row_.format(i))
            seqblock = seq[i : i + wrap].tolist()
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
            ".c3seq .label { font-size: %dpt ; text-align: right !important; "
            "color: black !important; padding: 0 4px; }" % font_size,
            "\n".join([".c3seq " + style for style in css]),
            "</style>",
            '<div class="c3seq">',
            "\n".join(table),
            f"<p><i>{summary}</i></p>",
            "</div>",
        ]
        return "\n".join(text)

    def __add__(self, other):
        """Adds two sequences (other can be a string as well)."""
        if hasattr(other, "moltype"):
            if self.moltype != other.moltype:
                raise ValueError(
                    f"MolTypes don't match: ({self.moltype},{other.moltype})"
                )
        other_seq = str(other)

        if not self.moltype.most_degen_alphabet().is_valid(str(other)):
            raise new_alphabet.AlphabetError(
                f"Invalid sequence characters in other for moltype={self.moltype.label}"
            )

        # If two sequences with the same name are being added together the name should not be None
        if type(other) == type(self):
            name = self.name if self.name == other.name else None
        else:
            name = None

        return self.__class__(
            moltype=self.moltype, seq=str(self) + other_seq, name=name
        )

    @property
    def annotation_offset(self):
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

        return self._seq.slice_record.parent_start

    @property
    def annotation_db(self):
        return self._annotation_db

    @annotation_db.setter
    def annotation_db(self, value: SupportsFeatures):
        # Without knowing the contents of the db we cannot
        # establish whether self.moltype is compatible, so
        # we rely on the user to get that correct
        # one approach to support validation might be to add
        # to the SupportsFeatures protocol a is_nucleic flag,
        # for both DNA and RNA. But if a user trys get_slice()
        # on a '-' strand feature, they will get a TypeError.
        # I think that's enough.
        self.replace_annotation_db(value, check=False)

    def replace_annotation_db(
        self, value: SupportsFeatures, check: bool = True
    ) -> None:
        """public interface to assigning the annotation_db

        Parameters
        ----------
        value
            the annotation db instance
        check
            whether to check value supports the feature interface

        Notes
        -----
        The check can be very expensive, so if you're confident set it to False
        """
        if value == self._annotation_db:
            return

        if check and value and not isinstance(value, SupportsFeatures):
            raise TypeError(f"{type(value)} does not satisfy SupportsFeatures")

        self._annotation_db = value

    def get_features(
        self,
        *,
        biotype: OptStr = None,
        name: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        allow_partial: bool = False,
    ):
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

        # note: offset is handled by absolute_position
        # we set include_boundary=False because start is inclusive indexing,
        # i,e., the start cannot be equal to the length of the view
        (
            # parent_id can differ from self.name if there was a
            # rename operation
            parent_id,
            *_,
            strand,
        ) = self.parent_coordinates()
        query_start = self._seq.slice_record.absolute_position(
            start, include_boundary=False
        )
        # we set include_boundary=True because stop is exclusive indexing,
        # i,e., the stop can be equal to the length of the view
        query_stop = self._seq.slice_record.absolute_position(
            stop, include_boundary=True
        )

        rev_strand = strand == -1
        if rev_strand:
            query_start, query_stop = query_stop, query_start

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

        for feature in self.annotation_db.get_features_matching(
            seqid=parent_id,
            name=name,
            biotype=biotype,
            start=query_start,
            stop=query_stop,
            allow_partial=allow_partial,
        ):
            # spans need to be converted from absolute to relative positions
            # DO NOT do adjustment in make_feature since that's user facing,
            # and we expect them to make a feature manually wrt to their
            # current view
            spans = array(feature.pop("spans"), dtype=int)
            for i, v in enumerate(spans.ravel()):
                rel_pos = self._seq.slice_record.relative_position(v)
                spans.ravel()[i] = rel_pos

            if rev_strand:
                # see above comment
                spans = len(self) - spans

            feature["spans"] = spans.tolist()
            yield self.make_feature(feature)

    def _relative_spans(self, spans):
        r_spans = []

        for s, e in spans:
            # reverse feature determined by absolute position
            reverse = s > e

            if reverse:
                start = self._seq.relative_position(s, stop=True)
                end = self._seq.relative_position(e)
            else:
                start = self._seq.relative_position(s)
                end = self._seq.relative_position(e, stop=True)

            r_spans += [(start, end)]

        return r_spans

    def make_feature(self, feature: FeatureDataType, *args) -> Feature:
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
        feature = dict(feature)
        seq_rced = self._seq.is_reversed
        spans = feature.pop("spans", None)
        revd = feature.pop("strand", None) == "-"
        feature["strand"] = "+" if revd == seq_rced else "-"

        vals = array(spans)
        pre = abs(vals.min()) if vals.min() < 0 else 0
        post = abs(vals.max() - len(self)) if vals.max() > len(self) else 0

        # we find the spans > len(self)
        new_spans = []
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
            spans = list(fmap.spans)
            if pre:
                spans.insert(0, LostSpan(pre))
            if post:
                spans.append(LostSpan(post))
            fmap = FeatureMap(spans=spans, parent_length=len(self))

        if seq_rced:
            fmap = fmap.nucleic_reversed()

        feature.pop("on_alignment", None)
        feature.pop("seqid", None)
        return Feature(parent=self, seqid=self.name, map=fmap, **feature)

    def annotate_from_gff(self, f: os.PathLike, offset=None):
        """copies annotations from a gff file to self,

        Parameters
        ----------
        f : path to gff annotation file.
        offset : Optional, the offset between annotation coordinates and sequence coordinates.
        """
        if isinstance(self.annotation_db, GenbankAnnotationDb):
            raise ValueError("GenbankAnnotationDb already attached")

        self.annotation_db = load_annotations(
            path=f, seqids=self.name, db=self.annotation_db
        )

        if offset:
            self.annotation_offset = offset

    def add_feature(
        self,
        *,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]],
        parent_id: OptStr = None,
        strand: OptStr = None,
        on_alignment: bool = False,
        seqid: OptStr = None,
    ) -> Feature:
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
        feature_data = FeatureDataType(
            seqid=self.name,
            **{n: v for n, v in locals().items() if n not in ("self", "seqid")},
        )

        self.annotation_db.add_feature(**feature_data)
        for discard in ("on_alignment", "parent_id"):
            feature_data.pop(discard)
        return self.make_feature(feature_data)

    def to_moltype(self, moltype: typing.Union[str, new_moltype.MolType]) -> "Sequence":
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
            raise ValueError(f"unknown moltype '{moltype}'")

        moltype = new_moltype.get_moltype(moltype)

        if moltype is self.moltype:
            return self

        if len(self.moltype.alphabet) == len(moltype.alphabet):
            # converting between dna and rna
            seq = (
                moltype.most_degen_alphabet()
                .array_to_bytes(array(self))
                .decode("utf-8")
            )

        else:
            # assume converting between ASCII/BYTES and nucleic acid
            seq = str(self)

        if not moltype.most_degen_alphabet().is_valid(seq):
            raise ValueError(
                f"Changing from old moltype={self.moltype.label!r} to new "
                f"moltype={moltype.label!r} is not valid for this data"
            )
        sv = SeqView(parent=seq, alphabet=moltype.most_degen_alphabet())
        new = self.__class__(moltype=moltype, seq=sv, name=self.name, info=self.info)
        new.annotation_db = self.annotation_db
        return new

    def copy_annotations(self, seq_db: SupportsFeatures) -> None:
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
            raise TypeError(
                f"type {type(seq_db)} does not match SupportsFeatures interface"
            )

        if not seq_db.num_matches(seqid=self.name):
            return

        if self.annotation_db and not self.annotation_db.compatible(seq_db):
            raise TypeError(f"type {type(seq_db)} != {type(self.annotation_db)}")

        if self.annotation_db is None:
            self.annotation_db = type(seq_db)()

        self.annotation_db.update(seq_db, seqids=self.name)

    def copy(self, exclude_annotations: bool = False, sliced: bool = True):
        """returns a copy of self

        Parameters
        -----------
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
        new = self.__class__(
            moltype=self.moltype,
            seq=data,
            name=self.name,
            info=self.info,
            annotation_offset=offset,
        )
        db = None if exclude_annotations else copy.deepcopy(self.annotation_db)
        new._annotation_db = db
        return new

    def with_masked_annotations(
        self,
        annot_types: StrORIterableStr,
        mask_char: str = None,
        shadow: bool = False,
        extend_query: bool = False,
    ):
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
            mask_char = (
                self.moltype.missing
                or max(self.moltype.ambiguities.items(), key=lambda x: len(x[1]))[0]
            )
        assert (
            mask_char in self.moltype.most_degen_alphabet()
        ), f"Invalid mask_char {mask_char}"

        annotations = []
        annot_types = [annot_types] if isinstance(annot_types, str) else annot_types
        for annot_type in annot_types:
            annotations += list(
                self.get_features(biotype=annot_type, allow_partial=True)
            )

        if not annotations:
            region = Feature(
                parent=self,
                seqid=self.name,
                name=None,
                biotype=None,
                map=FeatureMap.from_locations(locations=[], parent_length=len(self)),
                strand="+",
            )
        else:
            region = annotations[0].union(annotations[1:])

        if shadow:
            region = region.shadow()

        i = 0
        segments = []
        coords = region.map.get_coordinates()
        for b, e in coords:
            segments.extend((str(self[i:b]), mask_char * (e - b)))
            i = e
        segments.append(str(self[i:]))

        new = self.__class__(
            moltype=self.moltype, seq="".join(segments), name=self.name, info=self.info
        )
        new.annotation_db = self.annotation_db
        return new

    def gapped_by_map_segment_iter(
        self, segment_map: IndelMap, allow_gaps: bool = True, recode_gaps: bool = False
    ) -> typing.Iterator[str, str, str]:
        if not allow_gaps and not segment_map.complete:
            raise ValueError(f"gap(s) in map {segment_map}")

        for span in segment_map.spans:
            if span.lost:
                unknown = "?" if span.terminal or recode_gaps else "-"
                seg = unknown * span.length
            else:
                seg = str(self[span.start : span.end])

            yield seg

    def gapped_by_map_motif_iter(
        self, segment_map: IndelMap
    ) -> typing.Iterator[str, str, str]:
        for segment in self.gapped_by_map_segment_iter(segment_map):
            yield from segment

    def gapped_by_map(self, segment_map: IndelMap, recode_gaps: bool = False):
        segments = self.gapped_by_map_segment_iter(segment_map, True, recode_gaps)
        return self.__class__(
            moltype=self.moltype,
            seq="".join(segments),
            name=self.name,
            info=self.info,
        )

    def _mapped(self, segment_map: IndelMap):
        # Called by generic __getitem__
        if segment_map.num_spans == 1:
            seq = self._seq[segment_map.start : segment_map.end]
            offset = segment_map.start
        else:
            segments = self.gapped_by_map_segment_iter(segment_map, allow_gaps=False)
            seq = "".join(segments)
            offset = 0

        return self.__class__(
            moltype=self.moltype,
            seq=seq,
            name=self.name,
            info=self.info,
            annotation_offset=offset,
        )

    def __repr__(self):
        myclass = f"{self.__class__.__name__}"
        myclass = myclass.split(".")[-1]
        seq = f"{str(self)[:7]}... {len(self):,}" if len(self) > 10 else str(self)
        return f"{myclass}({seq})"

    def __getitem__(self, index):
        preserve_offset = False
        if hasattr(index, "get_slice"):
            if index.parent is not self:
                raise ValueError("cannot slice Feature not bound to self")
            return index.get_slice()

        if hasattr(index, "map"):
            index = index.map

        if isinstance(index, (FeatureMap, IndelMap)):
            new = self._mapped(index)
            # annotations have no meaning if disjoint slicing segments
            preserve_offset = index.num_spans == 1

        elif isinstance(index, slice) or _is_int(index):
            new = self.__class__(
                moltype=self.moltype,
                seq=self._seq[index],
                name=self.name,
                info=self.info,
            )
            stride = getattr(index, "step", 1) or 1
            preserve_offset = stride > 0

        if isinstance(index, (list, tuple)):
            raise TypeError("cannot slice using list or tuple")

        if self.annotation_db is not None and preserve_offset:
            new.replace_annotation_db(self.annotation_db, check=False)

        if _is_float(index):
            raise TypeError("cannot slice using float")

        if hasattr(self, "_repr_policy"):
            new._repr_policy.update(self._repr_policy)

        return new

    def __iter__(self):
        yield from iter(str(self))

    def get_name(self):
        """Return the sequence name -- should just use name instead."""
        return self.name

    def __len__(self):
        return len(self._seq)

    def get_type(self):
        """Return the sequence type as moltype label."""
        return self.moltype.label

    def resolved_ambiguities(self) -> list[set[str]]:
        """Returns a list of sets of strings."""
        ambigs = self.moltype.ambiguities
        return [set(ambigs.get(motif, motif)) for motif in str(self)]

    def iter_kmers(self, k: int, strict: bool = True) -> typing.Iterator[str]:
        """generates all overlapping k-mers.
        When strict is True, the characters in the k-mer must be
        a subset of the canonical characters for the moltype"""
        if k <= 0:
            raise ValueError(f"k must be an int > 0, not {k}")

        if not isinstance(k, int):
            raise ValueError(f"k must be an int, not {k}")

        canonical = set(self.moltype)
        seq = str(self)
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if not strict:
                yield kmer
            else:
                if all(char in canonical for char in kmer):
                    yield kmer

    def get_kmers(self, k: int, strict: bool = True) -> list[str]:
        """return all overlapping k-mers"""
        return list(self.iter_kmers(k, strict))

    def sliding_windows(self, window, step, start=None, end=None):
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
        start = [start, 0][start is None]
        end = [end, len(self) - window + 1][end is None]
        end = min(len(self) - window + 1, end)
        if start < end and len(self) - end >= window - 1:
            for pos in range(start, end, step):
                yield self[pos : pos + window]

    def get_in_motif_size(self, motif_length=1, warn=False):
        """returns sequence as list of non-overlapping motifs

        Parameters
        ----------
        motif_length
            length of the motifs
        warn
            whether to notify of an incomplete terminal motif
        """
        seq = self._seq
        if isinstance(seq, SeqViewABC):
            seq = str(self)
        if motif_length == 1:
            return seq

        length = len(seq)
        remainder = length % motif_length
        if remainder and warn:
            warnings.warn(
                f'Dropped remainder "{seq[-remainder:]}" from end of sequence'
            )
        return [
            seq[i : i + motif_length]
            for i in range(0, length - remainder, motif_length)
        ]

    def parse_out_gaps(self):
        """returns Map corresponding to gap locations and ungapped Sequence"""
        gap = re.compile(f"[{re.escape(self.moltype.gap)}]+")
        seq = str(self)
        gap_pos = []
        cum_lengths = []
        for match in gap.finditer(seq):
            pos = match.start()
            gap_pos.append(pos)
            cum_lengths.append(match.end() - pos)

        gap_pos = array(gap_pos)
        cum_lengths = array(cum_lengths).cumsum()
        gap_pos[1:] = gap_pos[1:] - cum_lengths[:-1]

        seq = self.__class__(
            moltype=self.moltype,
            seq=gap.sub("", seq),
            name=self.get_name(),
            info=self.info,
        )
        indel_map = IndelMap(
            gap_pos=gap_pos, cum_gap_lengths=cum_lengths, parent_length=len(seq)
        )
        seq.annotation_db = self.annotation_db
        return indel_map, seq

    def is_annotated(
        self, biotype: typing.Optional[typing.Union[str, tuple[str]]] = None
    ) -> bool:
        """returns True if sequence parent name has any annotations

        Parameters
        ----------
        biotype
            amend condition to return True only if the sequence is
            annotated with one of provided biotypes.
        """
        with contextlib.suppress(AttributeError):
            return (
                self.annotation_db.num_matches(seqid=self._seq.seqid, biotype=biotype)
                != 0
            )
        return False

    def annotate_matches_to(
        self, pattern: str, biotype: str, name: str, allow_multiple: bool = False
    ):
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
        self, *, biotype: typing.Optional[StrORIterableStr] = None
    ) -> dict:
        """returns a dict of drawables, keyed by type

        Parameters
        ----------
        biotype
            passed to get_features(biotype). Can be a single biotype or
            series. Only features matching this will be included.
        """
        # make sure the drawables are unique by adding to a set
        features = set(self.get_features(biotype=biotype, allow_partial=True))
        result = defaultdict(list)
        for f in features:
            result[f.biotype].append(f.get_drawable())
        return result

    def get_drawable(
        self,
        *,
        biotype: typing.Optional[StrORIterableStr] = None,
        width: int = 600,
        vertical: bool = False,
    ):
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
        top = 0
        space = 0.25
        annotes = []
        for feature_type in biotypes:
            new_bottom = top + space
            for i, annott in enumerate(drawables[feature_type]):
                annott.shift(y=new_bottom - annott.bottom)
                if i > 0:
                    annott._showlegend = False
                annotes.append(annott)

            top = annott.top

        top += space
        height = max((top / len(self)) * width, 300)
        xaxis = dict(range=[0, len(self)], zeroline=False, showline=True)
        yaxis = dict(range=[0, top], visible=False, zeroline=True, showline=True)

        if vertical:
            all_traces = [t.T.as_trace() for t in annotes]
            width, height = height, width
            xaxis, yaxis = yaxis, xaxis
        else:
            all_traces = [t.as_trace() for t in annotes]

        drawer = Drawable(
            title=self.name, traces=all_traces, width=width, height=height
        )
        drawer.layout.update(xaxis=xaxis, yaxis=yaxis)
        return drawer

    def parent_coordinates(self) -> tuple[str, int, int, int]:
        """returns seqid, start, stop, strand of this sequence on its parent

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
        strand = -1 if self._seq.is_reversed else 1
        return (
            self._seq.seqid,
            self._seq.slice_record.parent_start,
            self._seq.slice_record.parent_stop,
            strand,
        )


@register_deserialiser(get_object_provenance(Sequence))
def deserialise_sequence(data) -> Sequence:
    return Sequence.from_rich_dict(data)


class ProteinSequence(Sequence):
    """Holds the standard Protein sequence."""

    # constructed by PROTEIN moltype


@register_deserialiser(get_object_provenance(ProteinSequence))
def deserialise_protein_sequence(data) -> ProteinSequence:
    return ProteinSequence.from_rich_dict(data)


class ByteSequence(Sequence):
    """Holds the bytes sequence."""

    # constructed by BYTES moltype


@register_deserialiser(get_object_provenance(ByteSequence))
def deserialise_bytes_sequence(data) -> ByteSequence:
    return ByteSequence.from_rich_dict(data)


class ProteinWithStopSequence(Sequence):
    """Holds the standard Protein sequence, allows for stop codon."""

    # constructed by PROTEIN_WITH_STOP moltype


@register_deserialiser(get_object_provenance(ProteinWithStopSequence))
def deserialise_protein_with_stop_sequence(data) -> ProteinWithStopSequence:
    return ProteinWithStopSequence.from_rich_dict(data)


class NucleicAcidSequenceMixin:
    """Mixin class for DNA and RNA sequences."""

    def __str__(self):
        result = str(self._seq)
        if self._seq.is_reversed:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def can_pair(self, other) -> bool:
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

    def can_mispair(self, other) -> bool:
        """Returns True if any position in self could mispair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.

        Gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        return self.moltype.can_mispair(self, str(other))

    def must_pair(self, other) -> bool:
        """Returns True if all positions in self must pair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.
        """
        return not self.moltype.can_mispair(self, str(other))

    def complement(self):
        """Returns complement of self, using data from MolType.

        Always tries to return same type as item: if item looks like a dict,
        will return list of keys.
        """
        return self.__class__(
            moltype=self.moltype,
            seq=self.moltype.complement(bytes(self)),
            info=self.info,
        )

    def reverse_complement(self):
        """Converts a nucleic acid sequence to its reverse complement.
        Synonymn for rc."""
        return self.rc()

    def rc(self):
        """Converts a nucleic acid sequence to its reverse complement."""
        rc = self.__class__(
            moltype=self.moltype,
            seq=self._seq[::-1],
            name=self.name,
            info=self.info,
        )
        rc.annotation_db = self.annotation_db
        return rc

    def has_terminal_stop(
        self,
        gc: OptGeneticCodeType = None,
        strict: bool = False,
    ) -> bool:
        """Return True if the sequence has a terminal stop codon.

        Parameters
        ----------
        gc
            valid input to new_genetic_code.get_code(), a genetic code object, number
            or name
        strict
            If True, raises an exception if length not divisible by 3
        """
        from cogent3.core import new_genetic_code

        gc = new_genetic_code.get_code(gc)
        _, s = self.parse_out_gaps()

        divisible_by_3 = len(s) % 3 == 0
        if divisible_by_3:
            end3 = str(s[-3:])
            return gc.is_stop(end3)

        if strict:
            raise new_alphabet.AlphabetError(f"{self.name!r} length not divisible by 3")

        return False

    def trim_stop_codon(
        self,
        gc: OptGeneticCodeType = None,
        strict: bool = False,
    ):
        """Removes a terminal stop codon from the sequence

        Parameters
        ----------
        gc
            valid input to new_genetic_code.get_code(), a genetic code object, number
            or name
        strict
            If True, raises an exception if length not divisible by 3

        Notes
        -----
        If sequence contains gap characters, the result preserves the sequence
        length by adding gap characters at the end.
        """
        from cogent3.core import new_genetic_code

        if not self.has_terminal_stop(gc=gc, strict=strict):
            return self

        gc = new_genetic_code.get_code(gc)
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
        s = str(self)
        gaps = "".join(self.moltype.gaps)
        pattern = f"({'|'.join(gc['*'])})[{gaps}]*$"
        terminal_stop = re.compile(pattern)
        if match := terminal_stop.search(s):
            diff = len(s) - match.start()
            s = terminal_stop.sub("-" * diff, s)

        result = self.__class__(
            moltype=self.moltype, seq=s, name=self.name, info=self.info
        )
        result.annotation_db = self.annotation_db
        return result

    def get_translation(
        self,
        gc: OptGeneticCodeType = 1,
        incomplete_ok: bool = False,
        include_stop: bool = False,
        trim_stop: bool = True,
    ):
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

        from cogent3.core import new_genetic_code, new_moltype

        protein = new_moltype.get_moltype(
            "protein_with_stop" if include_stop else "protein"
        )
        gc = new_genetic_code.get_code(gc)

        if trim_stop:
            seq = self.trim_stop_codon(gc=gc, strict=not incomplete_ok)
        else:
            seq = self

        # since we are realising the view, reverse complementing will be
        # dealt with, so rc=False
        pep = gc.translate(array(seq), rc=False, incomplete_ok=incomplete_ok)

        if not include_stop and "*" in pep:
            raise new_alphabet.AlphabetError("stop codon in translation")

        if not incomplete_ok and "X" in pep:
            raise new_alphabet.AlphabetError(
                "Incomplete codon in translation, set incomplete_ok=True to "
                "allow translation"
            )
        return protein.make_seq(seq=pep, name=self.name)

    def to_rna(self):
        """Returns copy of self as RNA."""
        return self.to_moltype("rna")

    def to_dna(self):
        """Returns copy of self as DNA."""
        return self.to_moltype("dna")

    def strand_symmetry(self, motif_length=1):
        """returns G-test for strand symmetry"""
        counts = self.counts(motif_length=motif_length)
        ssym_pairs = self.moltype.strand_symmetric_motifs(motif_length=motif_length)

        motifs = []
        obs = []
        for plus, minus in sorted(ssym_pairs):
            row = array([counts[plus], counts[minus]], dtype=int)
            if row.max() == 0:  # we ignore motifs missing on both strands
                continue

            obs.append(row)
            motifs.append(plus)

        template = DictArrayTemplate(motifs, ["+", "-"])
        obs = template.wrap(obs)
        cat = CategoryCounts(obs)
        return cat.G_fit()


class DnaSequence(Sequence, NucleicAcidSequenceMixin):
    """Holds the standard DNA sequence."""

    # constructed by DNA moltype


@register_deserialiser(get_object_provenance(DnaSequence))
def deserialise_dna_sequence(data) -> DnaSequence:
    return DnaSequence.from_rich_dict(data)


class RnaSequence(Sequence, NucleicAcidSequenceMixin):
    """Holds the standard RNA sequence."""

    # constructed by RNA moltype


@register_deserialiser(get_object_provenance(RnaSequence))
def deserialise_rna_sequence(data) -> RnaSequence:
    return RnaSequence.from_rich_dict(data)


class SliceRecordABC(ABC):
    """Abstract base class for recording the history of operations to be applied
    to some underlying data. Provides slicing functionality for the underlying data.

    Parameters
    ----------
    start
        start of the slice (inclusive indexing)
    stop
        stop of the slice (exclusive indexing)
    step
        step of the slice
    offset
        can be set with any additional offset that exists before the start of
        the underlying data
    parent_len
        length of the underlying data (not including offset)
    """

    __slots__ = ("start", "stop", "step", "_offset")

    @property
    @abstractmethod
    def parent_len(self) -> int: ...

    @abstractmethod
    def _get_init_kwargs(self) -> dict:
        """return required arguments for construction that are unique to the
        subclass"""
        ...

    @abstractmethod
    def copy(self): ...

    # refactor: design
    # can we remove the need for this method on the ABC and inheriting
    @property
    @abstractmethod
    def _zero_slice(self): ...

    @property
    def offset(self) -> int:
        return self._offset

    @offset.setter
    def offset(self, value: int):
        value = value or 0
        self._offset = int(value)

    @property
    def parent_start(self) -> int:
        """returns the start on the parent plus strand

        Returns
        -------
        offset + start, taking into account whether reversed. Result
        is positive.
        """
        if self.is_reversed:
            # self.stop becomes the start, self.stop will be negative
            assert self.stop < 0, "expected stop on reverse strand SeqView < 0"
            start = self.stop + self.parent_len + 1
        else:
            start = self.start

        return self.offset + start

    @property
    def is_reversed(self):
        return self.step < 0

    @property
    def parent_stop(self) -> int:
        """returns the stop on the parent plus strand

        Returns
        -------
        offset + stop, taking into account whether reversed. Result
        is positive.
        """
        if self.is_reversed:
            # self.start becomes the stop, self.start will be negative
            assert self.start < 0, "expected start on reverse strand SeqView < 0"
            stop = self.start + self.parent_len + 1
        else:
            stop = self.stop
        return self.offset + stop

    def absolute_position(self, rel_index: int, include_boundary: bool = False):
        """Converts an index relative to the current view to be with respect
        to the coordinates of the original "Python sequence".

        Parameters
        ----------
        rel_index
            relative position with respect to the current view

        Returns
        -------
        the absolute index with respect to the coordinates of the original
        sequence (including offset if present).
        """
        if not self:
            return 0

        if rel_index < 0:
            raise IndexError("only positive indexing supported!")

        # _get_index return the absolute position relative to the underlying sequence
        seq_index, _, _ = self._get_index(rel_index, include_boundary=include_boundary)

        offset = self.offset
        return (
            offset + self.parent_len + seq_index + 1
            if self.is_reversed
            else offset + seq_index
        )

    def relative_position(self, abs_index: int, stop: bool = False):
        """converts an index on the original "Python sequence" into an index
        on this "view"

        Notes
        -----
        The returned value DOES NOT reflect python indexing. Importantly,
        negative values represent positions that precede the current view.
        """
        if not self:
            return 0

        if abs_index < 0:
            raise IndexError("Index must be +ve and relative to the + strand")

        if self.is_reversed:
            offset = self.offset

            if (
                tmp := (self.parent_len - abs_index + offset + self.start + 1)
            ) % self.step == 0 or stop:
                rel_pos = tmp // abs(self.step)
            else:
                rel_pos = (tmp // abs(self.step)) + 1

        else:
            offset = self.offset + self.start

            if (tmp := abs_index - offset) % self.step == 0 or stop:
                rel_pos = tmp // self.step
            else:
                rel_pos = (tmp // self.step) + 1

        return rel_pos

    def __len__(self):
        return abs((self.start - self.stop) // self.step)

    def __getitem__(self, segment: typing.Union[int, slice]):
        kwargs = self._get_init_kwargs()

        if _is_int(segment):
            start, stop, step = self._get_index(segment)
            return self.__class__(
                start=start,
                stop=stop,
                step=step,
                offset=self.offset,
                parent_len=self.parent_len,
                **kwargs,
            )

        if segment.start is segment.stop is segment.step is None:
            return self.copy()

        if len(self) == 0:
            return self

        if segment.start is not None and segment.start == segment.stop:
            return self._zero_slice

        slice_step = 1 if segment.step is None else segment.step

        if slice_step > 0:
            return self._get_slice(segment, slice_step, **kwargs)
        elif slice_step < 0:
            return self._get_reverse_slice(segment, slice_step, **kwargs)
        else:
            raise ValueError(
                f"{self.__class__.__name__} cannot be sliced with a step of 0"
            )

    def _get_index(self, val: int, include_boundary: bool = False):
        if len(self) == 0:
            raise IndexError(val)

        if val > 0 and include_boundary and val > len(self):
            raise IndexError(val)
        elif val > 0 and not include_boundary and val >= len(self):
            raise IndexError(val)
        elif val < 0 and include_boundary and abs(val) > (len(self) + 1):
            raise IndexError(val)
        elif val < 0 and not include_boundary and abs(val) > len(self):
            raise IndexError(val)

        if self.step > 0:
            if val >= 0:
                val = self.start + val * self.step
            else:
                val = self.start + len(self) * self.step + val * abs(self.step)

            return val, val + 1, 1

        elif self.step < 0:
            if val >= 0:
                val = self.start + val * self.step
            else:
                val = self.start + len(self) * self.step + val * self.step

            return val, val - 1, -1

    def _get_slice(self, segment: slice, step: int, **kwargs):
        slice_start = segment.start if segment.start is not None else 0
        slice_stop = segment.stop if segment.stop is not None else len(self)

        if self.step > 0:
            return self._get_forward_slice_from_forward_seqview_(
                slice_start, slice_stop, step, **kwargs
            )

        elif self.step < 0:
            return self._get_forward_slice_from_reverse_seqview_(
                slice_start, slice_stop, step, **kwargs
            )

    def _get_forward_slice_from_forward_seqview_(
        self, slice_start: int, slice_stop: int, step: int, **kwargs
    ):
        start = (
            self.start + slice_start * self.step
            if slice_start >= 0
            else max(
                self.start + len(self) * self.step + slice_start * self.step,
                self.start,
            )
        )
        if slice_stop > self.stop:
            stop = self.stop
        elif slice_stop >= 0:
            stop = self.start + slice_stop * self.step
        else:
            # "true stop" adjust for if abs(stop-start) % step != 0
            # "true stop" = self.start + len(self) * self.step
            stop = self.start + len(self) * self.step + slice_stop * self.step

        # if -ve, it's an invalid slice
        if start < 0 or stop < 0:
            return self._zero_slice

        # checking for zero-length slice
        if stop < start:
            return self._zero_slice
        if start > self.parent_len:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=min(self.stop, stop),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_forward_slice_from_reverse_seqview_(
        self, slice_start: int, slice_stop: int, step: int, **kwargs
    ):
        if slice_start >= 0:
            start = self.start + slice_start * self.step
        elif abs(slice_start) > len(self):
            start = self.start
        else:
            start = self.start + len(self) * self.step + slice_start * self.step

        if slice_stop >= 0:
            stop = self.start + slice_stop * self.step
        else:  # slice_stop < 0
            stop = self.start + len(self) * self.step + slice_stop * self.step

        # if +ve, it's an invalid slice
        if start >= 0 or stop >= 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=max(self.stop, stop),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_reverse_slice(self, segment: slice, step: int, **kwargs):
        slice_start = segment.start if segment.start is not None else -1
        slice_stop = segment.stop if segment.stop is not None else -len(self) - 1

        if self.step < 0:
            return self._get_reverse_slice_from_reverse_seqview_(
                slice_start, slice_stop, step, **kwargs
            )
        elif self.step > 0:
            return self._get_reverse_slice_from_forward_seqview_(
                slice_start, slice_stop, step, **kwargs
            )

    def _get_reverse_slice_from_forward_seqview_(
        self, slice_start: int, slice_stop: int, step: int, **kwargs
    ):
        # "true stop" adjust for if abs(stop-start) % step != 0
        # max possible start is "true stop" - step, because stop is not inclusive
        # "true stop" - step is converted to -ve index via subtracting len(self)
        if slice_start >= len(self):
            start = (self.start + len(self) * self.step - self.step) - self.parent_len
        elif slice_start >= 0:
            start = (self.start + slice_start * self.step) - self.parent_len
        else:
            start = (
                self.start
                + len(self) * self.step
                + slice_start * self.step
                - self.parent_len
            )

        if slice_stop >= self.parent_len:
            return self._zero_slice

        if slice_stop >= 0:
            stop = self.start + (slice_stop * self.step) - self.parent_len
        else:
            stop = (
                self.start
                + (len(self) * self.step)
                + (slice_stop * self.step)
                - self.parent_len
            )

        if start >= 0 or stop >= 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=max(stop, self.start - self.parent_len - 1),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_reverse_slice_from_reverse_seqview_(
        self, slice_start: int, slice_stop: int, step: int, **kwargs
    ):
        # max start is "true stop" + abs(step), because stop is not inclusive
        # "true stop" adjust for if abs(stop-start) % step != 0
        if slice_start >= len(self):
            start = (
                self.parent_len + self.start + len(self) * self.step + abs(self.step)
            )
        elif slice_start >= 0:
            start = self.parent_len + (self.start + slice_start * self.step)
        else:
            start = self.parent_len + (
                self.start + len(self) * self.step + slice_start * self.step
            )

        if slice_stop >= 0:
            stop = self.parent_len + (self.start + slice_stop * self.step)
            if stop <= self.parent_len + self.stop:
                return self._zero_slice
        else:
            stop = self.parent_len + (
                self.start + len(self) * self.step + slice_stop * self.step
            )
            if stop > self.parent_len + self.start:
                stop = self.parent_len + self.start + 1

        # if -ve, it's an invalid slice becomes zero
        # checking for zero-length slice
        if stop < start or start > self.parent_len or min(start, stop) < 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=stop,
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )


class SliceRecord(SliceRecordABC):
    __slots__ = "_parent_len"

    def __init__(
        self,
        *,
        parent_len: int,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
        offset: int = 0,
    ):
        if step == 0:
            raise ValueError("step cannot be 0")
        step = step if step is not None else 1
        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(parent_len, start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self._parent_len = parent_len
        self._offset = offset or 0

    @property
    def parent_len(self) -> int:
        return self._parent_len

    def _get_init_kwargs(self) -> dict:
        return {}

    def copy(self, sliced: bool = False):
        # todo: kath
        return self

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(start={self.start}, stop={self.stop}, step={self.step}, "
            f"parent_len={self.parent_len}, offset={self.offset})"
        )

    @property
    def _zero_slice(self):
        return self.__class__(start=0, stop=0, step=1, parent_len=self.parent_len)

    def to_rich_dict(self) -> dict:
        data = {"type": get_object_provenance(self), "version": __version__}
        data["init_args"] = {
            "parent_len": self.parent_len,
            "start": self.start,
            "stop": self.stop,
            "step": self.step,
            "offset": self.offset,
        }
        return data


class SeqViewABC(ABC):
    """
    An abstract base class for providing a view of a sequence.

    This class defines an interface for sequence views, by which operations
    performed on the sequence do not modify the original sequence data. Instead,
    modifications and operations are stored on the returned view of the sequence
    and can be realised by accessing the values of the view.
    """

    __slots__ = ()

    @property
    @abstractmethod
    def seqid(self) -> OptStr: ...

    @property
    @abstractmethod
    def parent_len(self) -> int: ...

    @property
    @abstractmethod
    def slice_record(self) -> SliceRecordABC: ...

    @property
    def offset(self) -> int:
        return self.slice_record.offset

    @property
    def is_reversed(self) -> bool:
        return self.slice_record.is_reversed

    @property
    @abstractmethod
    def str_value(self) -> str: ...

    @property
    @abstractmethod
    def array_value(self) -> array: ...

    @property
    @abstractmethod
    def bytes_value(self) -> bytes: ...

    @abstractmethod
    def copy(self, sliced: bool = False): ...

    @abstractmethod
    def to_rich_dict(self) -> dict: ...

    @abstractmethod
    def __str__(self) -> str: ...

    @abstractmethod
    def __array__(self, dtype=None, copy=None): ...

    @abstractmethod
    def __bytes__(self): ...

    @abstractmethod
    def __getitem__(self, segment: typing.Union[int, slice]) -> SeqViewABC: ...

    def __len__(self):
        return len(self.slice_record)

    def with_offset(self, offset: int):
        if self._slice_record.offset:
            raise ValueError(
                f"cannot set {offset=} on a SeqView with an offset {self._slice_record.offset=}"
            )

        init_kwargs = self._get_init_kwargs()
        init_kwargs["offset"] = offset
        return self.__class__(**init_kwargs)


class SeqView(SeqViewABC):
    """
    Provides a view of a sequence with support for slicing operations.

    This class represents a view of a sequence, allowing for efficient slicing
    without altering the original sequence data.

    Parameters
    ----------
    parent
        the original sequence data
    alphabet
        the alphabet object defining valid characters for the sequence
    seqid
        the name or identifier of the sequence
    parent_len
        the length of the sequence. Defaults to the length of the input sequence

    """

    # todo: kath,
    # update the docstring to reflect the new design

    __slots__ = (
        "parent",
        "alphabet",
        "_seqid",
        "_parent_len",
        "_slice_record",
    )

    def __init__(
        self,
        *,
        parent: typing.Union[str, "SeqsDataABC"],
        alphabet: new_alphabet.AlphabetABC,
        seqid: OptStr = None,
        parent_len: OptInt = None,
        slice_record: SliceRecordABC = None,
        offset: int = 0,
    ):
        self.alphabet = alphabet
        self.parent = parent
        self._seqid = seqid
        self._parent_len = parent_len or len(self.parent)
        self._slice_record = (
            slice_record
            if slice_record is not None
            else SliceRecord(parent_len=self._parent_len)
        )
        if offset and self._slice_record.offset:
            raise ValueError(
                f"cannot set {offset=} on a SeqView with an offset {self._slice_record.offset=}"
            )
        elif offset:
            self._slice_record.offset = offset

    @property
    def seqid(self) -> str:
        return self._seqid

    @property
    def slice_record(self) -> SliceRecordABC:
        return self._slice_record

    @property
    def parent_len(self) -> int:
        return self._parent_len

    def _get_init_kwargs(self):
        return {
            "parent": self.parent,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    @property
    def str_value(self):
        return self.parent[
            self.slice_record.start : self.slice_record.stop : self.slice_record.step
        ]

    @property
    def array_value(self):
        return self.alphabet.to_indices(self.str_value)

    @property
    def bytes_value(self):
        return self.str_value.encode("utf-8")

    def __str__(self) -> str:
        return self.str_value

    def __array__(self, dtype=None, copy=None) -> numpy.ndarray:
        arr = self.array_value
        if dtype is not None:
            arr = arr.astype(dtype)
        return arr

    def __bytes__(self) -> bytes:
        return self.bytes_value

    def __getitem__(self, segment: typing.Union[int, slice]) -> SeqViewABC:
        return self.__class__(
            parent=self.parent,
            seqid=self.seqid,
            alphabet=self.alphabet,
            parent_len=self.parent_len,
            slice_record=self.slice_record[segment],
        )

    def __repr__(self) -> str:
        seq_preview = (
            f"{self.parent[:10]}...{self.parent[-5:]}"
            if self.parent_len > 15
            else self.parent
        )
        return (
            f"{self.__class__.__name__}(seqid={self.seqid!r}, parent={seq_preview!r}, "
            f"slice_record={self.slice_record.__repr__()})"
        )

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict

        Notes
        -----
        This method will slice the underlying sequence to the start and stop values

        Warning
        -------
        This method is not intended to provide serialisation of this object,
        instead, it is intended for usage by an enclosing class.
        """
        # get the current state
        data = {"type": get_object_provenance(self), "version": __version__}
        data["init_args"] = self._get_init_kwargs()

        if self.slice_record.is_reversed:
            adj = self.parent_len + 1
            start, stop = self.slice_record.stop + adj, self.slice_record.start + adj

        else:
            start, stop = self.slice_record.start, self.slice_record.stop

        data["init_args"]["parent"] = self.parent[start:stop]
        new_sr = SliceRecord(
            parent_len=(stop - start),
            step=self.slice_record.step,
        )
        data["init_args"]["slice_record"] = new_sr.to_rich_dict()
        data["init_args"]["alphabet"] = self.alphabet.to_rich_dict()
        return data

    def copy(self, sliced: bool = False):
        """returns copy

        Parameters
        ----------
        sliced
            if True, the underlying sequence is truncated and the start/stop
            adjusted
        """
        if not sliced:
            return self.__class__(
                parent=self.parent,
                seqid=self.seqid,
                alphabet=self.alphabet,
                slice_record=self.slice_record.copy(),
                parent_len=self.parent_len,
            )
        data = self.to_rich_dict()
        return self.__class__(
            parent=data["init_args"]["parent"],
            seqid=self.seqid,
            alphabet=self.alphabet,
            slice_record=SliceRecord(**data["init_args"]["slice_record"]["init_args"]),
        )


@singledispatch
def _coerce_to_seqview(data, seqid, alphabet, offset) -> SeqViewABC:
    from cogent3.core.alignment import Aligned
    from cogent3.core.sequence import Sequence as old_Sequence
    from cogent3.core.sequence import SeqView as old_SeqView

    if isinstance(data, (Aligned, old_Sequence, old_SeqView)):
        return _coerce_to_seqview(str(data), seqid, alphabet, offset)
    raise NotImplementedError(f"{type(data)}")


@_coerce_to_seqview.register
def _(data: SeqViewABC, seqid, alphabet, offset) -> SeqViewABC:
    # we require the indexes of shared states in alphabets to be the same
    # SeqView has an alphabet but SeqViewABC does NOT because that is
    # more general and covers the case where the SeqsData collection has the
    # alphabet
    if hasattr(data, "alphabet"):
        n = min(len(data.alphabet), len(alphabet))
        if data.alphabet[:n] != alphabet[:n]:
            raise new_alphabet.AlphabetError(
                f"element order {data.alphabet=} != to that in {alphabet=} for {data=!r}"
            )

    if offset and data.offset:
        raise ValueError(
            f"cannot set {offset=} on a SeqView with an offset {data.offset=}"
        )
    elif offset:
        return data.with_offset(offset)
    return data


@_coerce_to_seqview.register
def _(data: Sequence, seqid, alphabet, offset) -> SeqViewABC:
    return _coerce_to_seqview(data._seq, seqid, alphabet, offset)


@_coerce_to_seqview.register
def _(data: str, seqid, alphabet, offset) -> SeqViewABC:
    return SeqView(parent=data, seqid=seqid, alphabet=alphabet, offset=offset)


@_coerce_to_seqview.register
def _(data: bytes, seqid, alphabet, offset) -> SeqViewABC:
    data = data.decode("utf8")
    return SeqView(parent=data, seqid=seqid, alphabet=alphabet, offset=offset)


@_coerce_to_seqview.register
def _(data: tuple, seqid, alphabet, offset) -> SeqViewABC:
    return _coerce_to_seqview("".join(data), seqid, alphabet, offset)


@_coerce_to_seqview.register
def _(data: list, seqid, alphabet, offset) -> SeqViewABC:
    return _coerce_to_seqview("".join(data), seqid, alphabet, offset)

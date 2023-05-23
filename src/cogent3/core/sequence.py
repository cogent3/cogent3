"""Contains classes that represent biological sequence data. These
provide generic biological sequence manipulation functions, plus functions
that are critical for the EVOLVE calculations.

WARNING: Do not import sequence classes directly! It is expected that you will
access them through the moltype module. Sequence classes depend on information
from the MolType that is _only_ available after MolType has been imported.

Sequences are intended to be immutable. This is not enforced by the code for
performance reasons, but don't alter the MolType or the sequence data after
creation.
"""
from __future__ import annotations

import contextlib
import copy
import json
import os
import re
import warnings

from collections import defaultdict
from functools import singledispatch, total_ordering
from operator import eq, ne
from random import shuffle
from typing import Generator, List, Optional, Tuple

from numpy import (
    arange,
    array,
    compress,
    logical_not,
    logical_or,
    nonzero,
    put,
    ravel,
    take,
    zeros,
)
from numpy.random import permutation

from cogent3._version import __version__
from cogent3.core.alphabet import AlphabetError
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    BasicAnnotationDb,
    FeatureDataType,
    GenbankAnnotationDb,
    GffAnnotationDb,
    SupportsFeatures,
    load_annotations,
)
from cogent3.core.genetic_code import get_code
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import LostSpan, Map
from cogent3.format.fasta import alignment_to_fasta
from cogent3.maths.stats.contingency import CategoryCounts
from cogent3.maths.stats.number import CategoryCounter
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.util.misc import (
    DistanceFromMatrix,
    bytes_to_string,
    get_object_provenance,
    get_setting_from_environ,
)
from cogent3.util.transform import for_seq, per_shortest
from cogent3.util.warning import (
    deprecated,
    deprecated_args,
    deprecated_callable,
)


ARRAY_TYPE = type(array(1))

# standard distance functions: left  because generally useful
frac_same = for_seq(f=eq, aggregator=sum, normalizer=per_shortest)
frac_diff = for_seq(f=ne, aggregator=sum, normalizer=per_shortest)


@total_ordering
class SequenceI(object):
    """Abstract class containing Sequence interface.

    Specifies methods that Sequence delegates to its MolType, and methods for
    detecting gaps.
    """

    # String methods delegated to self._seq -- remember to override if self._seq
    # isn't a string in your base class, but it's probably better to make
    # self._seq a property that contains the string.
    line_wrap = None  # used for formatting FASTA strings

    def __str__(self):
        """__str__ returns self._seq unmodified."""
        return str(self._seq)

    def to_fasta(self, make_seqlabel=None, block_size=60):
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
        return alignment_to_fasta({label: str(self)}, block_size=block_size)

    def to_rich_dict(self, exclude_annotations=False):
        """returns {'name': name, 'seq': sequence, 'moltype': moltype.label}"""
        info = {} if self.info is None else self.info
        if not info.get("Refs", None) is None and "Refs" in info:
            info.pop("Refs")

        info = info or None
        data = dict(
            name=self.name,
            seq=str(self),
            moltype=self.moltype.label,
            info=info,
            type=get_object_provenance(self),
            version=__version__,
        )
        if hasattr(self, "annotation_offset"):
            offset = self._seq.offset + self._seq.start
            data.update(dict(annotation_offset=offset))

        if (
            hasattr(self, "annotation_db")
            and self.annotation_db
            and not exclude_annotations
        ):
            data["annotation_db"] = self.annotation_db.to_rich_dict()

        return data

    def to_json(self):
        """returns a json formatted string"""
        return json.dumps(self.to_rich_dict())

    def translate(self, *args, **kwargs):  # pragma: no cover
        """returns the result of call str.translate

        Notes
        -----
        This is a string method, nothing to do with translating into a
        protein sequence.
        """
        from cogent3.util.warning import discontinued

        discontinued(
            "function",
            "cogent3.core.sequence.Sequence.translate",
            "2023.6",
            "Better if user just converts sequence to string.",
        )
        return str(self).translate(*args, **kwargs)

    def count(self, item):
        """count() delegates to self._seq."""
        return str(self).count(item)

    def counts(
        self,
        motif_length=1,
        include_ambiguity=False,
        allow_gap=False,
        exclude_unobserved=False,
    ):
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

        """
        try:
            data = str(self)
        except AttributeError:
            data = self._data

        not_array = isinstance(data, (SeqView, str))

        if motif_length == 1:
            counts = CategoryCounter(data)
        else:
            if len(data) % motif_length != 0:
                warnings.warn(
                    f"{self.name} length not divisible by {motif_length}, truncating"
                )
            limit = (len(data) // motif_length) * motif_length
            data = data[:limit]
            if not_array:
                counts = CategoryCounter(
                    data[i : i + motif_length] for i in range(0, limit, motif_length)
                )
            else:
                counts = CategoryCounter(
                    tuple(v) for v in data.reshape(limit // motif_length, motif_length)
                )
        if not not_array:
            for key in list(counts):
                indices = [key] if motif_length == 1 else key
                motif = self.alphabet.to_chars(indices).astype(str)
                motif = "".join(motif)
                counts[motif] = counts.pop(key)

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
        return self.__class__("".join(randomized_copy_list), info=self.info)

    def complement(self):
        """Returns complement of self, using data from MolType.

        Always tries to return same type as item: if item looks like a dict,
        will return list of keys.
        """
        return self.__class__(self.moltype.complement(self), info=self.info)

    def strip_degenerate(self):
        """Removes degenerate bases by stripping them out of the sequence."""
        return self.__class__(self.moltype.strip_degenerate(self), info=self.info)

    def strip_bad(self):
        """Removes any symbols not in the alphabet."""
        return self.__class__(self.moltype.strip_bad(self), info=self.info)

    def strip_bad_and_gaps(self):
        """Removes any symbols not in the alphabet, and any gaps."""
        return self.__class__(self.moltype.strip_bad_and_gaps(self), info=self.info)

    def rc(self):
        """Returns reverse complement of self w/ data from MolType.

        Always returns same type self.
        """
        return self.__class__(self.moltype.rc(self), info=self.info)

    def is_gapped(self):
        """Returns True if sequence contains gaps."""
        return self.moltype.is_gapped(self)

    def is_gap(self, char=None):
        """Returns True if char is a gap.

        If char is not supplied, tests whether self is gaps only.
        """
        if char is None:  # no char - so test if self is all gaps
            return len(self) == self.count_gaps()
        else:
            return self.moltype.is_gap(char)

    def is_degenerate(self):
        """Returns True if sequence contains degenerate characters."""
        return self.moltype.is_degenerate(self)

    def is_valid(self):
        """Returns True if sequence contains no items absent from alphabet."""
        return self.moltype.is_valid(self)

    def is_strict(self):
        """Returns True if sequence contains only monomers."""
        return self.moltype.is_strict(self)

    def first_gap(self):
        """Returns the index of the first gap in the sequence, or None."""
        return self.moltype.first_gap(self)

    def first_degenerate(self):
        """Returns the index of first degenerate symbol in sequence, or None."""
        return self.moltype.first_degenerate(self)

    def first_invalid(self):
        """Returns the index of first invalid symbol in sequence, or None."""
        return self.moltype.first_invalid(self)

    def first_non_strict(self):
        """Returns the index of first non-strict symbol in sequence, or None."""
        return self.moltype.first_non_strict(self)

    def disambiguate(self, method="strip"):
        """Returns a non-degenerate sequence from a degenerate one.

        method can be 'strip' (deletes any characters not in monomers or gaps)
        or 'random'(assigns the possibilities at random, using equal
        frequencies).
        """
        return self.__class__(self.moltype.disambiguate(self, method), info=self.info)

    def degap(self):
        """Deletes all gap characters from sequence."""
        return self.__class__(self.moltype.degap(self), name=self.name, info=self.info)

    def gap_indices(self):
        """Returns list of indices of all gaps in the sequence, or []."""
        return self.moltype.gap_indices(self)

    def gap_vector(self):
        """Returns vector of True or False according to which pos are gaps."""
        return self.moltype.gap_vector(self)

    def gap_maps(self):
        """Returns dicts mapping between gapped and ungapped positions."""
        return self.moltype.gap_maps(self)

    def count_gaps(self):
        """Counts the gaps in the specified sequence."""
        return self.moltype.count_gaps(self)

    def count_degenerate(self):
        """Counts the degenerate bases in the specified sequence."""
        return self.moltype.count_degenerate(self)

    def possibilities(self):
        """Counts number of possible sequences matching the sequence.

        Uses self.degenerates to decide how many possibilities there are at
        each position in the sequence.
        """
        return self.moltype.possibilities(self)

    def mw(self, method="random", delta=None):
        """Returns the molecular weight of (one strand of) the sequence.

        If the sequence is ambiguous, uses method (random or strip) to
        disambiguate the sequence.

        If delta is passed in, adds delta per strand (default is None, which
        uses the alphabet default. Typically, this adds 18 Da for terminal
        water. However, note that the default nucleic acid weight assumes
        5' monophosphate and 3' OH: pass in delta=18.0 if you want 5' OH as
        well.

        Note that this method only calculates the MW of the coding strand. If
        you want the MW of the reverse strand, add self.rc().mw(). DO NOT
        just multiply the MW by 2: the results may not be accurate due to
        strand bias, e.g. in mitochondrial genomes.
        """
        return self.moltype.mw(self, method, delta)

    def can_match(self, other):
        """Returns True if every pos in self could match same pos in other.

        Truncates at length of shorter sequence.
        gaps are only allowed to match other gaps.
        """
        return self.moltype.can_match(self, other)

    def can_mismatch(self, other):
        """Returns True if any position in self could mismatch with other.

        Truncates at length of shorter sequence.
        gaps are always counted as matches.
        """
        return self.moltype.can_mismatch(self, other)

    def must_match(self, other):
        """Returns True if all positions in self must match positions in other."""
        return self.moltype.must_match(self, other)

    def can_pair(self, other):
        """Returns True if self and other could pair.

        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.
        gaps are only allowed to pair with other gaps, and are counted as 'weak'
        (same category as GU and degenerate pairs).

        NOTE: second must be able to be reverse
        """
        return self.moltype.can_pair(self, other)

    def can_mispair(self, other):
        """Returns True if any position in self could mispair with other.

        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.
        gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        return self.moltype.can_mispair(self, other)

    def must_pair(self, other):
        """Returns True if all positions in self must pair with other.

        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.
        """
        return not self.moltype.can_mispair(self, other)

    def diff(self, other):
        """Returns number of differences between self and other.

        NOTE: truncates at the length of the shorter sequence. Case-sensitive.
        """
        return self.distance(other)

    def distance(self, other, function=None):
        """Returns distance between self and other using function(i,j).

        other must be a sequence.

        function should be a function that takes two items and returns a
        number. To turn a 2D matrix into a function, use
        cogent3.util.miscs.DistanceFromMatrix(matrix).

        NOTE: Truncates at the length of the shorter sequence.

        Note that the function acts on two _elements_ of the sequences, not
        the two sequences themselves (i.e. the behavior will be the same for
        every position in the sequences, such as identity scoring or a function
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

    def matrix_distance(self, other, matrix):
        """Returns distance between self and other using a score matrix.

        WARNING: the matrix must explicitly contain scores for the case where
        a position is the same in self and other (e.g. for a distance matrix,
        an identity between U and U might have a score of 0). The reason the
        scores for the 'diagonals' need to be passed explicitly is that for
        some kinds of distance matrices, e.g. log-odds matrices, the 'diagonal'
        scores differ from each other. If these elements are missing, this
        function will raise a KeyError at the first position that the two
        sequences are identical.
        """
        return self.distance(other, DistanceFromMatrix(matrix))

    def frac_same(self, other):
        """Returns fraction of positions where self and other are the same.

        Truncates at length of shorter sequence.
        Note that frac_same and frac_diff are both 0 if one sequence is empty.
        """
        return frac_same(self, other)

    def frac_diff(self, other):
        """Returns fraction of positions where self and other differ.

        Truncates at length of shorter sequence.
        Note that frac_same and frac_diff are both 0 if one sequence is empty.
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

    def frac_similar(self, other, similar_pairs):
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
        if first_nongap is None:  # sequence was all gaps
            result = self.__class__([missing for _ in len(self)], info=self.info)
        else:
            prefix = missing * first_nongap
            mid = str(self[first_nongap : last_nongap + 1])
            suffix = missing * (len(self) - last_nongap - 1)
            result = self.__class__(prefix + mid + suffix, info=self.info)
        return result

    def replace(self, oldchar, newchar):
        """return new instance with oldchar replaced by newchar"""
        NotImplemented

    def strand_symmetry(self, *args, **kwargs):
        raise TypeError("must be DNA or RNA moltype")

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
        wrap=60,
        limit=None,
        colors=None,
        font_size=12,
        font_family="Lucida Console",
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
            {character
            moltype.
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
            summary = f"{len(self)} (truncated to {limit if limit else len(self)}) {class_name}"
        else:
            summary = f"{len(self)} {class_name}"

        text = [
            "<style>",
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

        # If two sequences with the same name are being added together the name should not be None
        if type(other) == type(self):
            name = self.name if self.name == other.name else None
        else:
            name = None

        return self.__class__(str(self) + other_seq, name=name)


@total_ordering
class Sequence(SequenceI):
    """Holds the standard Sequence object. Immutable."""

    moltype = None  # connected to ACSII when moltype is imported

    def __init__(
        self,
        seq="",
        name=None,
        info=None,
        check=True,
        preserve_case=False,
        gaps_allowed=True,
        wildcards_allowed=True,
        annotation_offset=0,
    ):
        """Initialize a sequence.

        Parameters
        ----------
        seq
            the raw sequence string, default is ''
        name
            the sequence name
        check
            if True (the default), validates against the MolType
        annotation_offset
            integer indicating start position relative to annotations
        """

        if name is None and hasattr(seq, "name"):
            name = seq.name
        self.name = name
        orig_seq = seq

        checker = (
            (lambda x: self.moltype.verify_sequence(x, gaps_allowed, wildcards_allowed))
            if check
            else (lambda x: x)
        )

        self._seq = _coerce_seq(seq, preserve_case, checker)

        if not isinstance(info, InfoClass):
            try:
                info = InfoClass(info)
            except TypeError:
                info = InfoClass()
        if hasattr(orig_seq, "info"):
            with contextlib.suppress(Exception):
                info.update(orig_seq.info)
        self.info = info

        self._repr_policy = dict(num_pos=60)
        self._annotation_db = BasicAnnotationDb()
        self.annotation_offset = annotation_offset

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

        return self._seq.offset

    @annotation_offset.setter
    def annotation_offset(self, value):
        self._seq.offset = value

    @property
    def annotation_db(self):
        return self._annotation_db

    @annotation_db.setter
    def annotation_db(self, value):
        if value == self._annotation_db:
            return

        if value and not isinstance(value, SupportsFeatures):
            raise TypeError(f"{type(value)} does not satisfy SupportsFeatures")

        # Without knowing the contents of the db we cannot
        # establish whether self.moltype is compatible, so
        # we rely on the user to get that correct
        # one approach to support validation might be to add
        # to the SupportsFeatures protocol a is_nucleic flag,
        # for both DNA and RNA. But if a user trys get_slice()
        # on a '-' strand feature, they will get a TypError.
        # I think that's enough.
        if self._annotation_db and len(self._annotation_db):
            # we take the union
            self._annotation_db = self._annotation_db.union(value)
        else:
            self._annotation_db = value

    @deprecated_args(
        "2023.7",
        reason="consistent nomenclature",
        old_new=[("feature_type", "biotype")],
    )
    def get_features(
        self,
        *,
        biotype: Optional[str] = None,
        name: Optional[str] = None,
        start: Optional[int] = None,
        stop: Optional[int] = None,
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
            start, end positions to search between, relative to offset
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

        # sort start / stop
        start, stop = (start, stop) if start < stop else (stop, start)

        # note: offset is handled by absolute_position
        # we set include_boundary=False because start is inclusive indexing,
        # i,e., the start cannot be equal to the length of the view
        query_start = self._seq.absolute_position(start, include_boundary=False)
        # we set include_boundary=True because stop is exclusive indexing,
        # i,e., the stop can be equal to the length of the view
        query_end = self._seq.absolute_position(stop, include_boundary=True)

        reversed = query_end < query_start
        if reversed:
            query_start, query_end = query_end, query_start

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
            seqid=self.name,
            name=name,
            biotype=biotype,
            start=query_start,
            end=query_end,  # todo gah end should be stop
            allow_partial=allow_partial,
        ):
            # spans need to be converted from absolute to relative positions
            # DO NOT do adjustment in make_feature since that's user facing,
            # and we expect them to make a feature manually
            spans = array(feature.pop("spans"), dtype=int)
            for i, v in enumerate(spans.ravel()):
                rel_pos = self._seq.relative_position(v)
                spans.ravel()[i] = rel_pos

            if reversed:
                # see above comment
                spans = len(self) - spans

            feature["spans"] = spans.tolist()
            yield self.make_feature(feature)

    def _relative_spans(self, spans):
        r_spans = []

        for s, e in spans:
            # reverse feature determined by absolute position
            reverse = s > e

            # todo: kath, think about logic of stop=true for reverse locations
            # I think, we need to know the whether it is reverse so that the adjustment
            # for step (given stop=True) is in the correct direction.
            # However, we need to keep s > e for a reverse span because this information
            # is used in the Map?
            if reverse:
                start = self._seq.relative_position(s, stop=True)
                end = self._seq.relative_position(e)
            else:
                start = self._seq.relative_position(s)
                end = self._seq.relative_position(e, stop=True)

            r_spans += [(start, end)]

        return r_spans

    @deprecated_callable(
        "2023.7", reason="simpler name", new="<instance>.get_features()"
    )
    def get_features_matching(self, **kwargs):
        """use .get_features()"""
        return self.get_features(**kwargs)

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
        seq_rced = self._seq.reverse
        # todo gah check consistency of relationship between reversed and strand
        # i.e. which object has responsibility for transforming the strand value
        # (a string) into a bool?
        revd = feature.pop("reversed", None) or feature.pop("strand", None) == "-"
        spans = feature.pop("spans", None)

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

            if coord[0] == coord[1]:
                # would really to adjust pre and post to be up to the next span
                continue
            new_spans.append(new.tolist())

        fmap = Map(locations=new_spans, parent_length=len(self))
        if pre or post:
            # create a lost span to represent the segment missing from
            # the instance
            spans = fmap.spans
            if pre:
                spans.insert(0, LostSpan(pre))
            if post:
                spans.append(LostSpan(post))
            fmap = Map(spans=spans, parent_length=len(self))

        if revd and not seq_rced:
            # the sequence is on the plus strand, and the
            # feature coordinates are also for the plus strand
            # but their order needs to be changed to indicate
            # reverse complement is required
            fmap = fmap.reversed()
        elif seq_rced and not revd:
            # plus strand feature, but the sequence reverse complemented
            # so we need to nucleic-reverse the map
            fmap = fmap.nucleic_reversed()
        elif seq_rced:
            # sequence is rc'ed, the feature was minus strand of
            # original, so needs to be both nucleic reversed and
            # then reversed
            fmap = fmap.nucleic_reversed().reversed()

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
        if isinstance(self.annotation_db, GffAnnotationDb):
            # the db is updated directly
            self.annotation_db = load_annotations(
                path=f, seqids=self.name, db=self.annotation_db.db
            )
        elif isinstance(self.annotation_db, GenbankAnnotationDb):
            raise ValueError("GenbankAnnotationDb already attached")
        else:
            db = load_annotations(path=f, seqids=self.name)
            self.annotation_db = (
                db.update(self.annotation_db) if len(self.annotation_db) else db
            )

        if offset:
            self.annotation_offset = offset

    def add_feature(
        self,
        *,
        biotype: str,
        name: str,
        spans: List[Tuple[int, int]],
        parent_id: Optional[str] = None,
        strand: Optional[str] = None,
        on_alignment: bool = False,
        seqid: Optional[str] = None,
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

    def to_moltype(self, moltype):
        """returns copy of self with moltype seq

        Parameters
        ----------
        moltype : str
            molecular type
        """
        from cogent3 import get_moltype

        if not moltype:
            raise ValueError(f"unknown moltype '{moltype}'")

        moltype = get_moltype(moltype)
        s = moltype.coerce_str(self._seq.value)
        moltype.verify_sequence(s, gaps_allowed=True, wildcards_allowed=True)
        sv = SeqView(s)
        new = moltype.make_seq(sv, name=self.name, info=self.info)
        new.annotation_db = self.annotation_db
        return new

    def _seq_filter(self, seq):
        """Returns filtered seq; used to do DNA/RNA conversions."""
        return seq

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

        if self.annotation_db and not isinstance(seq_db, type(db)):
            raise TypeError(f"type {type(seq_db)} != {type(db)}")

        if not self.annotation_db:
            self.annotation_db = type(seq_db)()

        self.annotation_db.union(seq_db, seqids=self.name)

    def copy(self, exclude_annotations=False):
        """returns a copy of self"""
        offset = self._seq.offset + self._seq.start
        new = self.__class__(
            self._seq[:], name=self.name, info=self.info, annotation_offset=offset
        )
        db = None if exclude_annotations else copy.deepcopy(self.annotation_db)
        new._annotation_db = db
        return new

    def with_masked_annotations(
        self, annot_types, mask_char=None, shadow=False, extend_query=False
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
            ambigs = [(len(v), c) for c, v in list(self.moltype.ambiguities.items())]
            ambigs.sort()
            mask_char = ambigs[-1][1]
        assert mask_char in self.moltype, f"Invalid mask_char {mask_char}"

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
                map=Map(locations=[], parent_length=len(self)),
            )
        else:
            region = annotations[0].union(annotations[1:])

        if shadow:
            region = region.shadow()

        i = 0
        segments = []
        fmap = region.map.reversed() if region.map.reverse else region.map
        for b, e in fmap.get_coordinates():
            segments.extend((str(self[i:b]), mask_char * (e - b)))
            i = e
        segments.append(str(self[i:]))

        new = self.__class__(
            "".join(segments), name=self.name, check=False, info=self.info
        )
        new.annotation_db = copy.deepcopy(self.annotation_db)
        return new

    def gapped_by_map_segment_iter(self, map, allow_gaps=True, recode_gaps=False):
        for span in map.spans:
            if span.lost:
                if allow_gaps:
                    unknown = span.terminal or recode_gaps
                    seg = "-?"[unknown] * span.length
                else:
                    raise ValueError(f"gap(s) in map {map}")
            else:
                seg = str(self[span.start : span.end])
                if span.reverse:
                    complement = self.moltype.complement
                    seg = [complement(base) for base in seg[::-1]]
                    seg = "".join(seg)
            yield seg

    def gapped_by_map_motif_iter(self, map):
        for segment in self.gapped_by_map_segment_iter(map):
            yield from segment

    def gapped_by_map(self, map, recode_gaps=False):
        # todo gah do we propagate annotations here?
        segments = self.gapped_by_map_segment_iter(map, True, recode_gaps)
        new = self.__class__(
            "".join(segments), name=self.name, check=False, info=self.info
        )
        return new

    def _mapped(self, map):
        # Called by generic __getitem__
        segments = self.gapped_by_map_segment_iter(map, allow_gaps=False)
        return self.__class__("".join(segments), self.name, info=self.info)

    def __repr__(self):
        myclass = f"{self.__class__.__name__}"
        myclass = myclass.split(".")[-1]
        if len(self) > 10:
            seq = f"{str(self[:7])}... {len(self)}"
        else:
            seq = str(self)
        return f"{myclass}({seq})"

    def __getitem__(self, index):
        preserve_offset = False
        if hasattr(index, "map"):
            index = index.map

        if isinstance(index, Map):
            new = self._mapped(index)
            preserve_offset = not index.reverse

        elif isinstance(index, (int, slice)):
            new = self.__class__(
                self._seq[index], name=self.name, check=False, info=self.info
            )
            stride = getattr(index, "step", 1) or 1
            preserve_offset = stride > 0

        if self.annotation_db is not None and preserve_offset:
            new.annotation_db = self.annotation_db
            new.annotation_offset = self.annotation_offset

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

    def __str__(self):
        result = str(self._seq)
        if self._seq.reverse:
            with contextlib.suppress(TypeError):
                result = self.moltype.complement(result)
        return result

    def gettype(self):  # pragma: no cover
        """Return the sequence type."""
        deprecated("method", "gettype", "get_type", "2023.6", "pep8")
        return self.get_type()

    def get_type(self):
        """Return the sequence type as moltype label."""
        return self.moltype.label

    def resolveambiguities(self):  # pragma: no cover
        """Returns a list of tuples of strings."""
        deprecated(
            "method",
            "resolveambiguities",
            "resolved_ambiguities",
            "2023.6",
            "pep8",
        )
        return self.resolved_ambiguities()

    def resolved_ambiguities(self) -> List[Tuple[str]]:
        """Returns a list of tuples of strings."""
        ambigs = self.moltype.ambiguities
        return [ambigs[motif] for motif in self._seq]

    def iter_kmers(self, k: int, strict: bool = False) -> Generator[str, None, None]:
        """generates all overlapping k-mers.
        When strict is True, the characters in the k-mer must be
        a subset of the canonical characters for the moltype"""
        if k <= 0:
            raise ValueError(f"k must be an int > 0, not {k}")

        if not isinstance(k, int):
            raise ValueError(f"k must be an int, not {k}")

        canonical = set(self.moltype)
        for i in range(len(self) - k + 1):
            kmer = self[i : i + k]
            if not strict:
                yield kmer
            else:
                if all(char in canonical for char in kmer):
                    yield kmer

    def get_kmers(self, k: int, strict: bool = False) -> List[str]:
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

    def get_in_motif_size(self, motif_length=1, log_warnings=True):
        """returns sequence as list of non-overlapping motifs

        Parameters
        ----------
        motif_length
            length of the motifs
        log_warnings
            whether to notify of an incomplete terminal motif

        """
        seq = self._seq
        if isinstance(seq, SeqView):
            seq = str(self)
        if motif_length == 1:
            return seq

        length = len(seq)
        remainder = length % motif_length
        if remainder and log_warnings:
            warnings.warn(
                f'Dropped remainder "{seq[-remainder:]}" from end of sequence'
            )
        return [
            seq[i : i + motif_length]
            for i in range(0, length - remainder, motif_length)
        ]

    def parse_out_gaps(self):
        """returns Map corresponding to gap locations and ungapped Sequence"""
        gapless = []
        segments = []
        nongap = re.compile(f"([^{re.escape('-')}]+)")
        for match in nongap.finditer(str(self)):
            segments.append(match.span())
            gapless.append(match.group())
        map = Map(segments, parent_length=len(self)).inverse()
        seq = self.__class__(
            "".join(gapless), name=self.get_name(), info=self.info, preserve_case=True
        )
        seq.annotation_db = copy.deepcopy(self.annotation_db)
        return map, seq

    def replace(self, oldchar, newchar):
        """return new instance with oldchar replaced by newchar"""
        new = self._seq.replace(oldchar, newchar)
        result = self.__class__(new, name=self.name, info=self.info)
        result.annotation_db = copy.deepcopy(self.annotation_db)
        return result

    def is_annotated(self):
        """returns True if sequence has any annotations"""
        num = self.annotation_db.num_matches() if self.annotation_db else 0
        return num != 0

    def annotate_matches_to(self, pattern, annot_type, name, allow_multiple=False):
        """Adds an annotation at sequence positions matching pattern.

        Parameters
        ----------
        pattern : string
            The search string for which annotations are made. IUPAC ambiguities
            are converted to regex on sequences with the appropriate MolType.
        annot_type : string
            The type of the annotation (e.g. "domain").
        name : string
            The name of the annotation.
        allow_multiple : boolean
            If True, allows multiple occurrences of the input pattern. Otherwise
            only the first match is used.

        Returns
        -------
        Returns a list of Feature instances.
        """
        try:
            pattern = self.moltype.to_regex(seq=pattern)
        except ValueError:
            # assume already a regex
            pass

        pos = [m.span() for m in re.finditer(pattern, str(self))]
        if not pos:
            return []

        num_match = len(pos) if allow_multiple else 1
        return [
            self.add_feature(
                biotype=annot_type,
                name=f"{name}:{i}" if allow_multiple else name,
                spans=[pos[i]],
            )
            for i in range(num_match)
        ]

    def __add__(self, other):
        """Adds two sequences (other can be a string as well)"""
        new_seq = super(Sequence, self).__add__(other)
        new_seq.annotation_db = None
        # annotations are dropped in this case
        return new_seq

    def get_drawables(self):
        """returns a dict of drawables, keyed by type"""
        result = defaultdict(list)
        for f in self.get_features():
            result[f.biotype].append(f.get_drawable())
        return result

    def get_drawable(self, width=600, vertical=False):
        """returns Drawable instance"""
        from cogent3.draw.drawable import Drawable

        drawables = self.get_drawables()
        if not drawables:
            return None
        # we order by tracks
        top = 0
        space = 0.25
        annotes = []
        for feature_type in drawables:
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


class ProteinSequence(Sequence):
    """Holds the standard Protein sequence."""

    pass


class ProteinWithStopSequence(Sequence):
    """Holds the standard Protein sequence, allows for stop codon."""

    pass


class NucleicAcidSequence(Sequence):
    """Abstract base class for DNA and RNA sequences."""

    codon_alphabet = None  # will set in moltype

    def reverse_complement(self):
        """Converts a nucleic acid sequence to its reverse complement.
        Synonymn for rc."""
        return self.rc()

    def rc(self):
        """Converts a nucleic acid sequence to its reverse complement."""
        rc = self.__class__(
            self._seq[::-1], name=self.name, check=False, info=self.info
        )
        rc.annotation_db = copy.deepcopy(self.annotation_db)
        return rc

    def has_terminal_stop(self, gc=None, allow_partial=False):
        """Return True if the sequence has a terminal stop codon.

        Parameters
        ----------
        gc
            genetic code object
        allow_partial
            if True and the sequence length is not dividisble
            by 3, ignores the 3' terminal incomplete codon

        """
        gc = get_code(gc)
        codons = self._seq
        divisible_by_3 = len(codons) % 3 == 0
        end3 = self.__class__(self._seq[-3:]).degap()
        if not allow_partial and (not divisible_by_3 or len(end3) != 3):
            raise ValueError("seq length not divisible by 3")

        if not divisible_by_3:
            return False

        return codons and gc.is_stop(codons[-3:])

    def trim_stop_codon(self, gc=None, allow_partial=False):
        """Removes a terminal stop codon from the sequence

        Parameters
        ----------
        gc
            genetic code object
        allow_partial
            if True and the sequence length is not divisible
            by 3, ignores the 3' terminal incomplete codon

        """
        gc = get_code(gc)
        codons = self._seq
        divisible_by_3 = len(codons) % 3 == 0

        if not allow_partial and not divisible_by_3:
            raise ValueError("seq length not divisible by 3")

        if divisible_by_3 and codons and gc.is_stop(str(codons[-3:])):
            codons = codons[:-3]

        result = self.__class__(codons, name=self.name, info=self.info)
        result.annotation_db = copy.deepcopy(self.annotation_db)
        return result

    def get_translation(self, gc=None, incomplete_ok=False, include_stop=False):
        """translate to amino acid sequence

        Parameters
        ----------
        gc
            name or ID of genetic code
        incomplete_ok
            codons that are mixes of nucleotide and gaps converted to '?'.
            raises a ValueError if False
        include_stop
            allows stop codons in translation

        Returns
        -------
        sequence of PROTEIN moltype

        Raises
        ------
        AlphabetError if include_stop is False and a stop codon occurs
        """
        from cogent3.core.moltype import get_moltype

        protein = get_moltype("protein_with_stop" if include_stop else "protein")
        gc = get_code(gc)
        codon_alphabet = gc.get_alphabet(include_stop=include_stop).with_gap_motif()
        # translate the codons
        translation = []
        for posn in range(0, len(self) - 2, 3):
            orig_codon = str(self[posn : posn + 3])
            try:
                resolved = codon_alphabet.resolve_ambiguity(orig_codon)
            except AlphabetError:
                if not incomplete_ok or "-" not in orig_codon:
                    raise
                resolved = (orig_codon,)
            trans = []
            for codon in resolved:
                if codon == "---":
                    aa = "-"
                elif "-" in codon:
                    aa = "?"
                    if not incomplete_ok:
                        raise AlphabetError(f"incomplete codon {codon} in {self.name}")
                else:
                    aa = gc[codon]
                    if aa == "*" and not include_stop:
                        continue
                trans.append(aa)
            if not trans:
                raise ValueError(orig_codon)
            aa = protein.what_ambiguity(trans)
            translation.append(aa)

        return protein.make_seq(seq="".join(translation), name=self.name)

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


def _input_vals_pos_step(seqlen, start, stop, step):
    start = 0 if start is None else start
    if start > 0 and start >= seqlen:
        # start beyond seq is an empty slice
        return 0, 0, 1

    stop = seqlen if stop is None else stop
    if stop < 0 and abs(stop) >= seqlen:
        # finished slice before we started seq!
        return 0, 0, 1

    start = max(seqlen + start, 0) if start < 0 else start

    if stop > 0:
        stop = min(seqlen, stop)
    elif stop < 0:
        stop += seqlen

    if start >= stop:
        start = stop = 0
        step = 1

    return start, stop, step


def _input_vals_neg_step(seqlen, start, stop, step):
    # Note how Python reverse slicing works
    # we need to make sure the start and stop are both
    # negative, for example "abcd"[-1:-5:-1] returns "dcba"
    if start is None or start >= seqlen:  # set default
        start = -1  # Done
    elif start >= 0:  # convert to -ve index
        start = start - seqlen
    elif start < -seqlen:  # start is bounded by len(seq)
        return 0, 0, 1

    if stop is None:  # set default
        stop = -seqlen - 1
    elif stop >= 0:
        stop -= seqlen

    stop = max(stop, -seqlen - 1)  # stop should always be <= len(seq)

    # checking for zero-length slice
    return (0, 0, 1) if start < stop else (start, stop, step)


class SeqView:
    __slots__ = ("seq", "start", "stop", "step", "_offset")

    def __init__(self, seq, *, start: int = None, stop: int = None, step: int = None):
        if step == 0:
            raise ValueError("step cannot be 0")
        step = 1 if step is None else step

        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(len(seq), start, stop, step)
        self.seq = seq
        self.start = start
        self.stop = stop
        self._offset = 0
        self.step = step

    @property
    def offset(self) -> int:
        return self._offset

    @offset.setter
    def offset(self, value: int):
        self._offset = value or 0

    @property
    def reverse(self):
        return self.step < 0

    def absolute_position(self, rel_index: int, include_boundary=False):
        """Converts an index relative to the current view to be with respect to the coordinates of the sequence's annotations

        Parameters
        ----------
        rel_index    int
                relative position with respect to the current view

        Returns
        -------
        the absolute index with respect to the coordinates of the sequence's annotations

        """

        if rel_index < 0:
            raise IndexError("only positive indexing supported!")

        # _get_index return the absolute position relative to the underlying sequence
        seq_index, _, _ = self._get_index(rel_index, include_boundary=include_boundary)

        # add offset and handle reversed views, now absolute relative to annotation coordinates
        offset = self.offset
        if self.reverse:
            abs_index = offset + len(self.seq) + seq_index + 1
        else:
            abs_index = offset + seq_index

        return abs_index

    def relative_position(self, abs_index, stop=False):
        """converts an index relative to annotation coordinates to be with respect to the current sequence view

        NOTE: the returned value DOES NOT reflect python indexing. Importantly, negative values represent positions that
        precede the current view.
        """

        if abs_index < 0:
            raise IndexError("Index must be +ve and relative to the + strand")

        if self.reverse:
            offset = self.offset

            if (
                tmp := ((len(self.seq) - (abs_index - offset)) + self.start + 1)
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

    def _get_index(self, val, include_boundary=False):
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

    def _get_slice(self, segment, step):
        slice_start = segment.start if segment.start is not None else 0
        slice_stop = segment.stop if segment.stop is not None else len(self)

        if self.step > 0:
            return self._get_forward_slice_from_forward_seqview_(
                slice_start, slice_stop, step
            )

        elif self.step < 0:
            return self._get_forward_slice_from_reverse_seqview_(
                slice_start, slice_stop, step
            )

    def _get_forward_slice_from_forward_seqview_(self, slice_start, slice_stop, step):
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
            return _zero_slice

        # checking for zero-length slice
        if stop < start:
            return _zero_slice
        if start > len(self.seq):
            return _zero_slice

        return self.__class__(
            self.seq,
            start=start,
            stop=min(self.stop, stop),
            step=self.step * step,
        )

    def _get_forward_slice_from_reverse_seqview_(self, slice_start, slice_stop, step):
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
            return _zero_slice

        return self.__class__(
            self.seq,
            start=start,
            stop=max(self.stop, stop),
            step=self.step * step,
        )

    def _get_reverse_slice(self, segment, step):
        slice_start = segment.start if segment.start is not None else -1
        slice_stop = segment.stop if segment.stop is not None else -len(self) - 1

        if self.step < 0:
            return self._get_reverse_slice_from_reverse_seqview_(
                slice_start, slice_stop, step
            )
        elif self.step > 0:
            return self._get_reverse_slice_from_forward_seqview_(
                slice_start, slice_stop, step
            )

    def _get_reverse_slice_from_forward_seqview_(self, slice_start, slice_stop, step):
        # "true stop" adjust for if abs(stop-start) % step != 0
        # max possible start is "true stop" - step, because stop is not inclusive
        # "true stop" - step is converted to -ve index via subtracting len(self)
        seq_len = len(self.seq)
        if slice_start >= len(self):
            start = (self.start + len(self) * self.step - self.step) - seq_len
        elif slice_start >= 0:
            start = (self.start + slice_start * self.step) - seq_len
        else:
            start = (
                self.start
                + len(self) * self.step
                + slice_start * self.step
                - len(self.seq)
            )

        if slice_stop >= seq_len:
            return _zero_slice

        if slice_stop >= 0:
            stop = self.start + (slice_stop * self.step) - seq_len
        else:
            stop = (
                self.start
                + (len(self) * self.step)
                + (slice_stop * self.step)
                - seq_len
            )

        if start >= 0 or stop >= 0:
            return _zero_slice

        return self.__class__(
            self.seq,
            start=start,
            stop=max(stop, self.start - len(self.seq) - 1),
            step=self.step * step,
        )

    def _get_reverse_slice_from_reverse_seqview_(self, slice_start, slice_stop, step):
        # max start is "true stop" + abs(step), because stop is not inclusive
        # "true stop" adjust for if abs(stop-start) % step != 0
        seq_len = len(self.seq)
        if slice_start >= len(self):
            start = seq_len + self.start + len(self) * self.step + abs(self.step)
        elif slice_start >= 0:
            start = seq_len + (self.start + slice_start * self.step)
        else:
            start = seq_len + (
                self.start + len(self) * self.step + slice_start * self.step
            )

        if slice_stop >= 0:
            stop = seq_len + (self.start + slice_stop * self.step)
            if stop <= seq_len + self.stop:
                return _zero_slice
        else:
            stop = seq_len + (
                self.start + len(self) * self.step + slice_stop * self.step
            )
            if stop > seq_len + self.start:
                stop = seq_len + self.start + 1

        # if -ve, it's an invalid slice becomes zero
        # checking for zero-length slice
        if stop < start or start > seq_len or min(start, stop) < 0:
            return _zero_slice

        return self.__class__(
            self.seq,
            start=start,
            stop=stop,
            step=self.step * step,
        )

    def __getitem__(self, segment):
        if isinstance(segment, int):
            start, stop, step = self._get_index(segment)
            return self.__class__(self.seq, start=start, stop=stop, step=step)

        if len(self) == 0:
            return self

        slice_step = 1 if segment.step is None else segment.step

        if slice_step > 0:
            return self._get_slice(segment, slice_step)
        elif slice_step < 0:
            return self._get_reverse_slice(segment, slice_step)
        else:
            raise ValueError(
                f"{self.__class__.__name__} cannot be sliced with a step of 0"
            )

    @property
    def value(self):
        return self.seq[self.start : self.stop : self.step]

    def replace(self, old, new):
        new_seq = self.seq.replace(old, new)
        if len(old) == len(new):
            return self.__class__(
                new_seq, start=self.start, stop=self.stop, step=self.step
            )

        return self.__class__(new_seq)

    def __len__(self):
        return abs((self.start - self.stop) // self.step)

    def __iter__(self):
        return iter(self.value)

    def __str__(self) -> str:
        return self.value

    def __repr__(self) -> str:
        seq = f"{self.seq[:10]}...{self.seq[-5:]}" if len(self.seq) > 15 else self.seq
        return f"{self.__class__.__name__}(seq={seq!r}, start={self.start}, stop={self.stop}, step={self.step})"


_zero_slice = SeqView(seq="")


class DnaSequence(NucleicAcidSequence):
    """Holds the standard DNA sequence."""

    def _seq_filter(self, seq):
        """Converts U to T."""
        return seq.replace("u", "t").replace("U", "T")


class RnaSequence(NucleicAcidSequence):
    """Holds the standard RNA sequence."""

    def _seq_filter(self, seq):
        """Converts T to U."""
        return seq.replace("t", "u").replace("T", "U")


class ABSequence(Sequence):
    """Holds a two-state sequence, with characters of 'a', 'b'"""

    pass


class ByteSequence(Sequence):
    """Used for storing arbitrary bytes."""

    def __init__(
        self, seq="", name=None, info=None, check=False, preserve_case=True, **kwargs
    ):
        super().__init__(
            seq=seq,
            name=name,
            info=info,
            check=check,
            preserve_case=preserve_case,
            **kwargs,
        )


@total_ordering
class ArraySequenceBase(object):
    """Holds the information for a non-degenerate sequence. Mutable.

    A ArraySequence is an array of indices of symbols, where those symbols are
    defined by an Alphabet. This representation of Sequence is convenient for
    counting symbol frequencies or tuple frequencies, remapping data (e.g. for
    reverse-complement), looking up model parameters, etc. Its main drawback is
    that the sequences can no longer be treated as strings, and conversion
    to/from strings can be fairly time-consuming. Also, any symbol not in the
    Alphabet cannot be represented at all.

    A sequence can have a name, which will be used for output in formats
    such as FASTA.

    A sequence Class has an alphabet (which can be overridden in instances
    where necessary), a delimiter used for string conversions, a line_wrap
    for wrapping characters into lines for e.g. FASTA output.

    Note that a ArraySequence _must_ have an Alphabet, not a MolType,
    because it is often important to store just a subset of the possible
    characters (e.g. the non-degenerate bases) for modeling purposes.
    """

    alphabet = None  # REPLACE IN SUBCLASSES
    moltype = None  # REPLACE IN SUBCLASSES
    delimiter = ""  # Used for string conversions
    line_wrap = 80  # Wrap sequences at 80 characters by default.

    def __init__(self, data="", alphabet=None, name=None, info=None, check="ignored"):
        """Initializes sequence from data and alphabet.

        WARNING: Does not validate the data or alphabet for compatibility.
        This is for speed. Use is_valid() to check whether the data
        is consistent with the alphabet.

        WARNING: If data has name and/or info, gets ref to same object rather
        than copying in each case.
        """
        if name is None and hasattr(data, "name"):
            name = data.name
        if info is None and hasattr(data, "info"):
            info = data.info
        # set the label
        self.name = name
        # override the class alphabet if supplied
        if alphabet is not None:
            self.alphabet = alphabet
        # if we haven't already set self._data (e.g. in a subclass __init__),
        # guess the data type and set it here
        if not hasattr(self, "_data"):
            # if data is a sequence, copy its data and alphabet
            if isinstance(data, ArraySequence):
                self._data = data._data
                self.alphabet = data.alphabet
            # if it's an array
            elif type(data) == ARRAY_TYPE:
                self._data = data
            else:  # may be set in subclass init
                data = bytes_to_string(data)
                self._from_sequence(data)

        self.moltype = self.alphabet.moltype
        self.info = info
        self._repr_policy = dict(num_pos=60)

    def __getitem__(self, *args):
        """__getitem__ returns char or slice, as same class."""
        if len(args) == 1 and not isinstance(args[0], slice):
            result = array([self._data[args[0]]])
        else:
            result = self._data.__getitem__(*args)
        return self.__class__(result)

    def __lt__(self, other):
        """compares based on string"""
        return str(self) < other

    def __eq__(self, other):
        """compares based on string"""
        return str(self) == other

    def _from_sequence(self, data):
        """Fills self using the values in data, via the alphabet."""
        if self.alphabet:
            indices = self.alphabet.to_indices(data)
            self._data = array(indices, self.alphabet.array_type)
        else:
            self._data = array(data)

    def __str__(self):
        """Uses alphabet to convert self to string, using delimiter."""
        if hasattr(self.alphabet, "to_string"):
            return self.alphabet.to_string(self._data)
        else:
            return self.delimiter.join(map(str, self.alphabet.from_indices(self._data)))

    def __len__(self):
        """Returns length of data."""
        return len(self._data)

    def to_fasta(self, make_seqlabel=None, block_size=60):
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

        return alignment_to_fasta({label: str(self)}, block_size=block_size)

    def to_phylip(self, name_len=28, label_len=30):
        """Return string of self in one line for PHYLIP, no newline.

        Default: max name length is 28, label length is 30.
        """
        return str(self.name)[:name_len].ljust(label_len) + str(self)

    def is_valid(self):
        """Checks that no items in self are out of the alphabet range."""
        return self._data == self._data.clip(0, len(self.alphabet) - 1)

    def __iter__(self):
        """iter returns characters of self, rather than slices."""
        if hasattr(self.alphabet, "to_string"):
            return iter(self.alphabet.to_string(self._data))

        return iter(self.alpabet.from_indices(self._data))

    def tostring(self):
        """to_string delegates to self._data."""
        return self._data.tostring()

    def gaps(self):
        """Returns array containing 1 where self has gaps, 0 elsewhere.

        WARNING: Only checks for standard gap character (for speed), and
        does not check for ambiguous gaps, etc.
        """
        return self._data == self.alphabet.gap_index

    def nongaps(self):
        """Returns array contining 0 where self has gaps, 1 elsewhere.

        WARNING: Only checks for standard gap character (for speed), and
        does not check for ambiguous gaps, etc.
        """
        return self._data != self.alphabet.gap_index

    def regap(self, other, strip_existing_gaps=False):
        """Inserts elements of self into gaps specified by other.

        WARNING: Only checks for standard gap character (for speed), and
        does not check for ambiguous gaps, etc.
        """
        if strip_existing_gaps:
            s = self.degap()
        else:
            s = self
        c = self.__class__
        a = self.alphabet.gapped
        result = zeros(len(other), a.array_type) + a.gap_index
        put(result, nonzero(other.nongaps()), s._data)
        return c(result)

    def degap(self):
        """Returns ungapped copy of self, not changing alphabet."""
        if not hasattr(self.alphabet, "gap") or self.alphabet.gap is None:
            return self.copy()
        d = take(self._data, nonzero(logical_not(self.gap_array()))[0])
        return self.__class__(d, alphabet=self.alphabet, name=self.name, info=self.info)

    def copy(self):
        """Returns copy of self, always separate object."""
        return self.__class__(
            self._data.copy(), alphabet=self.alphabet, name=self.name, info=self.info
        )

    def __contains__(self, item):
        """Returns true if item in self (converts to strings)."""
        return item in str(self)

    def disambiguate(self, *args, **kwargs):
        """Disambiguates self using strings/moltype. Should recode if demand."""
        return self.__class__(self.moltype.disambiguate(str(self), *args, **kwargs))

    def distance(self, other, function=None, use_indices=False):
        """Returns distance between self and other using function(i,j).

        other must be a sequence.

        function should be a function that takes two items and returns a
        number. To turn a 2D matrix into a function, use
        cogent3.util.miscs.DistanceFromMatrix(matrix).

        use_indices: if False, maps the indices onto items (e.g. assumes
        function relates the characters). If True, uses the indices directly.

        NOTE: Truncates at the length of the shorter sequence.

        Note that the function acts on two _elements_ of the sequences, not
        the two sequences themselves (i.e. the behavior will be the same for
        every position in the sequences, such as identity scoring or a function
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
            # use identity scoring
            shortest = min(len(self), len(other))
            if not hasattr(other, "_data"):
                other = self.__class__(other)
            distance = (self._data[:shortest] != other._data[:shortest]).sum()
        else:
            distance = 0
            if use_indices:
                self_seq = self._data
                if hasattr(other, "_data"):
                    other_seq = other._data
            else:
                self_seq = self.alphabet.from_indices(self._data)
                if hasattr(other, "_data"):
                    other_seq = other.alphabet.from_indices(other._data)
                else:
                    other_seq = other
            for first, second in zip(self_seq, other_seq):
                distance += function(first, second)
        return distance

    def matrix_distance(self, other, matrix, use_indices=False):
        """Returns distance between self and other using a score matrix.

        if use_indices is True (default is False), assumes that matrix is
        an array using the same indices that self uses.

        WARNING: the matrix must explicitly contain scores for the case where
        a position is the same in self and other (e.g. for a distance matrix,
        an identity between U and U might have a score of 0). The reason the
        scores for the 'diagonals' need to be passed explicitly is that for
        some kinds of distance matrices, e.g. log-odds matrices, the 'diagonal'
        scores differ from each other. If these elements are missing, this
        function will raise a KeyError at the first position that the two
        sequences are identical.
        """
        return self.distance(other, DistanceFromMatrix(matrix))

    def shuffle(self):
        """Returns shuffled copy of self"""
        return self.__class__(permutation(self._data), info=self.info)

    def gap_array(self):
        """Returns array of 0/1 indicating whether each position is a gap."""
        gap_indices = []
        a = self.alphabet
        for c in self.moltype.gaps:
            if c in a:
                gap_indices.append(a.index(c))
        gap_vector = None
        for i in gap_indices:
            if gap_vector is None:
                gap_vector = self._data == i
            else:
                gap_vector = logical_or(gap_vector, self._data == i)
        return gap_vector

    def gap_indices(self):
        """Returns list of gap indices."""
        return list(self.gap_array().nonzero()[0])

    def frac_same_gaps(self, other):
        """Returns fraction of positions where gaps match other's gaps."""
        if not other:
            return 0
        self_gaps = self.gap_array()
        if hasattr(other, "gap_array"):
            other_gaps = other.gap_array()
        elif hasattr(other, "gap_vector"):
            other_gaps = array(other.gap_vector())
        else:
            other_gaps = array(self.moltype.gap_vector(other))
        min_len = min(len(self), len(other))
        self_gaps, other_gaps = self_gaps[:min_len], other_gaps[:min_len]
        return (self_gaps == other_gaps).sum() / float(min_len)


class ArraySequence(ArraySequenceBase, SequenceI):
    """ArraySequence provides an array-based implementation of Sequence.

    Use ArraySequenceBase if you need a stripped-down, fast implementation.
    ArraySequence implements everything that SequenceI implements.

    See docstrings for ArraySequenceBase and SequenceI for information about
    these respective classes.
    """

    def strip_bad(self):
        """Returns copy of self with bad chars excised"""
        valid_indices = self._data < len(self.alphabet)
        result = compress(valid_indices, self._data)
        return self.__class__(result, info=self.info)

    def strip_bad_and_gaps(self):
        """Returns copy of self with bad chars and gaps excised."""
        gap_indices = list(map(self.alphabet.index, self.moltype.gaps))
        valid_indices = self._data < len(self.alphabet)
        for i in gap_indices:
            valid_indices[self._data == i] = False
        result = compress(valid_indices, self._data)
        return self.__class__(result, info=self.info)

    def strip_degenerate(self):
        """Returns copy of self without degenerate symbols.

        NOTE: goes via string intermediate because some of the algorithms
        for resolving degenerates are complex. This could be optimized if
        speed becomes critical.
        """
        return self.__class__(self.moltype.strip_degenerate(str(self)), info=self.info)

    def count_gaps(self):
        """Returns count of gaps in self."""
        return self.gap_array().sum()

    def gap_vector(self):
        """Returns list of bool containing whether each pos is a gap."""
        return list(map(bool, self.gap_array()))

    def gap_maps(self):
        """Returns dicts mapping gapped/ungapped positions."""
        nongaps = logical_not(self.gap_array())
        indices = arange(len(self)).compress(nongaps)
        new_indices = arange(len(indices))
        return (
            dict(list(zip(new_indices, indices))),
            dict(list(zip(indices, new_indices))),
        )

    def first_gap(self):
        """Returns position of first gap, or None."""
        a = self.gap_indices()
        try:
            return a[0]
        except IndexError:
            return None

    def is_gapped(self):
        """Returns True of sequence contains gaps."""
        return len(self.gap_indices())

    def mw(self, *args, **kwargs):
        """Returns molecular weight.

        Works via string intermediate: could optimize using array of MW if
        speed becomes important.
        """
        return self.moltype.mw(str(self), *args, **kwargs)

    def frac_similar(self, other, similar_pairs):
        """Returns fraction of positions where self[i] is similar to other[i].

        similar_pairs must be a dict such that d[(i,j)] exists if i and j are
        to be counted as similar. Use PairsFromGroups in cogent3.util.misc to
        construct such a dict from a list of lists of similar residues.

        Truncates at the length of the shorter sequence.

        Note: current implementation re-creates the distance function each
        time, so may be expensive compared to creating the distance function
        using for_seq separately.

        Returns 0 if one sequence is empty.

        NOTE: goes via string intermediate, could optimize using array if
        speed becomes important. Note that form of similar_pairs input would
        also have to change.
        """
        if not self or not other:
            return 0.0

        return for_seq(f=lambda x, y: (x, y) in similar_pairs, normalizer=per_shortest)(
            str(self), str(other)
        )

    def replace(self, oldchar, newchar):
        """return new instance with oldchar replaced by newchar"""
        oldindex = self.alphabet.index(oldchar)
        newindex = self.alphabet.index(newchar)
        new = self._data.copy()
        new[new == oldindex] = newindex
        return self.__class__(new, name=self.name, info=self.info)


class ArrayNucleicAcidSequence(ArraySequence):
    """Abstract class defining ops for codons, translation, etc."""

    def to_codons(self):
        """Returns copy of self in codon alphabet. Assumes ungapped."""
        alpha_len = len(self.alphabet)
        return ArrayCodonSequence(
            alpha_len * (alpha_len * self._data[::3] + self._data[1::3])
            + self._data[2::3],
            name=self.name,
            alphabet=self.alphabet**3,
        )

    def complement(self):
        """Returns complement of sequence"""
        return self.__class__(
            self.alphabet._complement_array.take(self._data), info=self.info
        )

    def rc(self):
        """Returns reverse-complement of sequence"""
        comp = self.alphabet._complement_array.take(self._data)
        return self.__class__(comp[::-1], info=self.info)

    def to_rna(self):
        """Returns self as RNA"""
        return ArrayRnaSequence(self._data)

    def to_dna(self):
        """Returns self as DNA"""
        return ArrayDnaSequence(self._data)


class ArrayRnaSequence(ArrayNucleicAcidSequence):
    moltype = None  # set to RNA in moltype.py
    alphabet = None  # set to RNA.alphabets.degen_gapped in moltype.py

    def __init__(self, data="", *args, **kwargs):
        """Returns new ArrayRnaSequence, converting T -> U"""
        if hasattr(data, "upper"):
            data = data.upper().replace("T", "U")
        return super(ArrayRnaSequence, self).__init__(data, *args, **kwargs)


class ArrayDnaSequence(ArrayNucleicAcidSequence):
    moltype = None  # set to DNA in moltype.py
    alphabet = None  # set to DNA.alphabets.degen_gapped in moltype.py

    def __init__(self, data="", *args, **kwargs):
        """Returns new ArrayDnaSequence, converting U -> T"""
        if hasattr(data, "upper"):
            data = data.upper().replace("U", "T")
        return super(ArrayDnaSequence, self).__init__(data, *args, **kwargs)


class ArrayCodonSequence(ArraySequence):
    """Abstract base class for codon sequences, incl. string conversion."""

    SequenceClass = ArrayNucleicAcidSequence

    def __str__(self):
        """Joins triplets together as string."""
        return self.delimiter.join(map("".join, self.alphabet.from_indices(self._data)))

    def _from_string(self, s):
        """Reads from a raw string, rather than a DnaSequence."""
        s = s.upper().replace("U", "T")  # convert to uppercase DNA
        d = self.SequenceClass(s, alphabet=self.alphabet.sub_enumerations[0])
        self._data = d.to_codons()._data

    def __init__(self, data="", alphabet=None, name=None, info=None):
        """Override __init__ to handle init from string."""
        if isinstance(data, str):
            self._from_string(data)
        ArraySequence.__init__(self, data, alphabet, name, info=info)

    def to_codons(self):
        """Converts self to codons -- in practice, just returns self.

        Supports interface of other NucleicAcidSequences."""
        return self

    def to_dna(self):
        """Returns a ArrayDnaSequence from the data in self"""
        unpacked = self.alphabet.unpack_arrays(self._data)
        result = zeros((len(self._data), 3))
        for i, v in enumerate(unpacked):
            result[:, i] = v
        return ArrayDnaSequence(ravel(result), name=self.name)

    def to_rna(self):
        """Returns a ArrayDnaSequence from the data in self."""
        unpacked = self.alphabet.unpack_arrays(self._data)
        result = zeros((len(self._data), 3))
        for i, v in enumerate(unpacked):
            result[:, i] = v
        return ArrayRnaSequence(ravel(result), name=self.name)


class ArrayDnaCodonSequence(ArrayCodonSequence):
    """Holds non-degenerate DNA codon sequence."""

    alphabet = None  # set to DNA.alphabets.base.Triples in moltype.py
    SequenceClass = ArrayDnaSequence


class ArrayRnaCodonSequence(ArrayCodonSequence):
    """Holds non-degenerate DNA codon sequence."""

    alphabet = None  # set to RNA.alphabets.base.Triples in motype.py
    SequenceClass = ArrayRnaSequence

    def _from_string(self, s):
        """Reads from a raw string, rather than a DnaSequence."""
        s = s.upper().replace("T", "U")  # convert to uppercase DNA
        d = self.SequenceClass(s, alphabet=self.alphabet.sub_enumerations[0])
        self._data = d.to_codons()._data


class ArrayProteinSequence(ArraySequence):
    moltype = None  # set to PROTEIN in moltype.py
    alphabet = None  # set to PROTEIN.alphabets.degen_gapped in moltype.py


class ArrayProteinWithStopSequence(ArraySequence):
    moltype = None  # set to PROTEIN_WITH_STOP in moltype.py
    alphabet = None  # set to PROTEIN_WITH_STOP.alphabets.degen_gapped in moltype.py


@singledispatch
def _coerce_seq(data, preserve_case, checker):
    from cogent3.core.alignment import Aligned

    if isinstance(data, Aligned):
        return _coerce_seq(str(data), preserve_case, checker)
    raise NotImplementedError(f"{type(data)}")


@_coerce_seq.register
def _(data: SeqView, preserve_case, checker):
    return data


@_coerce_seq.register
def _(data: Sequence, preserve_case, checker):
    return _coerce_seq(str(data), preserve_case, checker)


@_coerce_seq.register
def _(data: ArraySequence, preserve_case, checker):
    return _coerce_seq(str(data), preserve_case, checker)


@_coerce_seq.register
def _(data: str, preserve_case, checker):
    if not preserve_case:
        data = data.upper()
    checker(data)
    return SeqView(data)


@_coerce_seq.register
def _(data: bytes, preserve_case, checker):
    if not preserve_case:
        data = data.upper()
    data = data.decode("utf8")
    checker(data)
    return SeqView(data)


@_coerce_seq.register
def _(data: tuple, preserve_case, checker):
    return _coerce_seq("".join(data), preserve_case, checker)


@_coerce_seq.register
def _(data: list, preserve_case, checker):
    return _coerce_seq("".join(data), preserve_case, checker)

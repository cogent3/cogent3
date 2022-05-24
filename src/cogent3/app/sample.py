from collections import defaultdict

from numpy import array
from numpy import random as np_random

from cogent3.core.alignment import Alignment, ArrayAlignment
from cogent3.core.genetic_code import get_code
from cogent3.core.moltype import get_moltype

from .composable import (
    ALIGNED_TYPE,
    SEQUENCE_TYPE,
    SERIALISABLE_TYPE,
    ComposableAligned,
    ComposableSeq,
    NotCompleted,
)
from .translate import get_fourfold_degenerate_sets


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


# TODO need a function to filter sequences based on divergence, ala divergent
# set.


def intersection(groups):
    """returns the intersection of all groups"""
    common = set(groups.pop())
    return common.intersection(*map(set, groups))


def union(groups):
    """returns the intersection of all groups"""
    union = set(groups.pop())
    union = union.union(*map(set, groups))
    return union


class concat:
    """Creates a concatenated alignment from a series. Returns an Alignment."""

    def __init__(self, join_seq="", intersect=True, moltype=None):
        """concatenate sequences from a series of alignments

        Parameters
        ----------
        join_seq : str
            splice sequences together using this string
        intersect : bool
            result contains only sequences present in all alignments. If False,
            missings sequences will be replaced by a sequence of '?'.
        moltype : str
            molecular type, must be either DNA or RNA
        """
        self._name_callback = {True: intersection}.get(intersect, union)
        self._intersect = intersect
        self._moltype = moltype
        self._join_seq = join_seq

    def concat(self, data):
        """returns an alignment

        Parameters
        ----------
        data
            series of alignment instances
        """
        if len(data) == 0:
            raise ValueError("no data")

        names = []
        for aln in data:
            if not (isinstance(aln, ArrayAlignment) or isinstance(aln, Alignment)):
                raise TypeError(f"{type(aln)} invalid for concat")
            names.append(aln.names)

        names = self._name_callback(names)
        collated = defaultdict(list)
        if self._moltype is None:
            self._moltype = aln.moltype

        for aln in data:
            if self._moltype and aln.moltype != self._moltype:
                # try converting
                aln = aln.to_moltype(self.moltype)

            if self._intersect:
                seqs = aln.take_seqs(names).to_dict()
            else:
                seqs = defaultdict(lambda: "?" * len(aln))
                seqs.update(aln.to_dict())

            for name in names:
                collated[name].append(seqs[name])

        combined = {n: self._join_seq.join(collated[n]) for n in names}
        aln = ArrayAlignment(data=combined, moltype=self._moltype)
        return aln

    __call__ = concat


class omit_degenerates(ComposableAligned):
    """Excludes alignment columns with degenerate conditions. Can accomodate
    reading frame. Returns modified Alignment."""

    _input_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(self, moltype=None, gap_is_degen=True, motif_length=1):
        """excludes degenerate characters from alignment

        Parameters
        ----------
        moltype : str
            molecular type, must be either DNA or RNA
        gap_is_degen : bool
            include gap character in degenerate character set
        motif_length : int
            sequences split into non-overlapping tuples of this size. If a
            tuple contains a degen character at any position the entire tuple
            is excluded
        """
        super(omit_degenerates, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
            assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self.moltype = moltype
        self._no_degen = omit_degenerates
        self._allow_gap = not gap_is_degen
        self._motif_length = motif_length
        self.func = self.filter_degenerates

    def filter_degenerates(self, aln):
        if self.moltype and aln.moltype != self.moltype:
            # try converting
            aln = aln.to_moltype(self.moltype)

        result = aln.no_degenerates(
            motif_length=self._motif_length, allow_gap=self._allow_gap
        )
        if not result:
            result = NotCompleted(
                "FAIL", self, "all columns contained degenerates", source=aln
            )

        return result


class omit_gap_pos(ComposableAligned):
    """Excludes gapped alignment columns meeting a threshold. Can accomodate
    reading frame. Returns modified Alignment."""

    _input_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(self, allowed_frac=0.99, motif_length=1, moltype=None):
        """
        Parameters
        ----------
        allowed_frac : float
            columns with a fraction of gap characters exceeding allowed_frac are
            excluded
        motif_length : int
            sequences split into non-overlapping tuples of this size.
        moltype : str
            molecular type, must be either DNA or RNA
        """
        super(omit_gap_pos, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
            assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self.moltype = moltype
        self._allowed_frac = allowed_frac
        self._motif_length = motif_length
        self.func = self.omit

    def omit(self, aln):
        if self.moltype and aln.moltype != self.moltype:
            # try converting
            aln = aln.to_moltype(self.moltype)

        result = aln.omit_gap_pos(
            allowed_gap_frac=self._allowed_frac, motif_length=self._motif_length
        )
        if not result:
            result = NotCompleted(
                "FAIL", self, "all columns exceeded gap threshold", source=aln
            )

        return result


class take_codon_positions(ComposableAligned):
    """Extracts the specified codon position(s) from an alignment.
    Returns an Alignment."""

    _input_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(
        self,
        *positions,
        fourfold_degenerate=False,
        gc="Standard Nuclear",
        moltype="dna",
    ):
        """selects the indicated codon positions from an alignment

        Parameters
        ----------
        positions
            either an integer (1, 2, 3), or a tuple of position numbers,
            e.g. 3 is third position, (1,2) is first and second codon position
        fourfold_degenerate : bool
            if True, returns third positions from four-fold degenerate codons.
            Overrides positions.
        gc
            identifer for a genetic code or a genetic code instance
        moltype : str
            molecular type, must be either DNA or RNA
        """
        super(take_codon_positions, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        assert moltype is not None
        moltype = get_moltype(moltype)

        assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self._moltype = moltype
        self._four_fold_degen = fourfold_degenerate
        self._fourfold_degen_sets = None

        if fourfold_degenerate:
            gc = get_code(gc)
            sets = get_fourfold_degenerate_sets(
                gc, alphabet=moltype.alphabet, as_indices=True
            )
            self._fourfold_degen_sets = sets
            self.func = self.take_fourfold_positions
            return

        assert (
            1 <= min(positions) <= 3 and 1 <= max(positions) <= 3
        ), "Invalid codon positions"

        by_index = True if len(positions) == 1 else False
        if by_index:
            positions = positions[0] - 1
            self.func = self.take_codon_position
        else:
            positions = tuple(p - 1 for p in sorted(positions))
            self.func = self.take_codon_positions

        self._positions = positions

    def take_fourfold_positions(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        fourfold_codon_sets = self._fourfold_degen_sets

        def ffold(x):
            x = set(tuple(e) for e in list(x))
            for codon_set in fourfold_codon_sets:
                if x <= codon_set:
                    return True
            return False

        new = aln.filtered(ffold, motif_length=3)
        return new[2::3]

    def take_codon_position(self, aln):
        if isinstance(aln, Alignment):
            indices = list(range(self._positions, len(aln), 3))
            result = aln.take_positions(indices)
        elif isinstance(aln, ArrayAlignment):
            result = aln[self._positions :: 3]
        return result

    def take_codon_positions(self, aln):
        """takes multiple positions"""
        length = len(aln)
        indices = [k for k in range(length) if k % 3 in self._positions]
        return aln.take_positions(indices)


class take_named_seqs(ComposableSeq):
    """Extracts (or everything but) named sequences. Returns a filtered
    sequences, alignment that satisified the condition, NotCompleted otherwise."""

    _input_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment", "SequenceCollection")

    def __init__(self, *names, negate=False):
        """selects named sequences from a collection

        Returns
        -------
        A new sequence collection, or False if not all the named sequences are
        in the collection.
        """
        super(take_named_seqs, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        self._names = names
        self._negate = negate
        self.func = self.take_seqs

    def take_seqs(self, data):
        try:
            data = data.take_seqs(self._names, negate=self._negate)
        except KeyError:
            missing = set(self._names) - set(data.names)
            msg = f"named seq(s) {missing} not in {data.names}"
            data = NotCompleted("FALSE", self, msg, source=data)
        return data


class take_n_seqs(ComposableSeq):
    """Selects n sequences from a collection. Chooses first n sequences, or selects randomly if specified.

    Returns original object type with the selected sequences, NotCompleted object otherwise.
    """

    _input_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment", "SequenceCollection")

    def __init__(self, number, random=False, seed=None, fixed_choice=True):
        """
        selects n sequences from a collection

        Parameters
        ----------
        number: int
            number of sequences to sample. If number of sequences in a collectionis < n, returns NotCompleted
            indicating a FAIL.
        random: bool
            Whether to choose the sequences randomly.
        seed: int
            Seed for the numpy random number generator.
        fixed_choice: bool
            sequence names selected from the first alignment are used for all others.

        Returns
        -------
        A new sequence collection, or NotCompleted if not insufficient sequences are in the collection.
        """
        super(take_n_seqs, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()

        if seed:
            np_random.seed(seed)

        self._names = None
        self._number = number
        self._random = random
        self._fixed_choice = fixed_choice
        self.func = self.take_seqs

    def _set_names(self, data):
        """set the names attribute"""
        if not self._random:
            self._names = data.names[: self._number]
            return

        self._names = np_random.choice(data.names, self._number, replace=False)

    def take_seqs(self, data):
        if len(data.names) < self._number:
            return NotCompleted("FALSE", self.take_seqs, "not enough sequences")

        if self._names is None or not self._fixed_choice:
            self._set_names(data)

        try:
            data = data.take_seqs(self._names)
        except KeyError:
            missing = set(self._names) - set(data.names)
            msg = f"named seq(s) {missing} not in {data.names}"
            data = NotCompleted("FALSE", self, msg, source=data)
        return data


class min_length(ComposableSeq):
    """Filters sequence collections / alignments by length. Returns the
    data if it satisfies the condition, NotCompleted otherwise."""

    _input_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment", "SequenceCollection")

    def __init__(self, length, motif_length=1, subtract_degen=True, moltype=None):
        """
        Parameters
        ----------
        length : int
            only alignments with this length returned, False otherwise
        motif_length : int
            length is converted to modulo motif_length
        subtract_degen : bool
            degenerate characters subtracted from sequence length calculation
        moltype
            molecular type, can be string or instance
        """
        super(min_length, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        if motif_length > 1:
            length = length // motif_length
        self._min_length = length
        self._motif_length = motif_length
        self.func = self.if_long_enough
        self._subtract_degen = subtract_degen
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype

    def if_long_enough(self, data):
        if self._moltype and self._moltype != data.moltype:
            data = data.to_moltype(self._moltype)

        if self._subtract_degen:
            if not hasattr(data.alphabet, "non_degen"):
                name = self.__class__.__name__
                msg = (
                    "%s(subtract_degen=True) requires DNA, RNA or PROTEIN " "moltype"
                ) % name
                raise ValueError(msg)

        lengths = data.get_lengths(
            allow_gap=not self._subtract_degen,
            include_ambiguity=not self._subtract_degen,
        )
        length, _ = min([(l, n) for n, l in lengths.items()])

        if length < self._min_length:
            msg = f"{length} < min_length {self._min_length}"
            data = NotCompleted("FALSE", self, msg, source=data)

        return data


class _GetStart:
    choose = np_random.choice

    def __init__(self, start):
        self._start = start
        self.func = {True: self._int}.get(type(start) == int, self._rand)

    def _int(self, length):
        return self._start

    def _rand(self, length):
        return self.choose(length)

    def __call__(self, length):
        return self.func(length)


class fixed_length(ComposableAligned):
    """Sample an alignment to a fixed length. Returns an Alignment of the
    specified length, or NotCompleted if alignment too short."""

    _input_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(
        self, length, start=0, random=False, seed=None, motif_length=1, moltype=None
    ):
        """
        Parameters
        ----------
        length : int
            only alignments with this length returned, False otherwise
        start
            integer starting position for truncation, or 'random' in which case
            a random start is chosen (within the possible range returning an
            alignment of the specified length). Overrides  `random`.
        random : bool
            random positions for the corresponding tuple are chosen.
        seed : int
            random number seed
        motif_length : int
            length of sequence units to consider. If not 1, length and start are
            converted (reduced) if necessary to be modulo motif_length
        moltype
            molecular type, can be string or instance
        """
        super(fixed_length, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        diff = length % motif_length
        if diff != 0:
            length -= diff
        assert length % motif_length == 0

        self._length = length
        self._motif_length = motif_length
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype
        if type(start) == str:
            assert start.lower().startswith("rand")
            random = False
        else:
            assert type(start) == int
            assert start >= 0
            diff = start % motif_length
            if diff != 0:
                start -= diff

        self._start = _GetStart(start)
        if seed:
            np_random.seed(seed)

        self.func = {False: self.truncated}.get(random, self.sample_positions)

    def truncated(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        if len(aln) < self._length:
            msg = f"{len(aln)} < min_length {self._length}"
            result = NotCompleted("FALSE", self.__class__.__name__, msg, source=aln)
        else:
            start = self._start(len(aln) - self._length)
            result = aln[start : start + self._length]

        return result

    def sample_positions(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        indices = array(range(len(aln)))
        if self._motif_length == 1:
            number = self._length
        else:
            number = self._length // self._motif_length

        pos = np_random.choice(
            indices.shape[0] // self._motif_length, number, replace=False
        )
        if self._motif_length == 1:
            result = indices[pos]
        else:
            dim = indices.shape[0] // self._motif_length
            indices.resize((dim, self._motif_length))
            result = indices[pos, :]

        result.sort(axis=0)
        result = aln.take_positions(result.flatten().tolist())
        return result


class omit_bad_seqs(ComposableAligned):
    """Eliminates sequences from Alignment based on gap fraction, unique gaps.
    Returns modified alignment."""

    _input_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(self, quantile=None, gap_fraction=1, moltype="dna"):
        """Returns an alignment without the sequences responsible for
        exceeding disallowed_frac.

        Parameters
        ----------
        quantile : float or None
            The number of gaps uniquely introduced by a sequence are counted.
            The value corresponding to quantile is determined and all sequences
            whose unique gap count is larger than this cutoff are excluded.
            If None, this condition is not applied.
        gap_fraction
            sequences whose proportion of gaps is >= this value are excluded, the
            default excludes sequences that are just gaps.
        moltype
            molecular type, can be string or instance
        """
        super(omit_bad_seqs, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        assert (
            moltype.label.lower() in "dna rna protein protein_with_stop"
        ), "moltype must be one of DNA, RNA or PROTEIN"
        self._quantile = quantile
        self._gap_fraction = gap_fraction
        self._moltype = moltype
        self.func = self.drop_bad_seqs

    def drop_bad_seqs(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        gaps_per_seq = aln.count_gaps_per_seq()
        length = len(aln)
        keep = [n for n, c in gaps_per_seq.items() if c / length < self._gap_fraction]
        result = aln.take_seqs(keep)
        if self._quantile is not None:
            result = result.omit_bad_seqs(quantile=self._quantile)
        return result


class omit_duplicated(ComposableSeq):
    """Removes redundant sequences, recording dropped sequences in
    seqs.info.dropped. Returns sequence collection with only unique sequences."""

    _input_types = (SEQUENCE_TYPE, SERIALISABLE_TYPE, ALIGNED_TYPE)
    _output_types = (SEQUENCE_TYPE, SERIALISABLE_TYPE, ALIGNED_TYPE)
    _data_types = ("ArrayAlignment", "Alignment", "SequenceCollection")

    def __init__(self, mask_degen=False, choose="longest", seed=None, moltype=None):
        """Returns unique sequences, adds 'dropped' key to seqs.info

        Parameters
        ----------
        mask_degen
            if True, degenerate characters are ignored
        choose
            choose a representative from sets of duplicated sequences.
            Valid values are None (all members of a duplicated set are excluded),
            'longest', 'random'.
        seed : int
            set random number seed. Only applied of choose=='random'
        moltype
            molecular type, can be string or instance
        """
        super(omit_duplicated, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )

        assert not choose or choose in "longestrandom"
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype
        if choose == "random" and seed:
            np_random.seed(seed)

        self._mask_degen = mask_degen
        if choose == "longest":
            self.func = self.choose_longest
        elif choose == "random":
            self.func = self.choose_random
        else:
            self.func = self.take_unique

    def choose_longest(self, seqs):
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        duplicates = seqs.get_identical_sets(mask_degen=self._mask_degen)
        lengths = seqs.get_lengths()
        excludes = []
        for group in duplicates:
            group_lengths = [(lengths[n], n) for n in group]
            group_lengths.sort(reverse=True)
            excludes.extend([n for l, n in group_lengths[1:]])

        seqs = seqs.take_seqs(excludes, negate=True)
        return seqs

    def choose_random(self, seqs):
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        duplicates = seqs.get_identical_sets(mask_degen=self._mask_degen)
        excludes = []
        for group in duplicates:
            chosen = np_random.choice([e for e in group])
            group.remove(chosen)
            excludes.extend(group)

        seqs = seqs.take_seqs(excludes, negate=True)
        return seqs

    def take_unique(self, seqs):
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        duplicates = seqs.get_identical_sets(mask_degen=self._mask_degen)
        names = set()
        for dupes in duplicates:
            names.update(dupes)
        seqs = seqs.take_seqs(names, negate=True)
        return seqs


class trim_stop_codons(ComposableSeq):
    """Removes terminal stop codons. Returns sequences / alignment."""

    _input_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _output_types = (SEQUENCE_TYPE, ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment", "SequenceCollection")

    def __init__(self, gc=1):
        """selects named sequences from a collection

        Parameters
        ----------
        gc
            identifier for a genetic code or a genetic code instance, defaults
            to standard genetic code

        Returns
        -------
        A new sequence collection, or False if not all the named sequences are
        in the collection.
        """
        super(trim_stop_codons, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        self._gc = gc
        self.func = self.trim_stops

    def trim_stops(self, data):
        data = data.trim_stop_codons(gc=self._gc)
        return data

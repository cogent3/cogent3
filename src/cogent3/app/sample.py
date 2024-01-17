from collections import defaultdict
from typing import List, Optional, Union

from numpy import array
from numpy import random as np_random

from cogent3.core.alignment import Alignment, ArrayAlignment
from cogent3.core.genetic_code import get_code
from cogent3.core.moltype import get_moltype

from .composable import NON_COMPOSABLE, NotCompleted, define_app
from .translate import get_fourfold_degenerate_sets
from .typing import AlignedSeqsType, SeqsCollectionType, SerialisableType


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


@define_app(app_type=NON_COMPOSABLE)
class concat:
    """Creates a concatenated alignment from a series."""

    def __init__(
        self, join_seq: str = "", intersect: bool = True, moltype: Optional[str] = None
    ):
        """
        Parameters
        ----------
        join_seq : str
            splice sequences together using this string
        intersect : bool
            result contains only sequences present in all alignments. If False,
            missings sequences will be replaced by a sequence of '?'.
        moltype : str
            molecular type, must be either DNA or RNA

        Examples
        --------

        Create an app to concatenate two alignments

        >>> from cogent3 import app_help, get_app, make_aligned_seqs
        >>> concat_alns = get_app("concat", moltype="dna")

        Create sample alignments with matching sequence names

        >>> aln1 = make_aligned_seqs({"s1": "AAA", "s2": "CAA", "s3": "AAA"}, moltype="dna")
        >>> aln2 = make_aligned_seqs({"s1": "GCG", "s2": "GGG", "s3": "GGT"}, moltype="dna")

        Concatenate alignments. By default, sequences without matching names in
        the corresponding alignment are omitted (intersect=True).

        >>> result = concat_alns([aln1, aln2])
        >>> print(result.to_pretty(name_order=["s1","s2","s3"]))
        s1    AAAGCG
        s2    C...G.
        s3    ....GT

        Create an app that includes missing sequences across alignments.
        Missing sequences are replaced by a sequence of "?".

        >>> concat_missing = get_app("concat", moltype="dna", intersect=False)
        >>> aln3 = make_aligned_seqs({"s4": "GCG", "s5": "GGG"}, moltype="dna")
        >>> result = concat_missing([aln1, aln3])
        >>> print(result.to_pretty(name_order=["s1","s2","s3","s4","s5"]))
        s1    AAA???
        s2    C.....
        s3    ......
        s4    ???GCG
        s5    ???GGG

        Create an app that delimits concatenated alignments with "N"

        >>> concat_delim = get_app("concat", join_seq="N", moltype="dna")
        >>> result = concat_delim([aln1, aln2])
        >>> print(result.to_pretty(name_order=["s1","s2","s3"]))
        s1    AAANGCG
        s2    C....G.
        s3    .....GT
        """
        self._name_callback = {True: intersection}.get(intersect, union)
        self._intersect = intersect
        self._moltype = moltype
        self._join_seq = join_seq

    def main(
        self, data: List[AlignedSeqsType]
    ) -> Union[SerialisableType, AlignedSeqsType]:
        """returns an alignment

        Parameters
        ----------
        data
            series of alignment instances
        """
        if not data:
            raise ValueError("no data")

        names = []
        for aln in data:
            if self._moltype is None:
                self._moltype = aln.moltype

            if not isinstance(aln, (ArrayAlignment, Alignment)):
                raise TypeError(f"{type(aln)} invalid for concat")
            names.append(aln.names)

        names = self._name_callback(names)
        collated = defaultdict(list)

        for aln in data:
            if self._moltype and aln.moltype != self._moltype:
                # try converting
                aln = aln.to_moltype(self._moltype)

            if self._intersect:
                seqs = aln.take_seqs(names).to_dict()
            else:
                seqs = defaultdict(lambda: "?" * len(aln))
                seqs |= aln.to_dict()

            for name in names:
                collated[name].append(seqs[name])

        combined = {n: self._join_seq.join(collated[n]) for n in names}
        aln = ArrayAlignment(data=combined, moltype=self._moltype)
        return aln


@define_app
class omit_degenerates:
    """Excludes alignment columns with degenerate characters. Can accomodate
    reading frame."""

    def __init__(
        self, moltype: str = None, gap_is_degen: bool = True, motif_length: int = 1
    ):
        """
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

        Examples
        --------
        Degenerate IUPAC base symbols represents a site position that can have
        multiple possible nucleotides. For example, "Y" represents
        pyrimidines where the site can be either "C" or "T".

        Note: In molecular evolutionary and phylogenetic analyses, the gap
        character "-" is considered to be any base "N".

        Create sample data with degenerate characters

        >>> from cogent3 import app_help, get_app, make_aligned_seqs
        >>> aln = make_aligned_seqs({"s1": "ACGA-GACG", "s2": "GATGATGYT"}, moltype="dna")

        Create an app that omits aligned columns containing a degenerate
        character from an alignment

        >>> app = get_app("omit_degenerates", moltype="dna")
        >>> result = app(aln)
        >>> print(result.to_pretty())
        s1    ACGAGAG
        s2    GATGTGT

        Create an app which omits degenerate characters, but retains gaps

        >>> app = get_app("omit_degenerates", moltype="dna", gap_is_degen=False)
        >>> result = app(aln)
        >>> print(result.to_pretty())
        s1    ACGA-GAG
        s2    GATGATGT

        Split sequences into non-overlapping tuples of length 2 and exclude
        any tuple that contains a degenerate character

        >>> app = get_app("omit_degenerates", moltype="dna", motif_length=2)
        >>> result = app(aln)
        >>> print(result.to_pretty())
        s1    ACGA
        s2    GATG

        A NotCompleted object (see https://cogent3.org/doc/app/not-completed.html)
        is returned if the moltype is not specified in the alignment or app

        >>> aln = make_aligned_seqs({"s1": "ACGA-GACG", "s2": "GATGATGYT"})
        >>> app = get_app("omit_degenerates")
        >>> result = app(aln)
        >>> result.message
        'Traceback...
        """
        if moltype:
            moltype = get_moltype(moltype)
            assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self._moltype = moltype
        self._allow_gap = not gap_is_degen
        self._motif_length = motif_length

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, aln: AlignedSeqsType) -> T:
        if self._moltype and aln.moltype != self._moltype:
            # try converting
            aln = aln.to_moltype(self._moltype)

        return aln.no_degenerates(
            motif_length=self._motif_length, allow_gap=self._allow_gap
        ) or NotCompleted("FAIL", self, "all columns contained degenerates", source=aln)


@define_app
class omit_gap_pos:
    """Excludes gapped alignment columns meeting a threshold. Can accomodate
    reading frame."""

    def __init__(
        self, allowed_frac: float = 0.99, motif_length: int = 1, moltype: str = None
    ):
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

        Examples
        --------

        Create a sample alignment and an app that excludes highly gapped sites.
        Sites with over 99% gaps are excluded by default

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> aln = make_aligned_seqs({"s1": "ACGA-GA-CG", "s2": "GATGATG-AT"})

        >>> app = get_app("omit_gap_pos", moltype="dna")
        >>> result = app(aln)
        >>> print(result.to_pretty(name_order=["s1", "s2"]))
        s1    ACGA-GACG
        s2    GATGATGAT

        Create an app that excludes all aligned sites with over 49% gaps

        >>> app = get_app("omit_gap_pos", allowed_frac=0.49, moltype="dna")
        >>> result = app(aln)
        >>> print(result.to_pretty(name_order=["s1", "s2"]))
        s1    ACGAGACG
        s2    GATGTGAT

        Create an app that excludes non-overlapping alignment columns of three
        sites that have over 17% gaps. Trailing columns that are less than
        three sites long are excluded

        >>> app = get_app("omit_gap_pos", allowed_frac=0.17, motif_length=3, moltype="dna")
        >>> result = app(aln)
        >>> print(result.to_pretty(name_order=["s1", "s2"]))
        s1    ACGA-G
        s2    GATGAT

        An error is returned if all sites are excluded

        >>> aln = make_aligned_seqs({"s1": "ACGA------", "s2": "----ATG-AT"})
        >>> app = get_app("omit_gap_pos", allowed_frac=0.17, motif_length=3, moltype="dna")
        >>> result = app(aln)
        >>> result.message
        'all columns exceeded gap threshold'
        """
        if moltype:
            moltype = get_moltype(moltype)
            assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self._moltype = moltype
        self._allowed_frac = allowed_frac
        self._motif_length = motif_length

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, aln: AlignedSeqsType) -> T:
        if self._moltype and aln.moltype != self._moltype:
            # try converting
            aln = aln.to_moltype(self._moltype)

        return aln.omit_gap_pos(
            allowed_gap_frac=self._allowed_frac, motif_length=self._motif_length
        ) or NotCompleted(
            "FAIL", self, "all columns exceeded gap threshold", source=aln
        )


@define_app
class take_codon_positions:
    """Extracts the specified codon position(s) from an alignment."""

    def __init__(
        self,
        *positions,
        fourfold_degenerate=False,
        gc="Standard",
        moltype="dna",
    ):
        """
        Parameters
        ----------
        positions
            either an integer (1, 2, 3), or a tuple of position numbers,
            e.g. 3 is third position, (1,2) is first and second codon position
        fourfold_degenerate : bool
            if True, returns third positions from four-fold degenerate codons.
            Overrides positions.
        gc
            identifier for a genetic code or a genetic code instance
        moltype : str
            molecular type, must be either DNA or RNA
        """
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
            self._func = self.take_fourfold_positions
            return

        assert (
            1 <= min(positions) <= 3 and 1 <= max(positions) <= 3
        ), "Invalid codon positions"

        by_index = len(positions) == 1
        if by_index:
            positions = positions[0] - 1
            self._func = self.take_codon_position
        else:
            positions = tuple(p - 1 for p in sorted(positions))
            self._func = self.take_codon_positions

        self._positions = positions

    def take_fourfold_positions(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        fourfold_codon_sets = self._fourfold_degen_sets

        def ffold(x):
            x = {tuple(e) for e in list(x)}
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

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, aln: AlignedSeqsType) -> T:
        return self._func(aln)


@define_app
class take_named_seqs:
    """Selects named sequences from a collection."""

    def __init__(self, *names, negate=False):
        """
        Parameters
        ----------
        *names
            series of sequence names
        negate
            if True, excludes the provided names from the result
        """
        self._names = names
        self._negate = negate

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, data: SeqsCollectionType) -> T:
        try:
            data = data.take_seqs(self._names, negate=self._negate)
        except KeyError:
            missing = set(self._names) - set(data.names)
            msg = f"named seq(s) {missing} not in {data.names}"
            data = NotCompleted("FALSE", self, msg, source=data)
        return data


@define_app
class take_n_seqs:
    """Selects n sequences from a collection. Chooses first n sequences, or
    selects randomly if specified."""

    def __init__(self, number, random=False, seed=None, fixed_choice=True):
        """
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
        if seed:
            np_random.seed(seed)

        self._names = None
        self._number = number
        self._random = random
        self._fixed_choice = fixed_choice

    def _set_names(self, data):
        """set the names attribute"""
        if not self._random:
            self._names = data.names[: self._number]
            return

        self._names = np_random.choice(data.names, self._number, replace=False).tolist()

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, data: SeqsCollectionType) -> T:
        """returns data with n sequences"""
        if len(data.names) < self._number:
            return NotCompleted("FALSE", self.main, "not enough sequences")

        if self._names is None or not self._fixed_choice:
            self._set_names(data)

        try:
            data = data.take_seqs(self._names)
        except KeyError:
            missing = set(self._names) - set(data.names)
            msg = f"named seq(s) {missing} not in {data.names}"
            data = NotCompleted("FALSE", self, msg, source=data)
        return data


@define_app
class min_length:
    """Filters sequence collections / alignments by length."""

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

        Examples
        --------
        Create a sample alignment. Alignments must have a moltype specified

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> aln = make_aligned_seqs({"s1": "ACGGT", "s2": "AC-GT"}, moltype="dna")

        Create the app to filter sequences with a minimum length of 3 sites

        >>> app = get_app("min_length", 3)
        >>> result = app(aln)
        >>> len(result)
        5

        When all sequences are shorter than the min length, returns a
        NotCompleted (see https://cogent3.org/doc/app/not-completed.html)

        >>> app = get_app("min_length", 7)
        >>> result = app(aln)
        >>> result.message
        '4 < min_length 7'

        If the moltype of the alignment is not specified, returns a NotCompleted

        >>> aln = make_aligned_seqs({"s1": "ACGGT", "s2": "AC-GT"})
        >>> app = get_app("min_length", 3)
        >>> result = app(aln)
        >>> result.message
        'Traceback...
        """
        if motif_length > 1:
            length = length // motif_length
        self._min_length = length
        self._motif_length = motif_length
        self._subtract_degen = subtract_degen
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, data: SeqsCollectionType) -> T:
        if self._moltype and self._moltype != data.moltype:
            data = data.to_moltype(self._moltype)

        if self._subtract_degen and not hasattr(data.alphabet, "non_degen"):
            raise ValueError(
                f"{self.__class__.__name__}(subtract_degen=True) requires DNA, RNA or PROTEIN "
                "moltype"
            )

        lengths = data.get_lengths(
            allow_gap=not self._subtract_degen,
            include_ambiguity=not self._subtract_degen,
        )
        length, _ = min((l, n) for n, l in lengths.items())

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


@define_app
class fixed_length:
    """Sample an alignment to a fixed length."""

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

        self._func = {False: self.truncated}.get(random, self.sample_positions)

    def truncated(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        if len(aln) < self._length:
            msg = f"{len(aln)} < min_length {self._length}"
            return NotCompleted("FALSE", self.__class__.__name__, msg, source=aln)
        else:
            start = self._start(len(aln) - self._length)
            return aln[start : start + self._length]

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

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, data: AlignedSeqsType) -> T:
        """return a fixed length alignment"""
        return self._func(data)


@define_app
class omit_bad_seqs:
    """Eliminates sequences from Alignment based on gap fraction, unique gaps."""

    def __init__(self, quantile=None, gap_fraction=1, moltype="dna"):
        """
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
        if moltype:
            moltype = get_moltype(moltype)
        assert (
            moltype.label.lower() in "dna rna protein protein_with_stop"
        ), "moltype must be one of DNA, RNA or PROTEIN"
        self._quantile = quantile
        self._gap_fraction = gap_fraction
        self._moltype = moltype

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, aln: AlignedSeqsType) -> T:
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        gaps_per_seq = aln.count_gaps_per_seq()
        length = len(aln)
        keep = [n for n, c in gaps_per_seq.items() if c / length < self._gap_fraction]
        result = aln.take_seqs(keep)
        if self._quantile is not None:
            result = result.omit_bad_seqs(quantile=self._quantile)
        return result


@define_app
class omit_duplicated:
    """Removes redundant sequences, recording dropped sequences in
    seqs.info.dropped."""

    def __init__(self, mask_degen=False, choose="longest", seed=None, moltype=None):
        """
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
        assert not choose or choose in "longestrandom"
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype
        if choose == "random" and seed:
            np_random.seed(seed)

        self._mask_degen = mask_degen
        if choose == "longest":
            self._func = self.choose_longest
        elif choose == "random":
            self._func = self.choose_random
        else:
            self._func = self.take_unique

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
            chosen = np_random.choice(list(group))
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

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, seqs: SeqsCollectionType) -> T:
        return self._func(seqs)


@define_app
class trim_stop_codons:
    """Removes terminal stop codons."""

    def __init__(self, gc=1):
        """
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
        self._gc = gc

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, data: SeqsCollectionType) -> T:
        data = data.trim_stop_codons(gc=self._gc)
        return data

from collections import defaultdict

from numpy import random as np_random, array

from cogent3.core.moltype import get_moltype
from cogent3.core.genetic_code import get_code
from cogent3.core.alignment import ArrayAlignment, Alignment
from .translate import get_fourfold_degenerate_sets
from .composable import ComposableSeq, ComposableAligned

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


# TODO need a function to filter sequences based on divergence, ala divergent
# set.

def intersection(groups):
    """returns the intersection of all groups"""
    common = set(groups.pop())
    intersect = common.intersection(*map(set, groups))
    return intersect


def union(groups):
    """returns the intersection of all groups"""
    union = set(groups.pop())
    union = union.union(*map(set, groups))
    return union


class concat:
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
        names = self._name_callback(list(aln.names for aln in data))
        collated = defaultdict(list)
        for aln in data:
            assert isinstance(aln, ArrayAlignment) or isinstance(aln,
                                                                 Alignment)
            if self._intersect:
                seqs = aln.take_seqs(names).todict()
            else:
                seqs = defaultdict(lambda: '?' * len(aln))
                seqs.update(aln.todict())

            for name in names:
                collated[name].append(seqs[name])

        combined = {n: self._join_seq.join(collated[n]) for n in names}
        aln = ArrayAlignment(data=combined, moltype=self._moltype)
        return aln

    __call__ = concat


class omit_degenerates(ComposableAligned):
    """returns alignment filtered by position"""

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
        super(omit_degenerates, self).__init__(input_type='aligned',
                                               output_type=('aligned',
                                                            'serialisable'))
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
            assert moltype.label.lower() in ('dna', 'rna'), "Invalid moltype"

        self.moltype = moltype
        self._no_degen = omit_degenerates
        self._allow_gap = not gap_is_degen
        self._motif_length = motif_length
        self.func = self.filter_degenerates

    def filter_degenerates(self, aln):
        if aln.moltype != self.moltype:
            # try converting
            aln = aln.to_type(moltype=self.moltype,
                              array_align=True)
        aln = aln.no_degenerates(motif_length=self._motif_length,
                                 allow_gap=self._allow_gap)
        return aln


class take_codon_positions(ComposableAligned):
    """returns the specified codon position(s) from an alignment"""

    def __init__(self, *positions, fourfold_degenerate=False,
                 gc='Standard Nuclear', moltype='dna'):
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
        super(take_codon_positions, self).__init__(input_type='aligned',
                                                   output_type=('aligned',
                                                                'serialisable'))
        self._formatted_params()
        assert moltype is not None
        moltype = get_moltype(moltype)

        assert moltype.label.lower() in ('dna', 'rna'), "Invalid moltype"

        self._moltype = moltype
        self._four_fold_degen = fourfold_degenerate
        self._fourfold_degen_sets = None

        if fourfold_degenerate:
            gc = get_code(gc)
            sets = get_fourfold_degenerate_sets(gc, alphabet=moltype.alphabet,
                                                as_indices=True)
            self._fourfold_degen_sets = sets
            self.func = self.take_fourfold_positions
            return

        assert 1 <= min(positions) <= 3 and 1 <= max(positions) <= 3, \
            'Invalid codon positions'

        by_index = True if len(positions) == 1 else False
        if by_index:
            positions = positions[0] - 1
            self.func = self.take_codon_position
        else:
            positions = tuple(p - 1 for p in sorted(positions))
            self.func = self.take_codon_positions

        self._positions = positions

    def take_fourfold_positions(self, aln):
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
        return aln[self._positions::3]

    def take_codon_positions(self, aln):
        '''takes multiple positions'''
        length = len(aln)
        indices = [k for k in range(length) if k % 3 in self._positions]
        return aln.take_positions(indices)


class take_named_seqs(ComposableSeq):
    def __init__(self, *names):
        """selects named sequences from a collection
        
        Returns
        -------
        A new sequence collection, or False if not all the named sequences are 
        in the collection.
        """
        super(take_named_seqs, self).__init__(input_type=('sequences',
                                                          'aligned'),
                                              output_type=('sequences',
                                                           'aligned',
                                                           'serialisable'))
        self._formatted_params()
        self._names = names
        self.func = self.take_seqs

    def take_seqs(self, data):
        try:
            data = data.take_seqs(self._names)
        except KeyError:
            data = False
        return data


class min_length(ComposableSeq):
    """filters sequence collections by length"""

    def __init__(self, length, motif_length=1, subtract_degen=True,
                 moltype=None):
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
        super(min_length, self).__init__(input_type=('sequences',
                                                     'aligned'),
                                         output_type=('sequences',
                                                      'aligned',
                                                      'serialisable'))
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
        if self._moltype:
            data = data.to_moltype(self._moltype)

        if self._subtract_degen:
            if not hasattr(data.alphabet, 'non_degen'):
                name = self.__class__.__name__
                msg = ('%s(subtract_degen=True) requires DNA, RNA or PROTEIN '
                       'moltype') % name
                raise ValueError(msg)

        lengths = data.get_lengths(allow_gap=not self._subtract_degen,
                                   include_ambiguity=not self._subtract_degen)
        length, _ = min([(l, n) for n, l in lengths.items()])

        if length < self._min_length:
            data = False

        return data


def _GetStart(start):
    choose = np_random.choice

    def _int(length):
        return start

    def _rand(length):
        return choose(length)

    func = {True: _int}.get(type(start) == int, _rand)
    return func


class fixed_length(ComposableAligned):
    """return alignments of a fixed length"""

    def __init__(self, length, start=0, random=False, seed=None,
                 motif_length=1, moltype=None):
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
        super(fixed_length, self).__init__(input_type='aligned',
                                           output_type=('aligned',
                                                        'serialisable'))
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
            assert start.lower().startswith('rand')
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
        if self._moltype:
            aln = aln.to_moltype(self._moltype)

        if len(aln) < self._length:
            result = False
        else:
            start = self._start(len(aln) - self._length)
            result = aln[start:start + self._length]

        return result

    def sample_positions(self, aln):
        if self._moltype:
            aln = aln.to_moltype(self._moltype)

        indices = array(range(len(aln)))
        if self._motif_length == 1:
            number = self._length
        else:
            number = self._length // self._motif_length

        pos = np_random.choice(indices.shape[0] // self._motif_length, number,
                               replace=False)
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
    def __init__(self, disallowed_frac=0.9, allowed_frac_bad_cols=0,
                 exclude_just_gap=True, moltype='dna'):
        """Returns an alignment without the sequences responsible for
        exceeding disallowed_frac.

        Parameters
        ----------
        disallowed_frac
            an aligned column with gap frac > disallowed_frac is considered a
            bad column
        allowed_frac_bad_cols : float
            the fraction of gaps induced by the sequence
        exclude_just_gap
            sequences that are just gaps are excluded
        moltype
            molecular type, can be string or instance
        """
        super(omit_bad_seqs, self).__init__(input_type='aligned',
                                            output_type=('aligned',
                                                         'serialisable'))
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        assert moltype.label.lower() in 'dna rna protein protein_with_stop', \
            'moltype must be one of DNA, RNA or PROTEIN'
        self._disallowed_frac = disallowed_frac
        self._allowed_frac_bad_col = allowed_frac_bad_cols
        self._exclude_just_gap = exclude_just_gap
        self._moltype = moltype
        self.func = self.drop_bad_seqs

    def drop_bad_seqs(self, aln):
        aln = aln.to_moltype(self._moltype)
        result = aln.omit_bad_seqs(disallowed_frac=self._disallowed_frac,
                                   allowed_frac_bad_cols=self._allowed_frac_bad_col,
                                   exclude_just_gap=self._exclude_just_gap)
        return result


class omit_duplicated(ComposableSeq):
    def __init__(self, mask_degen=False, choose='longest', seed=None,
                 moltype=None):
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
        super(omit_duplicated, self).__init__(input_type='sequences',
                                              output_type=('sequences',
                                                           'serialisable'))

        assert not choose or choose in 'longestrandom'
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self._moltype = moltype
        if choose == 'random' and seed:
            np_random.seed(seed)

        self._mask_degen = mask_degen
        if choose == 'longest':
            self.func = self.choose_longest
        elif choose == 'random':
            self.func = self.choose_random
        else:
            self.func = self.take_unique

    def choose_longest(self, seqs):
        if self._moltype:
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
        if self._moltype:
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
        if self._moltype:
            seqs = seqs.to_moltype(self._moltype)
        duplicates = seqs.get_identical_sets(mask_degen=self._mask_degen)
        names = set()
        for dupes in duplicates:
            names.update(dupes)
        seqs = seqs.take_seqs(names, negate=True)
        return seqs

#!/usr/bin/env python
"""Code for handling multiple sequence alignments. In particular:

    - SequenceCollection handles both aligned and unaligned sequences.
    - Alignment and its subclasses handle multiple sequence alignments, storing
      the raw sequences and a gap map. Useful for very long alignments, e.g.
      genomics data.
    - ArrayAlignment and its subclasses handle multiple sequence alignments as
      arrays of characters. Especially useful for short alignments that contain
      many sequences.

    WARNING: The various alignment objects try to guess the input type from the
    input, but this behavior has a few quirks. In particular, if the input is
    a sequence of two-item sequences (e.g. a list of two-character strings),
    each sequence will be unpacked and the first item will be used as the
    label, the second as the sequence. For example, Alignment(['AA','CC','AA'])
    produces an alignment of three 1-character strings labeled A, C and A
    respectively. The reason for this is that the common case is that you have
    passed in a stream of two-item label, sequence pairs. However, this can
    cause confusion when testing.
"""
import re
import json
from types import GeneratorType
from collections import defaultdict, Counter
from functools import total_ordering

from cogent3.core.annotation import Map, _Annotatable
import cogent3  # will use to get at cogent3.parse.fasta.MinimalFastaParser,
# which is a circular import otherwise.
from cogent3.format.alignment import save_to_filename
from cogent3.core.genetic_code import DEFAULT
from cogent3.core.info import Info as InfoClass
from cogent3.core.sequence import frac_same, ArraySequence
from cogent3.core.location import LostSpan, Span
from cogent3.maths.stats.util import Freqs
from cogent3.format.fasta import alignment_to_fasta
from cogent3.format.phylip import alignment_to_phylip
from cogent3.format.nexus import nexus_from_alignment
from cogent3.parse.gff import GffParser, parse_attributes
from numpy import nonzero, array, logical_or, logical_and, logical_not, \
    transpose, arange, zeros, ones, take, put, uint8, ndarray, vstack
from numpy.random import randint, permutation

from cogent3.util.dict_array import DictArrayTemplate
from cogent3.util.misc import bytes_to_string, get_object_provenance

from copy import copy, deepcopy
from cogent3.core.profile import Profile


__author__ = "Peter Maxwell and Rob Knight"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight", "Gavin Huttley",
               "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
               "Micah Hamady", "Daniel McDonald", "Jan Kosinski"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

_compression = re.compile(r"\.(gz|bz2)$")


class DataError(Exception):
    pass


# small number: 1-eps is almost 1, and is used for things like the
# default number of gaps to allow in a column.
eps = 1e-6

# factory functions for identifying whether character set conditions are satsified


def AllowedCharacters(chars, is_array=False, negate=False):
    """factory function for evaluating whether a sequence contains only
    the specified gap characters"""
    try:
        chars = set(chars)
    except TypeError:
        chars = set([chars])

    def just_chars(data):
        if is_array:
            data = set(data.flatten())
        else:
            data = set("".join(data))
        return data <= chars

    def not_chars(data):
        if is_array:
            data = set(data.flatten())
        else:
            data = set("".join(data))
        return not data & chars

    if negate:
        return not_chars

    return just_chars


def GapsOk(gap_chars, allowed_frac, motif_length=1, is_array=False, negate=False):
    """factory function for determining whether sequence contains a
    fraction of gaps <= allowed_frac"""
    try:
        gap_chars = set(gap_chars)
    except TypeError:
        gap_chars = set([gap_chars])

    def gap_frac_ok(data):
        length = len(data) / motif_length
        if is_array:
            data = Counter(data.flatten())
        else:
            data = Counter(list(data))

        num_gap = sum(data[g] for g in gap_chars)
        return num_gap / length <= allowed_frac

    def gap_frac_not_ok(data):
        """returns True when the frac of gaps is > allowed_frac"""
        data = array(data)
        length = len(data) / motif_length
        if is_array:
            data = Counter(data.flatten())
        else:
            data = Counter(list(data))

        num_gap = sum(data[g] for g in gap_chars)
        return (num_gap / length) >= allowed_frac

    if negate:
        return gap_frac_not_ok
    else:
        return gap_frac_ok


def assign_sequential_names(ignored, num_seqs, base_name='seq', start_at=0):
    """Returns list of num_seqs sequential, unique names.

    First argument is ignored; expect this to be set as a class attribute.
    """
    return ['%s_%s' % (base_name, i) for i in range(start_at, start_at + num_seqs)]


class SeqLabeler(object):
    """Allows flexible seq labeling in to_fasta()."""

    def __init__(self, aln, label_f=assign_sequential_names, **kwargs):
        """Initializes a new seq labeler."""
        self._aln = aln
        self._label_f = label_f
        self._map = dict(
            list(zip(aln.names, label_f(len(aln.names, **kwargs)))))

    def __call__(self, s):
        """Returns seq name from seq id"""
        return self._map[s.name]


def coerce_to_string(s):
    """Converts an arbitrary sequence into a string."""
    if isinstance(s, str):  # if it's a string, OK as is
        return s
    if isinstance(s, Aligned):  # if it's an Aligned object, convert to string
        return str(s)
    curr = str(s)  # if its string is the same length, return that
    if len(curr) == len(s):
        return curr
    try:
        return ''.join(s)  # assume it's a seq of chars
    except(TypeError, ValueError):
        # general case (slow, might not be correct)
        return ''.join(map(str, s))


def seqs_from_array(a, alphabet=None):
    """SequenceCollection from array of pos x seq: names are integers.

    This is an InputHandler for SequenceCollection. It converts an arbitrary
    array of numbers into Sequence objects, and leaves the sequences unlabeled.
    """
    return list(transpose(a)), None


def seqs_from_array_seqs(seqs, alphabet=None):
    """Alignment from ArraySequence objects: seqs -> array, names from seqs.

    This is an InputHandler for SequenceCollection. It converts a list of
    Sequence objects with _data and name properties into a SequenceCollection
    that uses those sequences.
    """
    return seqs, [s.name for s in seqs]


def seqs_from_generic(seqs, alphabet=None):
    """returns seqs, names
    """
    names = []
    for s in seqs:
        if hasattr(s, 'name'):
            names.append(s.name)
        else:
            names.append(None)
    return seqs, names


def seqs_from_fasta(seqs, alphabet=None):
    """SequenceCollection from FASTA-format string or lines.

    This is an InputHandler for SequenceCollection. It converts a FASTA-format
    string or collection of lines into a SequenceCollection object, preserving
    order..
    """
    if isinstance(seqs, str):
        seqs = seqs.splitlines()
    names, seqs = list(
        zip(*list(cogent3.parse.fasta.MinimalFastaParser(seqs))))
    return list(seqs), list(names)


def seqs_from_dict(seqs, alphabet=None):
    """SequenceCollection from dict of {label:seq_as_str}.

    This is an InputHandler for SequenceCollection. It converts a dict in
    which the keys are the names and the values are the sequences
    (sequence only, no whitespace or other formatting) into a
    SequenceCollection. Because the dict doesn't preserve order, the result
    will not necessarily be in alphabetical order."""
    names, seqs = list(map(list, list(zip(*list(seqs.items())))))
    seqs = [bytes_to_string(s) for s in seqs]
    return seqs, names


def seqs_from_kv_pairs(seqs, alphabet=None):
    """SequenceCollection from list of (key, val) pairs.

    This is an InputHandler for SequenceCollection. It converts a dict in
    which the keys are the names and the values are the sequences
    (sequence only, no whitespace or other formatting) into a
    SequenceCollection. Because the dict doesn't preserve order, the result
    will be in arbitrary order."""
    names, seqs = list(map(list, list(zip(*seqs))))
    seqs = [bytes_to_string(s) for s in seqs]
    return seqs, names


def seqs_from_aln(seqs, alphabet=None):
    """SequenceCollection from existing SequenceCollection object: copies data.

    This is relatively inefficient: you should really use the copy() method
    instead, which duplicates the internal data structures.
    """
    return seqs.seqs, seqs.names


def seqs_from_empty(obj, *args, **kwargs):
    """SequenceCollection from empty data: raise exception."""
    raise ValueError("Cannot create empty SequenceCollection.")


@total_ordering
class SequenceCollection(object):
    """base class for Alignment, but also just stores unaligned seqs.

    Handles shared functionality: detecting the input type, writing out the
    sequences as different formats, translating the sequences, chopping off
    stop codons, looking up sequences by name, etc.

    A SequenceCollection must support:

    - input handlers for different data types
    - seq_data: behaves like list of lists of chars, holds seq data
    - seqs: behaves like list of Sequence objects, iterable in name order
    - names: behaves like list of names for the sequence objects
    - named_seqs: behaves like dict of {name:seq}
    - moltype: specifies what kind of sequences are in the collection
    """
    InputHandlers = {'array': seqs_from_array,
                     'array_seqs': seqs_from_array_seqs,
                     'generic': seqs_from_generic,
                     'fasta': seqs_from_fasta,
                     'collection': seqs_from_aln,
                     'aln': seqs_from_aln,
                     'array_aln': seqs_from_aln,
                     'dict': seqs_from_dict,
                     'empty': seqs_from_empty,
                     'kv_pairs': seqs_from_kv_pairs,
                     }

    IsArray = set(['array', 'array_seqs'])

    DefaultNameFunction = assign_sequential_names

    def __init__(self, data, names=None, alphabet=None, moltype=None,
                 name=None, info=None, conversion_f=None, is_array=False,
                 force_same_data=False,
                 remove_duplicate_names=False, label_to_name=None,
                 suppress_named_seqs=False):
        """Initialize self with data and optionally info.

        We are always going to convert to characters, so Sequence objects
        in the collection will lose additional special attributes they have.
        This is somewhat inefficient, so it might be worth revisiting this
        decision later.

        The handling of sequence names requires special attention. Depending
        on the input data, we might get the names from the sequences themselves,
        or we might add them from Names that are passed in. However, the Names
        attribute controls the order that we examine the sequences in, so if
        it is passed in it should override the order that we got from the
        input data (e.g. you might pass in unlabeled sequences with the names
        ['b','a'], so that you want the first sequence to be called 'b' and
        the second to be called 'a', or you might pass in labeled sequences,
        e.g. as a dict, and the names ['b','a'], indicating that you want the
        sequence called b to be first and the sequence called a to be second
        despite the fact that they are in arbitrary order in the original
        input. In this second situation, it is imortant that the sequences not
        be relabeled.

        This is handled as followed. If the sequences are passed in using a
        method that does not carry the names with it, the Names that are passed
        in will be handed out to successive sequences. If the sequences are
        passed in using a method that does carry the names with it, the Names
        that are passed in will be used to order the sequences, but they will
        not be relabeled. Note that if you're passing in a data type that
        is already labeled (e.g. a list of Sequence objects) you _must_ have
        unique names beforehand.

        It's possible that this additional handling should be moved to a
        separate object; the motivation for having it on Alignment __init__
        is that it's easy for users to construct Alignment objects directly.

        Parameters:

        data:           Data to convert into a SequenceCollection

        Names:          Order of Names in the alignment. Should match the
                        names of the sequences (after processing by
                        label_to_name if present).

        alphabet:       Alphabet to use for the alignment (primarily important
                        for ArrayAlignment)

        moltype:        moltype to be applied to the Alignment and to each seq.

        name:           name of the SequenceCollection.

        info:           info object to be attached to the alignment itself.

        conversion_f:   Function to convert string into sequence.

        is_array:       True if input is an array, False otherwise.

        force_same_data: True if data will be used as the same object.

        remove_duplicate_names: True if duplicate names are to be silently
                        deleted instead of raising errors.

        label_to_name: if present, converts name into f(name).
        """

        # read all the data in if we were passed a generator
        if isinstance(data, GeneratorType):
            data = list(data)
        # set the name
        self.name = name
        # figure out alphabet and moltype
        self.alphabet, self.moltype = \
            self._get_alphabet_and_moltype(alphabet, moltype, data)
        if not isinstance(info, InfoClass):
            if info:
                info = InfoClass(info)
            else:
                info = InfoClass()
        self.info = info
        # if we're forcing the same data, skip the validation
        if force_same_data:
            self._force_same_data(data, names)
            curr_seqs = data
        # otherwise, figure out what we got and coerce it into the right type
        else:
            per_seq_names, curr_seqs, name_order = \
                self._names_seqs_order(conversion_f, data, names, is_array,
                                       label_to_name, remove_duplicate_names,
                                       alphabet=self.alphabet)
            self.names = list(name_order)

            # will take only the seqs and names that are in name_order
            if per_seq_names != name_order:
                good_indices = []
                for n in name_order:
                    good_indices.append(per_seq_names.index(n))
                if hasattr(curr_seqs, 'astype'):  # it's an array
                    # much faster to check than to raise exception in this case
                    curr_seqs = take(curr_seqs, good_indices, axis=0)
                else:
                    curr_seqs = [curr_seqs[i] for i in good_indices]
                per_seq_names = name_order

            # create named_seqs dict for fast lookups
            if not suppress_named_seqs:
                self.named_seqs = self._make_named_seqs(
                    self.names, curr_seqs)
        # Sequence objects behave like sequences of chars, so no difference
        # between seqs and seq_data. Note that this differs for Alignments,
        # so be careful which you use if writing methods that should work for
        # both SequenceCollections and Alignments.
        self._set_additional_attributes(curr_seqs)

    def __str__(self):
        """Returns self in FASTA-format, respecting name order."""
        return ''.join(['>%s\n%s\n' % (name, self.get_gapped_seq(name))
                        for name in self.names])

    def _make_named_seqs(self, names, seqs):
        """Returns named_seqs: dict of name:seq."""
        name_seq_tuples = list(zip(names, seqs))
        for n, s in name_seq_tuples:
            s.name = n
        return dict(name_seq_tuples)

    def _set_additional_attributes(self, curr_seqs):
        """Sets additional attributes based on current seqs: class-specific."""
        self.seq_data = curr_seqs
        self._seqs = curr_seqs
        try:
            self.seq_len = max(list(map(len, curr_seqs)))
        except ValueError:  # got empty sequence, for some reason?
            self.seq_len = 0

    def _force_same_data(self, data, names):
        """Forces dict that was passed in to be used as self.named_seqs"""
        self.named_seqs = data
        self.names = names or list(data.keys())

    def copy(self):
        """Returns deep copy of self."""
        result = self.__class__(self, moltype=self.moltype, info=self.info)
        return result

    def _get_alphabet_and_moltype(self, alphabet, moltype, data):
        """Returns alphabet and moltype, giving moltype precedence."""
        if type(moltype) == str:
            from cogent3.core.moltype import get_moltype
            moltype = get_moltype(moltype)

        if alphabet is None and moltype is None:
            if hasattr(data, 'moltype'):
                moltype = data.moltype
            elif hasattr(data, 'alphabet'):
                alphabet = data.alphabet
            # check for containers
            else:
                curr_item = self._get_container_item(data)
                if hasattr(curr_item, 'moltype'):
                    moltype = curr_item.moltype
                elif hasattr(curr_item, 'alphabet'):
                    alphabet = curr_item.alphabet
                else:
                    moltype = self.moltype  # will be BYTES by default
        if alphabet is not None and moltype is None:
            moltype = alphabet.moltype
        if moltype is not None and alphabet is None:
            try:
                alphabet = moltype.alphabets.degen_gapped
            except AttributeError:
                alphabet = moltype.alphabet
        return alphabet, moltype

    def _get_container_item(self, data):
        """Checks container for item with alphabet or moltype"""
        curr_item = None
        if hasattr(data, 'values'):
            curr_item = next(iter(data.values()))
        else:
            try:
                curr_item = next(iter(data))
            except:
                pass
        return curr_item

    def _strip_duplicates(self, names, seqs):
        """Internal function to strip duplicates from list of names"""
        if len(set(names)) == len(names):
            return set(), names, seqs
        # if we got here, there are duplicates
        unique_names = {}
        duplicates = {}
        fixed_names = []
        fixed_seqs = []
        for n, s in zip(names, seqs):
            if n in unique_names:
                duplicates[n] = 1
            else:
                unique_names[n] = 1
                fixed_names.append(n)
                fixed_seqs.append(s)
        if type(seqs) is ndarray:
            fixed_seqs = array(fixed_seqs, seqs.dtype)
        return duplicates, fixed_names, fixed_seqs

    def _names_seqs_order(self, conversion_f, data, name_order, is_array,
                          label_to_name, remove_duplicate_names, alphabet=None):
        """Internal function to figure out names, seqs, and name_order."""
        # figure out conversion function and whether it's an array
        if not conversion_f:
            input_type = self._guess_input_type(data)
            is_array = input_type in self.IsArray
            conversion_f = self.InputHandlers[input_type]
        # set seqs and names as properties
        if alphabet:
            seqs, names = conversion_f(data, alphabet=alphabet)
        else:
            seqs, names = conversion_f(data)
        if names and label_to_name:
            names = list(map(label_to_name, names))
        curr_seqs = self._coerce_seqs(seqs, is_array)

        # if no names were passed in as Names, if we obtained them from
        # the seqs we should use them, but otherwise we should use the
        # default names
        if name_order is None:
            if (names is None) or (None in names):
                per_seq_names = name_order = \
                    self.DefaultNameFunction(len(curr_seqs))
            else:  # got names from seqs
                per_seq_names = name_order = names
        else:
            # otherwise, names were passed in as Names: use this as the order
            # if we got names from the sequences, but otherwise assign the
            # names to successive sequences in order
            if (names is None) or (None in names):
                per_seq_names = name_order = name_order
            else:  # got names from seqs, so assume name_order is in Names
                per_seq_names = names
                name_order = name_order

        # check for duplicate names
        duplicates, fixed_names, fixed_seqs = \
            self._strip_duplicates(per_seq_names, curr_seqs)
        if duplicates:
            if remove_duplicate_names:
                per_seq_names, curr_seqs = fixed_names, fixed_seqs
                # if name_order doesn't have the same names as per_seq_names,
                # replace it with per_seq_names
                if (set(name_order) != set(per_seq_names)) or\
                        (len(name_order) != len(per_seq_names)):
                    name_order = per_seq_names
            else:
                raise ValueError("Some names were not unique. Duplicates are:\n" +
                                 str(sorted(duplicates.keys())))

        return per_seq_names, curr_seqs, name_order

    def _coerce_seqs(self, seqs, is_array):
        """Controls how seqs are coerced in _names_seqs_order.

        Override in subclasses where this behavior should differ.
        """
        if is_array:
            seqs = list(
                map(str, list(map(self.moltype.make_array_seq, seqs))))
        return list(map(self.moltype.make_seq, seqs))

    def _guess_input_type(self, data):
        """Guesses input type of data; returns result as key of InputHandlers.

        First checks whether data is an Alignment, then checks for some common
        string formats, then tries to do it based on string or array properties.

        Returns 'empty' if check fails, i.e. if it can't recognize the sequence
        as a specific type. Note that bad sequences are not guaranteed to
        return 'empty', and may be recognized as another type incorrectly.
        """
        if isinstance(data, ArrayAlignment):
            return 'array_aln'
        if isinstance(data, Alignment):
            return 'aln'
        if isinstance(data, SequenceCollection):
            return 'collection'
        if isinstance(data, dict):
            return 'dict'
        if isinstance(data, str):
            if data.startswith('>'):
                return 'fasta'
            else:
                return 'generic'
        first = None
        try:
            first = data[0]
        except (IndexError, TypeError):
            pass
        try:
            first = next(iter(data))
        except (IndexError, TypeError, StopIteration):
            pass
        if first is None:
            return 'empty'
        try:
            if isinstance(first, ArraySequence):  # model sequence base type
                return 'array_seqs'
            elif hasattr(first, 'dtype'):  # array object
                return 'array'
            elif isinstance(first, str) and first.startswith('>'):
                return 'fasta'
            else:
                try:
                    dict(data)
                    return 'kv_pairs'
                except (TypeError, ValueError):
                    pass
            return 'generic'
        except (IndexError, TypeError) as e:
            return 'empty'

    def __eq__(self, other):
        """first tests as dict, then as str"""
        c = self.named_seqs == other
        if not c:
            c = str(self) == str(other)

        return c

    def __ne__(self, other):
        """first tests as dict, then as str"""
        c = self.named_seqs != other
        if not c:
            c = str(self) != str(other)

        return c

    def __lt__(self, other):
        """cmp first tests as dict, then as str."""
        c = self.named_seqs < other
        if not c:
            return 0
        else:
            return str(self) < str(other)

    def keys(self):
        """keys uses self.Names, which defaults to known keys if None.

        Note: returns copy, not original.
        """
        return self.names[:]

    def values(self):
        """values returns values corresponding to self.Names."""
        return [self.named_seqs[n] for n in self.names]

    def items(self):
        """items returns (name, value) pairs."""
        return [(n, self.named_seqs[n]) for n in self.names]

    def iter_seqs(self, seq_order=None):
        """Iterates over values (sequences) in the alignment, in order.

        seq_order: list of keys giving the order in which seqs will be returned.
        Defaults to self.Names. Note that only these sequences will be
        returned, and that KeyError will be raised if there are sequences
        in order that have been deleted from the Alignment. If self.Names
        is None, returns the sequences in the same order as
        self.named_seqs.values().

        Use map(f, self.seqs()) to apply the constructor f to each seq. f must
        accept a single list as an argument.

        Always returns references to the same objects that are values of the
        alignment.
        """
        ns = self.named_seqs
        get = ns.__getitem__
        for key in seq_order or self.names:
            yield get(key)

    def _take_seqs(self):
        return list(self.iter_seqs())

    # access as attribute if using default order.
    seqs = property(_take_seqs)

    def take_seqs(self, seqs, negate=False, **kwargs):
        """Returns new Alignment containing only specified seqs.

        Note that the seqs in the new alignment will be references to the
        same objects as the seqs in the old alignment.
        """
        get = self.named_seqs.__getitem__
        result = {}
        if 'moltype' not in kwargs:
            kwargs['moltype'] = self.moltype

        if negate:
            # copy everything except the specified seqs
            negated_names = []
            row_lookup = dict.fromkeys(seqs)
            for r, row in list(self.named_seqs.items()):
                if r not in row_lookup:
                    result[r] = row
                    negated_names.append(r)
            seqs = negated_names  # remember to invert the list of names

        else:
            # copy only the specified seqs
            for r in seqs:
                result[r] = get(r)
        if result:
            return self.__class__(result, names=seqs, info=self.info, **kwargs)
        else:
            return {}  # safe value; can't construct empty alignment

    def get_seq_indices(self, f, negate=False):
        """Returns list of keys of seqs where f(row) is True.

        List will be in the same order as self.names, if present.
        """
        get = self.named_seqs.__getitem__
        # negate function if necessary
        if negate:
            def new_f(x): return not f(x)
        else:
            new_f = f
        # get all the seqs where the function is True
        return [key for key in self.names if new_f(get(key))]

    def take_seqs_if(self, f, negate=False, **kwargs):
        """Returns new Alignment containing seqs where f(row) is True.

        Note that the seqs in the new Alignment are the same objects as the
        seqs in the old Alignment, not copies.
        """
        # pass negate to get SeqIndices
        return self.take_seqs(self.get_seq_indices(f, negate), **kwargs)

    def iter_selected(self, seq_order=None, pos_order=None):
        """Iterates over elements in the alignment.

        seq_order (names) can be used to select a subset of seqs.
        pos_order (positions) can be used to select a subset of positions.

        Always iterates along a seq first, then down a position (transposes
        normal order of a[i][j]; possibly, this should change)..

        WARNING: Alignment.iter_selected() is not the same as alignment.iteritems()
        (which is the built-in dict iteritems that iterates over key-value
        pairs).
        """
        if pos_order:
            for row in self.iter_seqs(seq_order):
                for i in pos_order:
                    yield row[i]
        else:
            for row in self.iter_seqs(seq_order):
                for i in row:
                    yield i

    def get_items(self, items, negate=False):
        """Returns list containing only specified items.

        items should be a list of (row_key, col_key) tuples.
        """
        get = self.named_seqs.__getitem__
        if negate:
            # have to cycle through every item and check that it's not in
            # the list of items to return
            item_lookup = dict.fromkeys(list(map(tuple, items)))
            result = []
            for r in self.names:
                curr_row = get(r)
                for c in range(len(curr_row)):
                    if (r, c) not in items:
                        result.append(curr_row[c])
            return result
        # otherwise, just pick the selected items out of the list
        else:
            return [get(row)[col] for row, col in items]

    def item_indices_if(self, f, negate=False):
        """Returns list of (key,val) tuples where f(self.named_seqs[key][val]) is True"""
        get = self.named_seqs.__getitem__
        if negate:
            def new_f(x): return not f(x)
        else:
            new_f = f
        result = []
        for row_label in self.names:
            curr_row = get(row_label)
            for col_idx, item in enumerate(curr_row):
                if new_f(item):
                    result.append((row_label, col_idx))
        return result

    def items_if(self, f, negate=False):
        """Returns list of items where f(self.named_seqs[row][col]) is True."""
        return self.get_items(self.item_indices_if(f, negate))

    def get_identical_sets(self, mask_degen=False):
        """returns sets of names for sequences that are identical

        Arguments:
         - mask_degen: if True, degenerate characters are ignored"""
        # TODO handle case of not strict by building mask of degen positions
        # per seq
        if mask_degen and not hasattr(self.moltype, "alphabets"):
            UserWarning("in get_identical_sets, strict has no effect as moltype "
                        "has no degenerate characters")
            mask_degen = False
        elif mask_degen:
            degens = list(self.moltype.degenerates) + [self.moltype.gap]

        def reduced(seq, indices):
            s = "".join(seq[i] for i in range(len(seq)) if i not in indices)
            return s

        identical_sets = []
        mask_posns = {}
        seen = []
        # if strict, we do a sort and one pass through the list
        seqs = self.todict()
        if not mask_degen:
            seqs_names = [(s, n) for n, s in seqs.items()]
            seqs_names.sort()
            matched = None
            dupes = defaultdict(set)
            for i in range(len(seqs_names) - 1):
                if seqs_names[i][0] == seqs_names[i + 1][0]:
                    matched = seqs_names[i][1] if matched is None else matched
                    dupes[matched].update([seqs_names[i + 1][1], matched])
                else:
                    matched = None
            identical_sets = list(dupes.values())
            return identical_sets

        for i in range(len(self.names) - 1):
            n1 = self.names[i]
            if n1 in seen:
                continue

            seq1 = seqs[n1]
            if n1 not in mask_posns:
                pos = self.moltype.get_degenerate_positions(seq1,
                                                            include_gap=True)
                mask_posns[n1] = pos
                # TODO separately implement for ArrayAlignment

            group = set()
            for j in range(i + 1, len(self.names)):
                n2 = self.names[j]
                if n2 in seen:
                    continue

                seq2 = seqs[n2]
                if n2 not in mask_posns:
                    pos = self.moltype.get_degenerate_positions(seq2,
                                                                include_gap=True)
                    mask_posns[n2] = pos

                seq2 = seqs[n2]

                s1, s2 = seq1, seq2

                pos = mask_posns[n1] + mask_posns[n2]
                if pos:
                    s1 = reduced(seq1, pos)
                    s2 = reduced(seq2, pos)

                if s1 == s2:
                    seen.append(n2)
                    group.update([n1, n2])

            if group:
                identical_sets.append(group)

        return identical_sets

    def get_similar(self, target, min_similarity=0.0, max_similarity=1.0,
                    metric=frac_same, transform=None):
        """Returns new Alignment containing sequences similar to target.

        target: sequence object to compare to. Can be in the alignment.

        min_similarity: minimum similarity that will be kept. Default 0.0.

        max_similarity: maximum similarity that will be kept. Default 1.0.
        (Note that both min_similarity and max_similarity are inclusive.)

        metric: similarity function to use. Must be f(first_seq, second_seq).
        The default metric is fraction similarity, ranging from 0.0 (0%
        identical) to 1.0 (100% identical). The Sequence classes have lots
        of methods that can be passed in as unbound methods to act as the
        metric, e.g. frac_same_gaps.

        transform: transformation function to use on the sequences before
        the metric is calculated. If None, uses the whole sequences in each
        case. A frequent transformation is a function that returns a specified
        range of a sequence, e.g. eliminating the ends. Note that the
        transform applies to both the real sequence and the target sequence.

        WARNING: if the transformation changes the type of the sequence (e.g.
        extracting a string from an RnaSequence object), distance metrics that
        depend on instance data of the original class may fail.
        """
        if transform:
            target = transform(target)

        def m(x): return metric(target, x)

        if transform:
            def f(x):
                result = m(transform(x))
                return min_similarity <= result <= max_similarity
        else:
            def f(x):
                result = m(x)
                return min_similarity <= result <= max_similarity

        return self.take_seqs_if(f)

    def distance_matrix(self, calc='hamming'):
        """Returns pairwise distances between sequences.
        Parameters
        ----------
        name
            a pairwise distance calculator or name of one from pairwise_distance
        """
        from cogent3.evolve.pairwise_distance import get_calculator
        calculator = get_calculator(
            calc, moltype=self.moltype, alignment=self)
        calculator.run()
        result = calculator.get_pairwise_distances()
        return result

    def is_ragged(self):
        """Returns True if alignment has sequences of different lengths."""
        seqs = self.seqs  # Get all sequences in alignment
        length = len(seqs[0])  # Get length of first sequence
        for seq in seqs:
            # If lengths differ
            if length != len(seq):
                return True
        # lengths were all equal
        return False

    def to_phylip(self, generic_label=True, make_seqlabel=None):
        """
        Return alignment in PHYLIP format and mapping to sequence ids

        raises exception if invalid alignment

        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str
        """
        return alignment_to_phylip(self.todict())

    def to_rich_dict(self):
        """returns detailed content including info and moltype attributes"""
        data = {}
        moltype = self.moltype.label
        info = {} if self.info is None else self.info
        if not info.get('Refs', None) is None and 'Refs' in info:
            info.pop('Refs')

        info = info or None

        for seq in self.seqs:
            d = seq.to_rich_dict()
            for attr in ['moltype', 'type']:
                del(d[attr])
            data[seq.name] = d

        data = dict(seqs=data, moltype=moltype, info=info,
                    type=get_object_provenance(self))
        return data

    def to_json(self):
        """returns json formatted string"""
        return json.dumps(self.to_rich_dict())

    def to_fasta(self):
        """Return alignment in Fasta format

        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str
        """
        return alignment_to_fasta(self.todict())

    def to_nexus(self, seq_type, interleave_len=50):
        """
        Return alignment in NEXUS format and mapping to sequence ids

        **NOTE** Not that every sequence in the alignment MUST come from
            a different species!! (You can concatenate multiple sequences from
            same species together before building tree)

        seq_type: dna, rna, or protein

        Raises exception if invalid alignment
        """
        return nexus_from_alignment(self, seq_type,
                                    interleave_len=interleave_len)

    def get_int_map(self, prefix='seq_'):
        """Returns a dict with names mapped to enumerates integer names.

            - prefix: prefix for sequence label. Default = 'seq_'
            - int_keys is a dict mapping int names to sorted original names.
        """
        get = self.named_seqs.__getitem__
        int_keys = dict([(prefix + str(i), k) for i, k in
                         enumerate(sorted(self.named_seqs.keys()))])
        int_map = dict([(k, copy(get(v)))
                        for k, v in list(int_keys.items())])
        return int_map, int_keys

    @property
    def num_seqs(self):
        """Returns the number of sequences in the alignment."""
        return len(self.named_seqs)

    def copy_annotations(self, unaligned):
        """Copies annotations from seqs in unaligned to self, matching by name.

        Alignment programs like ClustalW don't preserve annotations,
        so this method is available to copy annotations off the unaligned
        sequences.

        unaligned should be a dictionary of Sequence instances.

        Ignores sequences that are not in self, so safe to use on larger dict
        of seqs that are not in the current collection/alignment.
        """
        for name, seq in list(unaligned.items()):
            if name in self.named_seqs:
                self.named_seqs[name].copy_annotations(seq)

    def annotate_from_gff(self, f):
        """Copies annotations from gff-format file to self.

        Matches by name of sequence. This method expects a file handle, not
        the name of a file.

        Skips sequences in the file that are not in self.
        """
        for (name, source, feature, start, end, score,
                strand, frame, attributes, comments) in GffParser(f):
            if name in self.named_seqs:
                self.named_seqs[name].add_feature(feature,
                                                  parse_attributes(
                                                      attributes),
                                                  [(start, end)])

    def get_gapped_seq(self, seq_name, recode_gaps=False, moltype=None):
        """Return a gapped Sequence object for the specified seqname.

        Note: return type may depend on what data was loaded into the
        SequenceCollection or Alignment.
        """
        return self.named_seqs[seq_name]

    def __add__(self, other):
        """Concatenates sequence data for same names"""
        aligned = isinstance(self, Alignment)

        if len(self.named_seqs) != len(other.named_seqs):
            raise ValueError(
                "Alignments don't have same number of sequences")

        concatenated = []
        for name in self.names:
            if name not in other.names:
                raise ValueError(
                    "Right alignment doesn't have a '%s'" % name)
            new_seq = self.named_seqs[name] + other.named_seqs[name]
            concatenated.append(new_seq)

        new = self.__class__(moltype=self.moltype,
                             data=list(zip(self.names, concatenated)), info=self.info)

        if aligned:
            left = [a for a in self._shifted_annotations(new, 0)
                    if a.map.end <= len(self)]
            right = [a for a in other._shifted_annotations(new, len(self))
                     if a.map.start >= len(self)]
            new.annotations = left + right
        return new

    def add_seqs(self, other, before_name=None, after_name=None):
        """Returns new object of class self with sequences from other added.

        By default the sequence is appended to the end of the alignment,
        this can be changed by using either before_name or after_name arguments.

        Arguments:
            - other: same class as self or coerceable to that class
            - before_name: str - [default:None] name of the sequence before
              which sequence is added
            - after_name: str - [default:None] name of the sequence after
              which sequence is added

        If both before_name and after_name are specified, the seqs will be
        inserted using before_name.
        """
        assert not isinstance(other, str), "Must provide a series of seqs " +\
            "or an alignment"
        self_seq_class = self.seqs[0].__class__
        try:
            combined = self.seqs + other.seqs
        except AttributeError:
            combined = self.seqs + list(other)

        for seq in combined:
            assert seq.__class__ == self_seq_class,\
                "Seq classes different: Expected %s, Got %s" % \
                (seq.__class__, self_seq_class)

        combined_aln = self.__class__(data=combined, info=self.info)

        if before_name is None and after_name is None:
            return combined_aln

        if (before_name and before_name not in self.names) \
                or \
                (after_name and after_name not in self.names):
            name = before_name or after_name
            raise ValueError("The alignment doesn't have a sequence named '{0}'"
                             .format(name))

        if before_name is not None:  # someone might have seqname of int(0)
            index = self.names.index(before_name)
        elif after_name is not None:
            index = self.names.index(after_name) + 1

        names_before = self.names[:index]
        names_after = self.names[index:]
        new_names = combined_aln.names[len(self.names):]

        aln_new = combined_aln.take_seqs(new_names)
        if len(names_before) > 0:
            aln_before = self.take_seqs(names_before)
            combined_aln = aln_before
            combined_aln = combined_aln.add_seqs(aln_new)
        else:
            combined_aln = aln_new

        if len(names_after) > 0:
            aln_after = self.take_seqs(names_after)
            combined_aln = combined_aln.add_seqs(aln_after)

        return combined_aln

    def write(self, filename=None, format=None, **kwargs):
        """Write the alignment to a file, preserving order of sequences.

        Arguments:
        - filename: name of the sequence file
        - format: format of the sequence file

        If format is None, will attempt to infer format from the filename
        suffix.
        """

        if filename is None:
            raise DataError('no filename specified')

        # need to turn the alignment into a dictionary
        align_dict = {}
        for seq_name in self.names:
            align_dict[seq_name] = str(self.named_seqs[seq_name])

        if format is None and '.' in filename:
            # allow extension to work if provided
            r = _compression.search(filename)
            if r:
                format = filename[:r.start()].split('.')[-1]
            else:
                format = filename.split('.')[-1]

        if 'order' not in kwargs:
            kwargs['order'] = self.names
        save_to_filename(align_dict, filename, format, **kwargs)

    def __len__(self):
        """len of SequenceCollection returns length of longest sequence."""
        return self.seq_len

    def get_translation(self, gc=None, **kwargs):
        """Returns a new alignment object with the DNA sequences translated,
        using the current codon moltype, into an amino acid sequence.
        """
        translated = []
        aligned = isinstance(self, Alignment)
        kwargs['moltype'] = cogent3.PROTEIN
        # do the translation
        try:
            for seqname in self.names:
                if aligned:
                    seq = self.get_gapped_seq(seqname)
                else:
                    seq = self.named_seqs[seqname]
                pep = seq.get_translation(gc)
                translated.append((seqname, pep))
            return self.__class__(translated, info=self.info, **kwargs)
        except AttributeError as msg:
            raise AttributeError("%s -- %s" %
                                 (msg, "Did you set a DNA moltype?"))

    def get_seq(self, seqname):
        """Return a sequence object for the specified seqname.
        """
        return self.named_seqs[seqname]

    def todict(self):
        """Returns the alignment as dict of names -> strings.

        Note: returns strings, NOT Sequence objects.
        """
        align_dict = {}

        for seq_name in self.names:
            align_dict[seq_name] = str(self.named_seqs[seq_name])

        return align_dict

    def get_ambiguous_positions(self):
        """Returns dict of seq:{position:char} for ambiguous chars.

        Used in likelihood calculations.
        """
        result = {}
        for name in self.names:
            result[name] = ambig = {}
            for (i, motif) in enumerate(self.get_gapped_seq(name)):
                if self.moltype.is_ambiguity(motif):
                    ambig[i] = motif
        return result

    def degap(self, **kwargs):
        """Returns copy in which sequences have no gaps."""
        new_seqs = []
        aligned = isinstance(self, Alignment)
        for seq_name in self.names:
            if aligned:
                seq = self.named_seqs[seq_name].data
            else:
                seq = self.named_seqs[seq_name]
            new_seqs.append((seq_name, seq.degap()))
        result = SequenceCollection(moltype=self.moltype, data=new_seqs,
                                    info=self.info, **kwargs)
        return result

    def with_modified_termini(self):
        """Changes the termini to include termini char instead of gapmotif.

        Useful to correct the standard gap char output by most
        alignment programs when aligned sequences have different ends.
        """
        seqs = []
        for name in self.names:
            seq = self.named_seqs[name].with_termini_unknown()
            seqs.append((name, seq))
        return self.__class__(moltype=self.moltype, data=seqs, info=self.info)

    def has_terminal_stops(self, gc=None, allow_partial=False):
        """Returns True if any sequence has a terminal stop codon.

        Arguments:
            - gc: genetic code object
            - allow_partial: if True and the sequence length is not divisible
              by 3, ignores the 3' terminal incomplete codon
        """
        stops = []
        aligned = isinstance(self, Alignment)
        for seq_name in self.names:
            if aligned:
                seq = self.named_seqs[seq_name].data
            else:
                seq = self.named_seqs[seq_name]
            stops.append(seq.has_terminal_stop(
                gc=gc, allow_partial=allow_partial))
        return max(stops)

    def trim_stop_codons(self, gc=None, allow_partial=False, **kwargs):
        """Removes any terminal stop codons from the sequences

        Arguments:
            - gc: genetic code object
            - allow_partial: if True and the sequence length is not divisible
              by 3, ignores the 3' terminal incomplete codon
        """
        new_seqs = []
        aligned = isinstance(self, Alignment)

        new_length = 0
        for seq_name in self.names:
            old_seq = self.named_seqs[seq_name]
            if not aligned:
                new_seq = old_seq.trim_stop_codon(gc=gc,
                                                  allow_partial=allow_partial)
                new_seqs.append((seq_name, new_seq))
                continue

            new_seq = old_seq.data.trim_stop_codon(gc=gc,
                                                   allow_partial=allow_partial)

            diff = len(old_seq.data._seq) - len(new_seq._seq)
            if diff and not old_seq.map.spans[-1].lost:
                new_length = max(new_length, (len(old_seq) - diff))

            # calc lengths of gaps up to last span, add to raw seq length
            seq_length = sum([len(s)
                              for s in old_seq.map.spans[:-1] if s.lost])
            seq_length += len(new_seq._seq)

            new_length = max(new_length, seq_length)

            new_seqs.append([seq_name, diff, old_seq.map, new_seq])

        if aligned:
            if new_length == 0:  # all seqs ended in a gap
                new_length = len(self)

            self_length = len(self)
            aln_diff = self_length - new_length
            assert aln_diff >= 0

            for i, data in enumerate(new_seqs):
                seq_name, diff, old_map, new_seq = data
                if not aln_diff == diff == 0:
                    # duplicate the spans
                    spans = [deepcopy(s) for s in old_map.spans]

                if diff == 0:  # seq unchanged
                    if aln_diff == 0:  # aln length unchanged
                        new_map = old_map
                    else:  # aln length changed
                        assert spans[-1].lost
                        # shrink/remove this last gap
                        spans[-1].length -= aln_diff
                        assert spans[-1].length >= 0
                        if spans[-1].length == 0:
                            del(spans[-1])
                        new_map = Map(spans=spans, parent_length=new_length)
                    seq = Aligned(new_map, new_seq)
                    new_seqs[i] = (seq_name, seq)
                    continue

                # seq length has changed, establish whether we need to adjust
                # the last one and/or to add a lost span
                index = -2
                if spans[-1].lost:
                    # grow terminal gap
                    spans[-1].length = spans[-1].length - aln_diff + diff
                elif aln_diff == 0:
                    spans.append(LostSpan(diff))
                else:
                    index = -1

                spans[index].end -= diff
                spans[index].length -= diff

                new_map = Map(spans=spans, parent_length=new_length)
                seq = Aligned(new_map, new_seq)
                new_seqs[i] = (seq_name, seq)

        return self.__class__(moltype=self.moltype, data=new_seqs, info=self.info, **kwargs)

    def get_seq_names(self):
        """Return a list of sequence names."""
        return self.names[:]

    def get_lengths(self, include_ambiguity=False, allow_gap=False):
        """returns {name: seq length, ...}

        Arguments:
            - include_ambiguity: if True, motifs containing ambiguous characters
              from the seq moltype are included. No expansion of those is attempted.
            - allow_gaps: if True, motifs containing a gap character are included.
        """
        counts = self.counts_per_seq(motif_length=1,
                                     include_ambiguity=include_ambiguity,
                                     allow_gap=allow_gap)
        lengths = counts.row_sum()
        return lengths

    def counts_per_seq(self, motif_length=1, include_ambiguity=False, allow_gap=False):
        """returns dict of counts of motifs per sequence

            only non-overlapping motifs are counted.

            Arguments:
            - motif_length: number of characters per tuple.
            - include_ambiguity: if True, motifs containing ambiguous characters
              from the seq moltype are included. No expansion of those is attempted.
            - allow_gaps: if True, motifs containing a gap character are included.
            """
        # this is overridden for Alignments, so just rely on the sequence counts
        # method
        counts = {}
        for seq in self.seqs:
            c = seq.counts(motif_length=motif_length,
                           include_ambiguity=include_ambiguity,
                           allow_gap=allow_gap)
            counts[seq.name] = c
        return counts

    
    def counts(self, motif_length=1, include_ambiguity=False, allow_gap=False):
        """returns dict of counts of motifs

            only non-overlapping motifs are counted.

            Arguments:
            - motif_length: number of elements per character.
            - include_ambiguity: if True, motifs containing ambiguous characters
              from the seq moltype are included. No expansion of those is attempted.
            - allow_gaps: if True, motifs containing a gap character are included.
            """
        per_seq = self.counts_per_seq(motif_length=motif_length,
                                            include_ambiguity=include_ambiguity,
                                            allow_gap=allow_gap)
        counts = Counter()
        for vals in per_seq.values():
            counts.update(vals)
        return counts
    
    
    def get_motif_probs(self, alphabet=None, include_ambiguity=False,
                        exclude_unobserved=False, allow_gap=False, pseudocount=0):
        """Return a dictionary of motif probs, calculated as the averaged
        frequency across sequences.

        Arguments:
            - include_ambiguity: if True resolved ambiguous codes are
              included in estimation of frequencies, default is False.
            - exclude_unobserved: if True, motifs that are not present in
              the alignment are excluded from the returned dictionary,
              default is False.
            - allow_gap: allow gap motif
        """
        if alphabet is None:
            alphabet = self.moltype.alphabet
            if allow_gap:
                alphabet = alphabet.gapped

        counts = {}
        for seq_name in self.names:
            sequence = self.named_seqs[seq_name]
            motif_len = alphabet.get_motif_len()
            if motif_len > 1:
                posns = list(range(0, len(sequence) +
                                   1 - motif_len, motif_len))
                sequence = [sequence[i:i + motif_len] for i in posns]
            for motif in sequence:
                if not allow_gap:
                    if self.moltype.gap in motif:
                        continue

                if motif in counts:
                    counts[motif] += 1
                else:
                    counts[motif] = 1

        probs = {}
        if not exclude_unobserved:
            for motif in alphabet:
                probs[motif] = pseudocount

        for (motif, count) in list(counts.items()):
            motif_set = alphabet.resolve_ambiguity(motif)
            if len(motif_set) > 1:
                if include_ambiguity:
                    count = float(count) / len(motif_set)
                else:
                    continue
            for motif in motif_set:
                probs[motif] = probs.get(motif, pseudocount) + count

        total = float(sum(probs.values()))
        for motif in probs:
            probs[motif] /= total

        return probs

    def _make_gaps_ok(self, allowed_gap_frac):
        """Makes the gaps_ok function used by omit_gap_pos and omit_gap_seqs.

        Need to make the function because if it's a method of Alignment, it
        has unwanted 'self' and 'allowed_gap_frac' parameters that impede the
        use of map() in take_seqs_if.

        WARNING: may not work correctly if component sequences have gaps that
        are not the Alignment gap character. This is because the gaps are
        checked at the position level (and the positions are lists), rather than
        at the sequence level. Working around this issue would probably cause a
        significant speed penalty.
        """
        def gaps_ok(seq):
            seq_len = len(seq)
            try:
                num_gaps = seq.count_gaps()
            except AttributeError:
                num_gaps = len(
                    list(filter(self.moltype.gaps.__contains__, seq)))
            return num_gaps / seq_len <= allowed_gap_frac

        return gaps_ok

    def omit_gap_seqs(self, allowed_gap_frac=0):
        """Returns new alignment with seqs that have <= allowed_gap_frac.

        allowed_gap_frac should be a fraction between 0 and 1 inclusive.
        Default is 0.
        """
        gaps_ok = self._make_gaps_ok(allowed_gap_frac)

        return self.take_seqs_if(gaps_ok)

    def omit_gap_runs(self, allowed_run=1):
        """Returns new alignment where all seqs have runs of gaps <=allowed_run.

        Note that seqs with exactly allowed_run gaps are not deleted.
        Default is for allowed_run to be 1 (i.e. no consecutive gaps allowed).

        Because the test for whether the current gap run exceeds the maximum
        allowed gap run is only triggered when there is at least one gap, even
        negative values for allowed_run will still let sequences with no gaps
        through.
        """
        def ok_gap_run(x):
            try:
                is_gap = x.alphabet.gaps.__contains__
            except AttributeError:
                is_gap = self.moltype.gaps.__contains__
            curr_run = max_run = 0
            for i in x:
                if is_gap(i):
                    curr_run += 1
                    if curr_run > allowed_run:
                        return False
                else:
                    curr_run = 0
            # can only get here if max_run was never exceeded (although this
            # does include the case where the sequence is empty)
            return True

        return self.take_seqs_if(ok_gap_run)

    def to_moltype(self, moltype):
        """returns copy of self with moltype seqs"""
        make_seq = moltype.make_seq
        data = [make_seq(s, s.name) for s in self.seqs]
        new = self.__class__(data=data, moltype=moltype, name=self.name,
                             info=self.info)
        if isinstance(self, _Annotatable) and self.annotations:
            new.annotations = self.annotations[:]
        return new

    def to_dna(self):
        """returns copy of self as an alignment of DNA moltype seqs"""
        make_seq = cogent3.DNA.make_seq
        data = [make_seq(s, s.name) for s in self.seqs]
        new = self.__class__(data=data, moltype=cogent3.DNA, name=self.name,
                             info=self.info)
        if isinstance(self, _Annotatable) and self.annotations:
            new.annotations = self.annotations[:]
        return new

    def to_rna(self):
        """returns copy of self as an alignment of RNA moltype seqs"""
        make_seq = cogent3.RNA.make_seq
        data = [make_seq(s, s.name) for s in self.seqs]
        new = self.__class__(data=data, moltype=cogent3.RNA, name=self.name,
                             info=self.info)
        if isinstance(self, _Annotatable) and self.annotations:
            new.annotations = self.annotations[:]
        return new

    def to_protein(self):
        """returns copy of self as an alignment of PROTEIN moltype seqs"""
        make_seq = cogent3.PROTEIN.make_seq
        data = [make_seq(s, s.name) for s in self.seqs]
        new = self.__class__(data=data, moltype=cogent3.PROTEIN, name=self.name,
                             info=self.info)
        if isinstance(self, _Annotatable) and self.annotations:
            new.annotations = self.annotations[:]
        return new

    def rc(self):
        """Returns the reverse complement alignment"""
        seqs = [self.named_seqs[name].rc() for name in self.names]
        rc = self.__class__(data=seqs, names=self.names[
                            :], name=self.name, info=self.info)
        if isinstance(self, _Annotatable) and self.annotations:
            self._annotations_nucleic_reversed_on(rc)
        return rc

    def reversecomplement(self):
        """Returns the reverse complement alignment. A synonymn for rc."""
        return self.rc()

    def pad_seqs(self, pad_length=None, **kwargs):
        """Returns copy in which sequences are padded to same length.

            pad_length: Length all sequences are to be padded to.  Will pad
                to max sequence length if pad_length is None or less than max
                length.
        """
        # get max length
        max_len = max([len(s) for s in self.seqs])
        # If a pad_length was passed in, make sure it is valid
        if pad_length is not None:
            pad_length = int(pad_length)
            if pad_length < max_len:
                raise ValueError("pad_length must be at greater or equal to maximum sequence length: %s"
                                 % (str(max_len)))
        # pad_length is max sequence length.
        else:
            pad_length = max_len

        # Get new sequence list
        new_seqs = []
        aligned = isinstance(self, Alignment)
        # for each sequence, pad gaps to end
        for seq_name in self.names:
            if aligned:
                seq = self.named_seqs[seq_name].data
            else:
                seq = self.named_seqs[seq_name]
            padded_seq = seq + '-' * (pad_length - len(seq))
            new_seqs.append((seq_name, padded_seq))

        # return new SequenceCollection object
        return SequenceCollection(moltype=self.moltype, data=new_seqs, **kwargs)


@total_ordering
class Aligned(object):
    """One sequence in an alignment, a map between alignment coordinates and
    sequence coordinates"""

    def __init__(self, map, data, length=None):
        # Unlike the normal map constructor, here we take a list of pairs of
        # alignment coordinates, NOT a list of pairs of sequence coordinates
        if isinstance(map, list):
            map = Map(map, parent_length=length).inverse()
        self.map = map
        self.data = data
        if hasattr(data, 'info'):
            self.info = data.info
        if hasattr(data, 'name'):
            self.name = data.name

    def _get_moltype(self):
        return self.data.moltype
    moltype = property(_get_moltype)

    def copy(self, memo=None, _nil=[], constructor='ignored'):
        """Returns a shallow copy of self

        WARNING: cogent3.core.sequence.Sequence does NOT implement a copy method,
        as such, the data member variable of the copied object will maintain
        reference to the original object.

        WARNING: cogent3.core.location.Map does NOT implement a copy method, as
        such, the data member variable of the copied object will maintain
        reference to the original object.
        """
        return self.__class__(self.map, self.data)

    def __repr__(self):
        return '%s of %s' % (repr(self.map), repr(self.data))

    def with_termini_unknown(self):
        return self.__class__(self.map.with_termini_unknown(), self.data)

    def copy_annotations(self, other):
        self.data.copy_annotations(other)

    def annotate_from_gff(self, f):
        self.data.annotate_from_gff(f)

    def add_feature(self, *args, **kwargs):
        self.data.add_feature(*args, **kwargs)

    def __str__(self):
        """Returns string representation of aligned sequence, incl. gaps."""
        return str(self.get_gapped_seq())

    def __lt__(self, other):
        """Compares based on string representations."""
        return str(self) < str(other)

    def __eq__(self, other):
        """Compares based on string representations."""
        return str(self) == str(other)

    def __hash__(self):
        return hash(self.data)

    def __ne__(self, other):
        """Compares based on string representations."""
        return str(self) != str(other)

    def __iter__(self):
        """Iterates over sequence one motif (e.g. char) at a time, incl. gaps"""
        return self.data.gapped_by_map_motif_iter(self.map)

    def get_gapped_seq(self, recode_gaps=False, moltype=None):
        """Returns sequence as an object, including gaps."""
        return self.data.gapped_by_map(self.map, recode_gaps)

    def __len__(self):
        # these make it look like Aligned should be a subclass of Map,
        # but then you have to be careful with __getitem__, __init__ and
        # inverse.
        return len(self.map)

    def __add__(self, other):
        if self.data is other.data:
            (map, seq) = (self.map + other.map, self.data)
        else:
            seq = self.get_gapped_seq() + other.get_gapped_seq()
            (map, seq) = seq.parse_out_gaps()
        return Aligned(map, seq)

    def __getitem__(self, slice):
        return Aligned(self.map[slice], self.data)

    def rc(self):
        return Aligned(self.map.reversed(), self.data)

    def to_rna(self):
        return Aligned(self.map, self.data.to_rna())

    def to_dna(self):
        return Aligned(self.map, self.data.to_dna())

    def to_rich_dict(self):
        data = self.data.to_rich_dict()
        data['seq'] = str(self)
        return data

    def to_json(self):
        """returns a json formatted string"""
        return json.dumps(self.to_rich_dict())

    def get_tracks(self, policy):
        policy = policy.at(self.map.inverse())
        return self.data.get_tracks(policy)

    def remapped_to(self, map):
        result = Aligned(map[self.map.inverse()].inverse(), self.data)
        return result

    def get_annotations_matching(self, alignment, *args):
        for annot in self.data.get_annotations_matching(*args):
            yield annot.remapped_to(alignment, self.map.inverse())

    def gap_vector(self):
        """Returns gap_vector of GappedSeq, for omit_gap_pos."""
        return self.get_gapped_seq().gap_vector()

    def _masked_annotations(self, annot_types, mask_char, shadow):
        """returns a new aligned sequence with regions defined by align_spans
        and shadow masked."""
        new_data = self.data.with_masked_annotations(
            annot_types, mask_char, shadow)
        # we remove the mask annotations from self and new_data
        return self.__class__(self.map, new_data)


class AlignmentI(object):
    """Alignment interface object. Contains methods shared by implementations.

    Note that subclasses should inherit both from AlignmentI and from
    SequenceCollection (typically).

    Alignments are expected to be immutable once created. No mechanism is
    provided for maintaining reference consistency if data in the alignment
    are modified.

    An Alignment is expected to be able to generate the following:
    - seqs:         Sequence objects in the alignment, can turn themselves into
                    strings. These are usually thought of as "rows" in an
                    alignment.
    - positions:    Vectors representing data in each position in the alignment
                    These are usually thought of as "columns" in an alignment.
    - seq_data:      Vectors representing data in each sequence in the alignment,
                    not necessarily guaranteed to turn themselves into a string
    - names:        List of names of sequences in the alignment. Used for
                    display order. A cheap way to omit or reorder sequences is
                    to modify the list of names.
    - named_seqs:    Dict of name -> seq object, used for lookup.
    - moltype:      moltype of the alignment.
    """
    default_gap = '-'  # default gap character for padding
    gap_chars = dict.fromkeys('-?')  # default gap chars for comparisons

    def iter_positions(self, pos_order=None):
        """Iterates over positions in the alignment, in order.

        pos_order refers to a list of indices (ints) specifying the column
        order. This lets you rearrange positions if you want to (e.g. to pull
        out individual codon positions).

        Note that self.iter_positions() always returns new objects, by default
        lists of elements. Use map(f, self.iter_positions) to apply the
        constructor or function f to the resulting lists (f must take a single
        list as a parameter).

        Will raise IndexError if one of the indices in order exceeds the
        sequence length. This will always happen on ragged alignments:
        assign to self.seq_len to set all sequences to the same length.
        """
        get = self.named_seqs.__getitem__
        pos_order = pos_order or range(self.seq_len)
        seq_order = self.names
        for pos in pos_order:
            yield [get(seq)[pos] for seq in seq_order]

    positions = property(iter_positions)

    def take_positions(self, cols, negate=False):
        """Returns new Alignment containing only specified positions.

        By default, the seqs will be lists, but an alternative constructor
        can be specified.

        Note that take_positions will fail on ragged positions.
        """
        make_seq = self.moltype.make_seq
        result = {}
        # if we're negating, pick out all the positions except specified
        # indices
        if negate:
            col_lookup = dict.fromkeys(cols)
            for name, seq in list(self.named_seqs.items()):
                result[name] = make_seq([seq[i] for i in range(len(seq))
                                         if i not in col_lookup])
        # otherwise, just get the requested indices
        else:
            for name, seq in list(self.named_seqs.items()):
                result[name] = make_seq([seq[i] for i in cols])
        return self.__class__(result, names=self.names, info=self.info)

    def get_position_indices(self, f, native=False, negate=False):
        """Returns list of column indices for which f(col) is True.

        f : callable
          function that returns true/false given an alignment position
        native : boolean
          if True, and ArrayAlignment, f is provided with slice of array
          otherwise the string is used
        negate : boolean
          if True, not f() is used
        """
        # negate f if necessary
        if negate:
            def new_f(x): return not f(x)
        else:
            new_f = f

        if native and isinstance(self, ArrayAlignment):
            result = [i for i in range(
                self.seq_len) if new_f(self.array_seqs[:, i])]
        else:
            result = [i for i, col in enumerate(
                self.positions) if new_f(col)]

        return result

    def take_positions_if(self, f, negate=False):
        """Returns new Alignment containing cols where f(col) is True.
        """
        return self.take_positions(self.get_position_indices(f, negate=negate))

    def iupac_consensus(self, alphabet=None):
        """Returns string containing IUPAC consensus sequence of the alignment.
        """
        if alphabet is None:
            alphabet = self.moltype
        consensus = []
        degen = alphabet.degenerate_from_seq
        for col in self.positions:
            consensus.append(degen(coerce_to_string(col)))
        return coerce_to_string(consensus)

    def column_freqs(self, constructor=Freqs):
        """Returns list of Freqs with item counts for each column.
        """
        return list(map(constructor, self.positions))

    def column_probs(self, constructor=Freqs):
        """Returns FrequencyDistribuutions w/ prob. of each item per column.

        Implemented as a list of normalized Freqs objects.
        """
        freqs = self.column_freqs(constructor)

        for fd in freqs:
            fd.normalize()
        return freqs

    def majority_consensus(self, transform=None, constructor=Freqs):
        """Returns list containing most frequent item at each position.

        Optional parameter transform gives constructor for type to which result
        will be converted (useful when consensus should be same type as
        originals).
        """
        col_freqs = self.column_freqs(constructor)

        consensus = [freq.Mode for freq in col_freqs]
        if transform == str:
            return coerce_to_string(consensus)
        elif transform:
            return transform(consensus)
        else:
            return consensus

    def uncertainties(self, good_items=None):
        """Returns Shannon uncertainty at each position.

        Usage: information_list = alignment.information(good_items=None)

        If good_items is supplied, deletes any symbols that are not in
        good_items.
        """
        uncertainties = []
        # calculate column probabilities if necessary
        if hasattr(self, 'PositionumnProbs'):
            probs = self.PositionumnProbs
        else:
            probs = self.column_probs()
        # calculate uncertainty for each column
        for prob in probs:
            # if there's a list of valid symbols, need to delete everything
            # else
            if good_items:
                prob = prob.copy()  # do not change original
                # get rid of any symbols not in good_items
                for symbol in list(prob.keys()):
                    if symbol not in good_items:
                        del prob[symbol]
                # normalize the probabilities and add to the list
                prob.normalize()
            uncertainties.append(prob.Uncertainty)
        return uncertainties

    def get_pwm(self):
        """Returns a position specific weight matrix for the alignment."""
        result = []
        for col in self.positions:
            f = Freqs(col)
            result.append([f.get(c, 0) for c in self.moltype])
        darr = DictArrayTemplate(len(self), list(self.moltype))
        result = darr.wrap(result)
        return result

    def _get_freqs(self, index=None):
        """Gets array of freqs along index 0 (= positions) or 1 (= seqs).

        index: if 0, will calculate the frequency of each symbol in each
        position (=column) in the alignment. Will return 2D array where the
        first index is the position, and the second index is the index of the
        symbol in the alphabet. For example, for the TCAG DNA Alphabet,
        result[3][0] would store the count of T at position 3 (i.e. the 4th
        position in the alignment.

        if 1, does the same thing except that the calculation is performed for
        each sequence, so the 2D array has the sequence index as the first
        index, and the symbol index as the second index. For example, for the
        TCAG DNA Alphabet, result[3][0] would store the count of T in the
        sequence at index 3 (i.e. the 4th sequence).

        First an DenseAligment object is created, next the calculation is done
        on this object. It is important that the ArrayAlignment is initialized
        with the same moltype and Alphabet as the original Alignment.
        """
        da = ArrayAlignment(self, moltype=self.moltype, alphabet=self.alphabet)
        return da._get_freqs(index)

    def get_seq_freqs(self):
        """Returns Profile of counts: seq by character.

        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-sequence (i.e. per row).
        """
        return Profile(self._get_freqs(0), self.alphabet)

    def get_pos_freqs(self):
        """Returns Profile of counts: position by character.

        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-position (i.e. per column).
        """
        return Profile(self._get_freqs(1), self.alphabet)

    def no_degenerates(self, motif_length=1, allow_gap=False):
        """returns new alignment without degenerate characters

        Arguments:
        - motif_length: sequences are segmented into units of this size
        - allow_gaps: whether gaps are to be treated as a degenerate
          character (default, most evolutionary modelling treats gaps as
          N) or not.
        """
        try:
            chars = list(self.moltype.alphabet.non_degen)
        except AttributeError:
            msg = "Invalid MolType (no degenerate characters), "\
                "create the alignment using DNA, RNA or PROTEIN"
            raise ValueError(msg)

        is_array = isinstance(self, ArrayAlignment)
        alpha = self.alphabet
        if allow_gap:
            chars.extend(self.moltype.gap)
            alpha = self.moltype.alphabet.gapped

        if is_array:
            chars = list(map(alpha.index, chars))

        predicate = AllowedCharacters(chars, is_array=is_array)
        new = self.filtered(predicate, motif_length=motif_length)
        return new

    def omit_gap_pos(self, allowed_gap_frac=1 - eps, motif_length=1):
        """Returns new alignment where all cols (motifs) have <= allowed_gap_frac gaps.

        Argments:
        - allowed_gap_frac: specifies proportion of gaps is allowed in each
          column (default is just < 1, i.e. only cols with at least one gap
          character are preserved). Set to 1 - e-6 to exclude strictly gapped
          columns.
        - motif_length: set's the "column" width, e.g. setting to 3
          corresponds to codons. A motif that includes a gap at any position
          is included in the counting. Default is 1.
        """
        is_array = isinstance(self, ArrayAlignment)
        try:
            alpha = self.moltype.alphabets.degen_gapped
        except:
            alpha = self.moltype.alphabet

        gaps = list(self.moltype.gaps)
        if is_array:
            gaps = list(map(alpha.index, gaps))

        gaps_ok = GapsOk(gaps, allowed_gap_frac, is_array=is_array)
        # if we're not deleting the 'naughty' seqs that contribute to the
        # gaps, it's easy...
        result = self.filtered(gaps_ok, motif_length=motif_length)
        return result

    def omit_bad_seqs(self, disallowed_frac=0.9, allowed_frac_bad_cols=0, exclude_just_gap=True):
        """Returns new alignment without the sequences responsible for exceeding disallowed_frac.

        disallowed_frac: used to identify columns with excessive gaps
        allowed_frac_bad_cols: the fraction of gaps induced by the sequence
        exclude_just_gap: sequences that are just gaps are excluded
        """
        is_array = isinstance(self, ArrayAlignment)
        alpha = self.alphabet

        gaps = list(self.moltype.gaps)
        if is_array:
            gaps = list(map(alpha.index, gaps))

        excess_gaps = GapsOk(gaps, disallowed_frac,
                             is_array=is_array, negate=True)
        gapped_cols = self.get_position_indices(excess_gaps, native=True)

        seqs_to_delete = {}
        if exclude_just_gap:
            # any seqs that are all gaps?
            length = len(self)
            for name in self.names:
                num_gaps = self.count_gaps(name)
                if num_gaps == length:
                    seqs_to_delete[name] = True

        # we treat array and string data the same, likely inefficient
        is_gap = list(self.moltype.gaps).__contains__
        bad_cols_per_seq = Counter()
        for name, seq in list(self.named_seqs.items()):
            if name in seqs_to_delete:
                continue

            seq = str(seq)
            for col in gapped_cols:
                element = seq[col]
                if not is_gap(element):
                    bad_cols_per_seq[name] += 1

        # figure out which of the seqs cause too many gaps
        get = self.named_seqs.__getitem__
        for name, count in list(bad_cols_per_seq.items()):
            if count / len(get(name)) >= allowed_frac_bad_cols:
                seqs_to_delete[name] = True

        if not seqs_to_delete:
            return self

        good_seqs = self.take_seqs(seqs_to_delete, negate=True)
        return good_seqs

    def matching_ref(self, ref_name, gap_fraction, gap_run):
        """Returns new alignment with seqs well aligned with a reference.

        gap_fraction = fraction of positions that either have a gap in the
            template but not in the seq or in the seq but not in the template
        gap_run = number of consecutive gaps tolerated in query relative to
            sequence or sequence relative to query
        """
        template = self.named_seqs[ref_name]
        gap_filter = make_gap_filter(template, gap_fraction, gap_run)
        return self.take_seqs_if(gap_filter)

    def sample(self, n=None, with_replacement=False, motif_length=1,
               randint=randint, permutation=permutation):
        """Returns random sample of positions from self, e.g. to bootstrap.

        Arguments:
            - n: the number of positions to sample from the alignment.
              Default is alignment length
            - with_replacement: boolean flag for determining if sampled
              positions
            - random_series: a random number generator with
              .randint(min,max) .random() methods


        Notes:
            By default (resampling all positions without replacement), generates
            a permutation of the positions of the alignment.

            Setting with_replacement to True and otherwise leaving parameters
            as defaults generates a standard bootstrap resampling of the
            alignment.
            """
        population_size = len(self) // motif_length
        if not n:
            n = population_size
        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            assert n <= population_size, (n, population_size, motif_length)
            locations = permutation(population_size)[:n]
        positions = [(loc * motif_length, (loc + 1) * motif_length)
                     for loc in locations]
        sample = Map(positions, parent_length=len(self))
        return self.gapped_by_map(sample, info=self.info)

    def sliding_windows(self, window, step, start=None, end=None):
        """Generator yielding new Alignments of given length and interval.

        Arguments:
            - window: The length of each returned alignment.
            - step: The interval between the start of the successive
              alignment objects returned.
            - start: first window start position
            - end: last window start position
        """
        start = [start, 0][start is None]
        end = [end, len(self) - window + 1][end is None]
        end = min(len(self) - window + 1, end)
        if start < end and len(self) - end >= window - 1:
            for pos in range(start, end, step):
                yield self[pos:pos + window]

    def _get_raw_pretty(self, name_order):
        """returns dict {name: seq, ...} for pretty print"""
        if name_order is not None:
            assert set(name_order) == set(self.names), "names don't match"

        names = name_order or self.names
        output = defaultdict(list)
        names = name_order or self.names
        num_seqs = len(names)

        seqs = []
        for name in names:
            seq = str(self.named_seqs[name])
            seqs.append(seq)

        positions = list(zip(*seqs))

        for position in positions:
            ref = position[0]
            output[names[0]].append(ref)
            for seq_num in range(1, num_seqs):
                val = '.' if position[seq_num] == ref else position[seq_num]
                output[names[seq_num]].append(val)

        return names, output

    def _repr_html_(self):
        # we put the longest sequence first
        html = self.to_html(longest_ref=True, limit=2000)
        return html

    def to_html(self, name_order=None, interleave_len=60, limit=None,
                longest_ref=True, colors=None,
                font_size=12, font_family='Lucida Console'):
        '''returns html with embedded styles for sequence colouring

        Arguments:
            - name_order: order of names for display.
            - interleave_len: maximum number of printed bases, defaults to
              alignment length
            - limit: truncate alignment to this length
            - longest_ref: If True, the longest sequence (excluding gaps and
              ambiguities) is selected as the reference.
            - colors: {character: color, ...}. Defaults to those defined by
              moltype.
            - font_size: in points. Affects labels and sequence and line spacing
              (proportional to value)
            - font_family: string denoting font family

        To display in jupyter notebook:

            >>> from IPython.core.display import HTML
            >>> HTML(aln.to_html())
        '''
        css, styles = self.moltype.get_css_style(colors=colors,
                                                 font_size=font_size,
                                                 font_family=font_family)

        if longest_ref and name_order is None:
            length_names = []
            for s in self.seqs:
                try:
                    l = len(s) - s.count_degenerate() - s.count_gaps()
                except AttributeError:
                    # We have the Aligned class, and Aligned.data is ungapped
                    nd = s.data.count_degenerate()
                    l = len(s.data) - nd
                length_names.append((l, s.name))
            length_names.sort(reverse=True)
            ref = length_names[0][1]
            name_order = list(self.names)
            name_order.remove(ref)
            name_order.insert(0, ref)

        if limit is None:
            names, output = self._get_raw_pretty(name_order)
        else:
            names, output = self[:limit]._get_raw_pretty(name_order)

        gaps = ''.join(self.moltype.gaps)
        refname = names[0]
        refseq = output[refname]
        seqlen = len(refseq)
        start_gap = re.search('^[%s]+' % gaps, ''.join(refseq))
        end_gap = re.search('[%s]+$' % gaps, ''.join(refseq))
        ref_colours = []
        start = 0 if start_gap is None else start_gap.end()
        end = len(refseq) if end_gap is None else end_gap.start()
        seq_style = []
        template = '<span class="%s">%%s</span>'
        styled_seqs = defaultdict(list)
        for i in range(seqlen):
            char = refseq[i]
            if i < start or i >= end:
                style = 'terminal_ambig_%s' % self.moltype.label
            else:
                style = styles[char]

            seq_style.append(template % style)
            styled_seqs[refname].append(seq_style[-1] % char)

        for name in names:
            if name == refname:
                continue

            seq = []
            for i, c in enumerate(output[name]):
                if c == '.':
                    s = seq_style[i] % c
                else:
                    s = template % (styles[c])
                    s = s % c
                seq.append(s)

            styled_seqs[name] = seq

        # make a html table
        seqs = array([styled_seqs[n] for n in names], dtype='O')
        table = ['<table>']
        seq_ = '<td>%s</td>'
        label_ = '<td class="label">%s</td>'
        for i in range(0, seqlen, interleave_len):
            seqblock = seqs[:, i: i + interleave_len].tolist()
            for n, s in zip(names, seqblock):
                s = ''.join(s)
                row = ''.join([label_ % n, seq_ % s])
                table.append('<tr>%s</tr>' % row)
            table.append('<tr style="background-color:gray" '
                         'class="blank_row"><td></td><td></td></tr>')
        table.append('</table>')
        if limit and limit < len(self):
            summary = ('%s x %s (truncated to %s) %s '
                       'alignment') % (len(self.names), len(self), limit,
                                       self.moltype.label)
        else:
            summary = ('%s x %s %s '
                       'alignment') % (len(self.names), len(self),
                                       self.moltype.label)

        text = ['<style>',
                'tr { line-height: %dpt ; }' % int(font_size / 4),
                '.blank_row{ line-height: %dpt !important; '
                'opacity: 0.10; }' % font_size,
                'td { border: none !important; text-align: left !important; }',
                '.label { font-size: %dpt ; text-align: right !important; '
                'color: black !important; }' % font_size,
                css,
                '</style>',
                '<body>',
                '\n'.join(table),
                '<p><i>%s</i></p>' % summary,
                '</body>']
        return '\n'.join(text)

    def to_pretty(self, name_order=None, interleave_len=None):
        """returns a string representation of the alignment in pretty print format

        Arguments:
            - name_order: order of names for display.
            - interleave_len: maximum number of printed bases, defaults to alignment length"""
        names, output = self._get_raw_pretty(name_order=name_order)
        label_width = max(list(map(len, names)))
        name_template = '{:>%d}' % label_width
        display_names = dict([(n, name_template.format(n)) for n in names])

        def make_line(label, seq): return "%s    %s" % (label, seq)
        if interleave_len is None:
            result = [make_line(display_names[n], ''.join(output[n]))
                      for n in names]
            return '\n'.join(result)

        align_length = len(self)
        result = []
        for start in range(0, align_length, interleave_len):
            for n in names:
                result.append(make_line(display_names[n], ''.join(
                    output[n][start: start + interleave_len])))

            result.append('')

        if not result[-1]:
            del(result[-1])

        return '\n'.join(result)
    
    def counts_per_seq(self, motif_length=1, include_ambiguity=False, allow_gap=False):
        """returns dict of counts of motifs per sequence
    
            only non-overlapping motifs are counted.
    
            Arguments:
            - motif_length: number of elements per character.
            - include_ambiguity: if True, motifs containing ambiguous characters
              from the seq moltype are included. No expansion of those is attempted.
            - allow_gaps: if True, motifs containing a gap character are included.
            """
        is_array = isinstance(self, ArrayAlignment)
        counts = {}
        for name in self.names:
            if is_array:
                seq = self.named_seqs[name]
            else:
                seq = self.get_gapped_seq(name)
            c = seq.counts(motif_length=motif_length,
                           include_ambiguity=include_ambiguity,
                           allow_gap=allow_gap)
            counts[name] = c
        return counts

    def count_gaps(self, seq_name):
        """returns number of gaps for seq_name"""
        is_array = isinstance(self, ArrayAlignment)
        if not is_array:
            seq = self.named_seqs[seq_name]
            return len(seq.map.gaps())

        seq_index = self.names.index(seq_name)
        seq = self.array_seqs[seq_index]
        alpha = self.alphabet
        gaps = list(map(alpha.index, self.moltype.gaps))
        counts = Counter(seq)
        num_gaps = sum(counts[g] for g in gaps)
        return num_gaps

    def variable_positions(self, include_gap_motif=True):
        """Return a list of variable position indexes.

        Arguments:
            - include_gap_motif: if False, sequences with a gap motif in a
              column are ignored."""
        seqs = [str(self.named_seqs[n]) for n in self.names]
        seq1 = seqs[0]
        positions = list(zip(*seqs[1:]))
        result = []
        for (position, (motif1, column)) in enumerate(zip(seq1, positions)):
            for motif in column:
                if motif != motif1:
                    if include_gap_motif:
                        result.append(position)
                        break
                    elif motif != '-' and motif1 != '-':
                        result.append(position)
                        break

        return result

    def to_type(self, array_align=False, moltype=None, alphabet=None):
        """returns alignment of type indicated by array_align

        Parameters
        ----------
        array_align: bool
          if True, returns as ArrayAlignment. Otherwise as "standard"
          Alignment class. Conversion to ArrayAlignment loses annotations.

        moltype : MolType instance
          overrides self.moltype

        alphabet : Alphabet instance
          overrides self.alphabet

        If array_align would result in no change (class is same as self),
        returns self
        """
        klass = ArrayAlignment if array_align else Alignment
        if isinstance(self, klass) and (moltype is None or
                                        moltype == self.moltype):
            return self

        data = self.todict()
        if moltype is None:
            # Alignment and ArrayAlignment have different default moltypes
            moltype_default = self.moltype == self.__class__.moltype
            if moltype_default:
                moltype = ArrayAlignment.moltype if array_align else Alignment.moltype
            else:
                moltype = self.moltype
        new = klass(data=data, moltype=moltype,
                    info=self.info, names=self.names)
        return new


def _one_length(seqs):
    """raises ValueError if seqs not all same length"""
    seq_lengths = set(len(s) for s in seqs)
    if len(seq_lengths) != 1:
        raise ValueError("not all sequences have same length")


def aln_from_array(a, array_type=None, alphabet=None):
    """Alignment from array of pos x seq: no change, names are integers.

    This is an InputHandler for Alignment. It converts an arbitrary array
    of numbers without change, but adds successive integer names (0-based) to
    each sequence (i.e. column) in the input a. Data type of input is
    unchanged.
    """
    if array_type is None:
        result = a.copy()
    else:
        result = a.astype(array_type)
    return transpose(result), None


def aln_from_array_seqs(seqs, array_type=None, alphabet=None):
    """Alignment from ArraySequence objects: seqs -> array, names from seqs.

    This is an InputHandler for Alignment. It converts a list of Sequence
    objects with _data and label properties into the character array Alignment
    needs. All sequences must be the same length.

    WARNING: Assumes that the ArraySeqs are already in the right alphabet. If
    this is not the case, e.g. if you are putting sequences on a degenerate
    alphabet into a non-degenerate alignment or you are putting protein
    sequences into a DNA alignment, there will be problems with the alphabet
    mapping (i.e. the resulting sequences may be meaningless).

    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    data, names = [], []
    for s in seqs:
        data.append(s._data)
        names.append(s.name)

    _one_length(data)

    result = array(data)
    if array_type:
        result = result.astype(array_type)
    return result, names


def aln_from_generic(data, array_type=None, alphabet=None):
    """Alignment from generic seq x pos data: sequence of sequences of chars.

    This is an InputHandler for Alignment. It converts a generic list (each
    item in the list will be mapped onto an Array object, with character
    transformations, all items must be the same length) into a numpy array,
    and assigns sequential integers (0-based) as names.

    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    result = array(list(map(alphabet.to_indices, data)))
    names = []
    for d in data:
        if hasattr(d, 'name'):
            names.append(d.name)
        else:
            names.append(None)
    if array_type:
        result = result.astype(array_type)
    return result, names


def aln_from_collection(seqs, array_type=None, alphabet=None):
    """Alignment from SequenceCollection object, or its subclasses."""
    names = seqs.names
    data = [seqs.named_seqs[i] for i in names]
    result = array(list(map(alphabet.to_indices, data)))
    if array_type:
        result = result.astype(array_type)
    return result, names


def aln_from_fasta(seqs, array_type=None, alphabet=None):
    """Alignment from FASTA-format string or lines.

    This is an InputHandler for Alignment. It converts a FASTA-format string
    or collection of lines into an Alignment object. All sequences must be the
    same length.

    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    if isinstance(seqs, bytes):
        seqs = seqs.decode("utf-8")
    if isinstance(seqs, str):
        seqs = seqs.splitlines()
    return aln_from_array_seqs([ArraySequence(s, name=l, alphabet=alphabet)
                                for l, s in cogent3.parse.fasta.MinimalFastaParser(seqs)], array_type)


def aln_from_dict(aln, array_type=None, alphabet=None):
    """Alignment from dict of {label:seq_as_str}.

    This is an InputHandler for Alignment. It converts a dict in which the
    keys are the names and the values are the sequences (sequence only, no
    whitespace or other formatting) into an alignment. Because the dict
    doesn't preserve order, the result will be in alphabetical order."""
    for n in aln:
        try:
            aln[n] = aln[n].upper()
        except AttributeError:
            pass

    names, seqs = list(zip(*sorted(aln.items())))
    seqs = [bytes_to_string(s) for s in seqs]
    _one_length(seqs)
    result = array(list(map(alphabet.to_indices, seqs)), array_type)
    return result, list(names)


def aln_from_kv_pairs(aln, array_type=None, alphabet=None):
    """Alignment from sequence of (key, value) pairs.

    This is an InputHandler for Alignment. It converts a list in which the
    first item of each pair is the label and the second item is the sequence
    (sequence only, no whitespace or other formatting) into an alignment.
    Because the dict doesn't preserve order, the result will be in arbitrary
    order."""
    names, seqs = list(zip(*aln))
    seqs = [bytes_to_string(s) for s in seqs]
    _one_length(seqs)
    result = array(list(map(alphabet.to_indices, seqs)), array_type)
    return result, list(names)


def aln_from_array_aln(aln, array_type=None, alphabet=None):
    """Alignment from existing ArrayAlignment object: copies data.

    Retrieves data from positions field. Uses copy(), so array data type
    should be unchanged.
    """
    if array_type is None:
        result = aln.array_positions.copy()
    else:
        result = aln.array_positions.astype(array_type)
    return transpose(result), aln.names[:]


def aln_from_empty(obj, *args, **kwargs):
    """Alignment from empty data: raise exception."""
    raise ValueError("Cannot create empty alignment.")

# Implementation of Alignment base class


class ArrayAlignment(AlignmentI, SequenceCollection):
    """Holds a dense array representing a multiple sequence alignment.

    An Alignment is _often_, but not necessarily, an array of chars. You might
    want to use some other data type for the alignment if you have a large
    number of symbols. For example, codons on an ungapped DNA alphabet has
    4*4*4=64 entries so can fit in a standard char data type, but tripeptides
    on the 20-letter ungapped protein alphabet has 20*20*20=8000 entries so
    can _not_ fit in a char and values will wrap around (i.e. you will get an
    unpredictable, wrong value for any item whose index is greater than the
    max value, e.g. 255 for uint8), so in this case you would need to use
    UInt16, which can hold 65536 values. DO NOT USE SIGNED DATA TYPES FOR YOUR
    ALIGNMENT ARRAY UNLESS YOU LOVE MISERY AND HARD-TO-DEBUG PROBLEMS.

    Implementation: aln[i] returns position i in the alignment.

    aln.positions[i] returns the same as aln[i] -- usually, users think of this
    as a 'column', because alignment editors such as Clustal typically display
    each sequence as a row so a position that cuts across sequences is a
    column.

    aln.seqs[i] returns a sequence, or 'row' of the alignment in standard
    terminology.

    WARNING: aln.seqs and aln.positions are different views of the same array,
    so if you change one you will change the other. This will no longer be
    true if you assign to seqs or positions directly, so don't do it. If you
    want to change the data in the whole array, always assign to a slice so
    that both views update: aln.seqs[:] = x instead of aln.seqs = x. If you
    get the two views out of sync, you will get all sorts of exceptions. No
    validation is performed on aln.seqs and aln.positions for performance
    reasons, so this can really get you into trouble.

    Alignments are immutable, though this is not enforced. If you change the
    data after the alignment is created, all sorts of bad things might happen.

    Class properties:
    alphabet: should be an Alphabet object. Must provide mapping between items
    (possibly, but not necessarily, characters) in the alignment and indices
    of those characters in the resulting Alignment object.

    SequenceType: Constructor to use when building sequences. Default: Sequence

    InputHandlers: dict of {input_type:input_handler} where input_handler is
    from the InputHandlers above and input_type is a result of the method
    self._guess_input_type (should always be a string).

    Creating a new array will always result in a new object unless you use
    the force_same_object=True parameter.

    WARNING: Rebinding the names attribute in a ArrayAlignment is not
    recommended because not all methods will use the updated name order. This
    is because the original sequence and name order are used to produce data
    structures that are cached for efficiency, and are not updated if you
    change the names attribute.

    WARNING: ArrayAlignment strips off info objects from sequences that have
    them, primarily for efficiency.
    """
    moltype = None  # will be set to BYTES on moltype import
    alphabet = None  # will be set to BYTES.alphabet on moltype import

    InputHandlers = {'array': aln_from_array,
                     'array_seqs': aln_from_array_seqs,
                     'generic': aln_from_generic,
                     'fasta': aln_from_fasta,
                     'array_aln': aln_from_array_aln,
                     'aln': aln_from_collection,
                     'collection': aln_from_collection,
                     'dict': aln_from_dict,
                     'kv_pairs': aln_from_kv_pairs,
                     'empty': aln_from_empty,
                     }

    def __init__(self, *args, **kwargs):
        """Returns new ArrayAlignment object. Inherits from SequenceCollection.
        """
        kwargs['suppress_named_seqs'] = True
        super(ArrayAlignment, self).__init__(*args, **kwargs)
        self.array_positions = transpose(
            self.seq_data.astype(self.alphabet.array_type))
        self.array_seqs = transpose(self.array_positions)
        self.seq_data = self.array_seqs
        self.seq_len = len(self.array_positions)
        self._type = self.moltype.gettype()

    def _force_same_data(self, data, names):
        """Forces array that was passed in to be used as selfarray_positions"""
        if isinstance(data, ArrayAlignment):
            data = data._positions
        self.array_positions = data
        self.names = names or self.DefaultNameFunction(len(data[0]))

    def _get_positions(self):
        """Override superclass positions to return positions as symbols."""
        return list(map(self.alphabet.from_indices, self.array_positions))

    positions = property(_get_positions)

    def _get_named_seqs(self):
        if not hasattr(self, '_named_seqs'):
            seqs = list(map(self.alphabet.to_string, self.array_seqs))
            if self.moltype:
                seqs = [self.moltype.make_seq(
                    s, preserve_case=True) for s in seqs]
            self._named_seqs = self._make_named_seqs(self.names, seqs)
        return self._named_seqs

    named_seqs = property(_get_named_seqs)

    def keys(self):
        """Supports dict-like interface: returns names as keys."""
        return self.names

    def values(self):
        """Supports dict-like interface: returns seqs as Sequence objects."""
        return [self.alphabet.moltype.make_array_seq(i, alphabet=self.alphabet)
                for i in self.array_seqs]

    def items(self):
        """Supports dict-like interface; returns (name, seq) pairs."""
        return list(zip(list(self.keys()), list(self.values())))

    def __iter__(self):
        """iter(aln) iterates over positions, returning array slices.

        Each item in the result is be a position ('column' in standard
        terminology) within the alignment, with the sequneces in the same
        order as in the names.

        The result shares data with the original array, so if you change
        the result you change the Alignment.
        """
        return iter(self.positions)

    def __getitem__(self, item):
        if not isinstance(item, slice):
            data = self.array_seqs[:, item]
            data = vstack(data)
        else:
            data = self.array_seqs[:, item]
        return self.__class__(data.T, list(map(str, self.names)), self.alphabet,
                              conversion_f=aln_from_array, info=self.info)

    def _coerce_seqs(self, seqs, is_array):
        """Controls how seqs are coerced in _names_seqs_order.

        Override in subclasses where this behavior should differ.
        """
        return seqs

    def get_sub_alignment(self, seqs=None, pos=None, invert_seqs=False,
                          invert_pos=False):
        """Returns subalignment of specified sequences and positions.

        seqs and pos can be passed in as lists of sequence indices to keep
        or positions to keep.

        invert_seqs: if True (default False), gets everything _except_ the
        specified sequences.

        invert_pos: if True (default False), gets everything _except_ the
        specified positions.

        Unlike most of the other code that gets things out of an alignment,
        this method returns a new alignment that does NOT share data with the
        original alignment.
        """
        # figure out which positions to keep, and keep them
        if pos is not None:
            if invert_pos:
                pos_mask = ones(len(self.array_positions))
                put(pos_mask, pos, 0)
                pos = nonzero(pos_mask)[0]
            data = take(self.array_positions, pos, axis=0)
        else:
            data = self.array_positions
        # figure out which sequences to keep, and keep them
        if seqs is not None:
            if invert_seqs:
                seq_mask = ones(len(self.array_seqs))
                put(seq_mask, seqs, 0)
                seqs = nonzero(seq_mask)[0]
            data = take(data, seqs, 1)
            names = [self.names[i] for i in seqs]
        else:
            names = self.names
        return self.__class__(data, list(map(str, names)), self.alphabet,
                              conversion_f=aln_from_array, info=self.info)

    def __str__(self):
        """Returns FASTA-format string.

        Should be able to handle joint alphabets, e.g. codons.
        """
        result = []
        names = list(map(str, self.names))
        max_label_length = max(list(map(len, names))) + 1
        seq2str = self.alphabet.from_indices
        for l, s in zip(self.names, self.array_seqs):
            result.append('>' + str(l) + '\n' + ''.join(seq2str(s)))
        return '\n'.join(result) + '\n'

    def __repr__(self):
        seqs = []
        limit = 10
        delimiter = ''
        for (count, name) in enumerate(self.names):
            if count == 3:
                seqs.append('...')
                break
            elts = list(str(self.named_seqs[name])[:limit + 1])
            if len(elts) > limit:
                elts.append('...')
            seqs.append("%s[%s]" % (name, delimiter.join(elts)))
        seqs = ', '.join(seqs)

        return "%s x %s %s alignment: %s" % (len(self.names),
                                             self.seq_len, self._type, seqs)
    def _get_freqs(self, index=None):
        """Gets array of freqs along index 0 (= positions) or 1 (= seqs).

        index: if 0, will calculate the frequency of each symbol in each
        position (=column) in the alignment. Will return 2D array where the
        first index is the position, and the second index is the index of the
        symbol in the alphabet. For example, for the TCAG DNA Alphabet,
        result[3][0] would store the count of T at position 3 (i.e. the 4th
        position in the alignment.

        if 1, does the same thing except that the calculation is performed for
        each sequence, so the 2D array has the sequence index as the first
        index, and the symbol index as the second index. For example, for the
        TCAG DNA Alphabet, result[3][0] would store the count of T in the
        sequence at index 3 (i.e. the 4th sequence).
        """
        if index:
            a = self.array_positions
        else:
            a = self.array_seqs
        count_f = self.alphabet.counts
        return array(list(map(count_f, a)))

    def get_pos_freqs(self):
        """Returns Profile of counts: position by character.

        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-position (i.e. per column).
        """
        return Profile(self._get_freqs(1), self.alphabet)

    def get_seq_entropy(self):
        """Returns array containing Shannon entropy for each seq in self.

        Uses the profile object from get_seq_freqs (see docstring) to calculate
        the per-symbol entropy in each sequence in the alignment, i.e. the
        uncertainty about each symbol in each sequence (or row). This can be
        used to, for instance, filter low-complexity sequences.
        """
        p = self.get_seq_freqs()
        p.normalize_positions()
        return p.row_uncertainty()

    def get_pos_entropy(self):
        """Returns array containing Shannon entropy for each pos in self.

        Uses the profile object from get_pos_freqs (see docstring) to calculate
        the per-symbol entropy in each position in the alignment, i.e. the
        uncertainty about each symbol at each position (or column). This can
        be used to, for instance, detect the level of conservation at each
        position in an alignment.
        """
        p = self.get_pos_freqs()
        p.normalize_positions()
        return p.row_uncertainty()

    def iupac_consensus(self, alphabet=None):
        """Returns string containing IUPAC consensus sequence of the alignment.
        """
        if alphabet is None:
            alphabet = self.moltype
        consensus = []
        degen = alphabet.degenerate_from_seq
        for col in self.positions:
            col = alphabet.make_array_seq(col,
                                              alphabet=alphabet.alphabets.degen_gapped)
            consensus.append(degen(str(col)))
        return coerce_to_string(consensus)

    def _make_gaps_ok(self, allowed_gap_frac):
        """Makes the gaps_ok function used by omit_gap_pos and omit_gap_seqs.

        Need to make the function because if it's a method of Alignment, it
        has unwanted 'self' and 'allowed_gap_frac' parameters that impede the
        use of map() in take_seqs_if.

        WARNING: may not work correctly if component sequences have gaps that
        are not the Alignment gap character. This is because the gaps are
        checked at the column level (and the positions are lists), rather than
        at the row level. Working around this issue would probably cause a
        significant speed penalty.
        """
        def gaps_ok(seq):
            seq_len = len(seq)
            if hasattr(seq, 'count_gaps'):
                num_gaps = seq.count_gaps()
            elif hasattr(seq, 'count'):
                num_gaps = seq.count(self.alphabet.gap)
            else:
                num_gaps = sum(seq == self.alphabet.gap_index)
            return num_gaps / seq_len <= allowed_gap_frac

        return gaps_ok

    def column_freqs(self, constructor=Freqs):
        """Returns list of Freqs with item counts for each column.
        """
        return list(map(constructor, self.positions))

    def sample(self, n=None, with_replacement=False, motif_length=1,
               randint=randint, permutation=permutation):
        """Returns random sample of positions from self, e.g. to bootstrap.

        Arguments:
            - n: the number of positions to sample from the alignment.
              Default is alignment length
            - with_replacement: boolean flag for determining if sampled
              positions
            - randint and permutation: functions for random integer in a
              specified range, and permutation, respectively.


        Notes:
            By default (resampling all positions without replacement), generates
            a permutation of the positions of the alignment.

            Setting with_replacement to True and otherwise leaving parameters
            as defaults generates a standard bootstrap resampling of the
            alignment.
            """
        population_size = len(self) // motif_length
        if not n:
            n = population_size
        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            assert n <= population_size, (n, population_size, motif_length)
            locations = permutation(population_size)[:n]
        # check if we need to convert coords for multi-width motifs
        if motif_length > 1:
            locations = (locations * motif_length).repeat(motif_length)
            wrapped_locations = locations.reshape((n, motif_length))
            wrapped_locations += arange(motif_length)
        positions = take(self.array_positions, locations, 0)
        result = self.__class__(positions.T, moltype=self.moltype,
                                force_same_data=True,
                                info=self.info, names=self.names)
        return result

    def filtered(self, predicate, motif_length=1, **kwargs):
        """The alignment positions where predicate(column) is true.

        Arguments:
            - predicate: a callback function that takes an tuple of motifs and
              returns True/False
            - motif_length: length of the motifs the sequences should be split
              into, eg. 3 for filtering aligned codons."""
        length = self.seq_len
        if length % motif_length != 0:
            raise ValueError("aligned length not divisible by "
                             "motif_length=%d" % motif_length)

        num_motifs = length // motif_length
        shaped = self.array_seqs.reshape(
            (self.num_seqs, num_motifs, motif_length))
        indices = []
        elements = arange(motif_length)
        for i in range(num_motifs):
            if not predicate(shaped[:, i]):
                continue

            start = motif_length * i
            indices.extend(elements + start)

        if not indices:
            return None

        positions = self.array_seqs.take(indices, axis=1)
        result = self.__class__(positions, force_same_data=True,
                                moltype=self.moltype,
                                info=self.info, names=self.names)
        return result

    def get_gapped_seq(self, seq_name, recode_gaps=False, moltype=None):
        """Return a gapped Sequence object for the specified seqname.

        Note: return type may depend on what data was loaded into the
        SequenceCollection or Alignment.
        """
        s = self.named_seqs[seq_name]
        moltype = self.moltype if moltype is None else moltype
        if recode_gaps:
            non_ambig = list(moltype)
            ambig = moltype.degenerate_from_seq(non_ambig)
            for gapchar in moltype.gaps:
                s = s.replace(gapchar, ambig)
        return s

    def trim_stop_codons(self, gc=DEFAULT, allow_partial=False, **kwargs):
        """Removes any terminal stop codons from the sequences

        Arguments:
            - gc: genetic code object
            - allow_partial: if True and the sequence length is not divisible
              by 3, ignores the 3' terminal incomplete codon
        """
        if len(self) % 3 != 0 and not allow_partial:
            raise ValueError("alignment length not divisible by 3")

        stops = gc['*']
        get_index = self.alphabet.degen.index
        stop_indices = set(tuple(map(get_index, stop)) for stop in stops)
        new_data = self.array_seqs.copy()

        gap_indices = set(get_index(gap) for gap in self.moltype.gaps)
        gap_index = get_index(self.moltype.gap)

        trim_length = len(self)
        seqs_with_stops = 0
        for seq_num, seq in enumerate(new_data):
            # reverse find first non-gap character
            nondegen_index = None
            for i in range(len(self) - 1, -1, -1):
                e = seq[i]
                if e not in gap_indices:
                    nondegen_index = i + 1
                    if nondegen_index % 3 != 0 and not allow_partial:
                        raise ValueError(
                            "'%s' length not divisible by 3" % self.names[seq_num])
                    break

            if nondegen_index is None or nondegen_index - 3 < 0:
                continue

            # slice last three valid positions and see if stop
            end_codon = tuple(seq[nondegen_index - 3: nondegen_index])
            if end_codon in stop_indices:
                seqs_with_stops += 1
                seq[nondegen_index - 3: nondegen_index] = gap_index
                trim_length = min([trim_length, nondegen_index - 3])

        result = self.__class__(new_data.T,
                                moltype=self.moltype,
                                names=self.names,
                                info=self.info)
        # this is an ugly hack for rather odd standard behaviour
        # we find the last alignment column to have not just gap chars
        # and trim up to that
        for i in range(len(result) - 1, -1, -1):
            col = set(result.array_seqs[:, i])
            if not col <= gap_indices:
                break
        if i != len(result):
            result = result[:i + 1]

        return result

    def get_degapped_relative_to(self, name):
        """Remove all columns with gaps in sequence with given name.

        Returns Alignment object of the same class.
        Note that the seqs in the new Alignment are always new objects.

        Arguments:
            - name: sequence name
        """
        if name not in self.names:
            raise ValueError("The alignment doesn't have a sequence named '{0}'"
                             .format(name))

        gapindex = self.alphabet.index(self.alphabet.gap)
        seqindex = self.names.index(name)
        indices = self.array_seqs[seqindex] != gapindex
        new = self.array_seqs[:, indices]

        return self.__class__(new.T, names=self.names,
                              moltype=self.moltype, info=self.info)

    def add_from_ref_aln(self, ref_aln, before_name=None, after_name=None):
        """
        Insert sequence(s) to self based on their alignment to a reference
        sequence. Assumes the first sequence in ref_aln.names[0] is the
        reference.

        By default the sequence is appended to the end of the alignment,
        this can be changed by using either before_name or after_name
        arguments.

        Returns Alignment object of the same class.

        Arguments:
            - ref_aln: reference alignment (Alignment object/series) of
              reference sequence and sequences to add.
              New sequences in ref_aln (ref_aln.names[1:] are sequences to add.
              If series is used as ref_aln, it must have the structure
              [['ref_name', SEQ], ['name', SEQ]]
            - before_name: name of the sequence before which
              sequence is added
            - after_name: name of the sequence after which sequence is added
              If both before_name and after_name are specified seqs will be
              inserted using before_name.

        Example:
        Aln1:
        -AC-DEFGHI (name: seq1)
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3)

        Aln2:
        ACDEFGHI   (name: seq1)
        KL--MNPR   (name: seqX)
        KLACMNPR   (name: seqY)
        KL--MNPR   (name: seqZ)

        Out:
        -AC-DEFGHI (name: seq1)
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3)
        -KL---MNPR (name: seqX)
        -KL-ACMNPR (name: seqY)
        -KL---MNPR (name: seqZ)
        """
        # delegates to Alignment
        new = self.to_type(array_align=False)
        new = new.add_from_ref_aln(ref_aln, before_name=before_name,
                                   after_name=after_name)
        result = new.to_type(array_align=True, moltype=self.moltype)
        return result

    def replace_seqs(self, seqs, aa_to_codon=True):
        """Returns new alignment with same shape but with data taken from seqs.

        Arguments:
            - aa_to_codon: If True (default) aligns codons from protein
              alignment, or, more generally, substituting in codons from a set
              of protein sequences (not necessarily aligned). For this reason,
              it takes characters from seqs three at a time rather than one at
              a time (i.e. 3 characters in seqs are put in place of 1 character
              in self). If False, seqs must be the same lengths.

        If seqs is an alignment, any gaps in it will be ignored.
        """
        if not hasattr(seqs, 'named_seqs'):
            seqs = SequenceCollection(seqs)

        seqs = seqs.degap()

        self_gapindex = self.alphabet.index(self.alphabet.gap)
        chars_indices = seqs.alphabet.to_indices
        seq_gapindex = chars_indices(seqs.moltype.gap)[0]
        if aa_to_codon:
            scale = 3
        else:
            scale = 1

        new_dim = (len(self), scale)
        aln_length = len(self) * scale
        new_seqarr = zeros((self.num_seqs, aln_length),
                           self.array_seqs.dtype)
        assert set(self.names) == set(seqs.names), "names don't match"
        new = array([seq_gapindex] * aln_length,
                    self.array_seqs.dtype).flatten()
        for seqindex, name in enumerate(self.names):
            # convert seq to indices and make an array
            indices = chars_indices(seqs.named_seqs[name])
            orig = array(indices, self.array_seqs.dtype).flatten()
            new[:] = seq_gapindex  # reset each time through
            if scale != 1:
                if len(orig) % scale != 0:
                    raise ValueError("%s length not divisible by %s" %
                                     (name, len(orig)))

                orig.resize((len(orig) // scale, scale))
                new.resize(new_dim)

            nongap = self.array_seqs[seqindex] != self_gapindex
            try:
                new[nongap] = orig
            except ValueError:
                if nongap.sum() != orig.shape[0]:
                    raise ValueError("%s has incorrect length" % name)
                raise

            new.resize((len(self) * scale,))
            new_seqarr[seqindex] = new

        return self.__class__(new_seqarr.T, names=self.names,
                              moltype=seqs.moltype, info=self.info)

    def get_identical_sets(self, mask_degen=False):
        """returns sets of names for sequences that are identical

        Arguments:
         - mask_degen: if True, degenerate characters are ignored"""
        if mask_degen and not hasattr(self.moltype, "alphabets"):
            UserWarning("in get_identical_sets, strict has no effect as moltype "
                        "has no degenerate characters")
            mask_degen = False

        if not mask_degen:
            seqs_names = list(zip(self.array_seqs.tolist(), self.names))
            seqs_names.sort()
            matched = None
            dupes = defaultdict(set)
            for i in range(self.num_seqs - 1):
                if seqs_names[i][0] == seqs_names[i + 1][0]:
                    matched = seqs_names[i][1] if matched is None else matched
                    dupes[matched].update([seqs_names[i + 1][1], matched])
                else:
                    matched = None
            identical_sets = list(dupes.values())
            return identical_sets

        # we get the indexes range for non-degenerate characters
        indices = [self.alphabet.index(c) for c in self.moltype]
        # make sure they're all consecutive
        diffs = set([indices[i - 1] - indices[i]
                     for i in range(1, len(indices))])
        assert diffs == set([-1]), diffs
        start, end = min(indices), max(indices)
        # identify non-masked positions
        degen_masks = {}
        for i in range(self.num_seqs):
            s1 = self.array_seqs[i]
            degen_masks[i] = logical_or(s1 <= start, s1 <= end)

        identical_sets = []
        seen = set()
        for i in range(self.num_seqs - 1):
            if i in seen:
                continue

            s1 = self.array_seqs[i]
            m1 = degen_masks[i]
            group = set()
            for j in range(i + 1, self.num_seqs):
                if j in seen:
                    continue

                s2 = self.array_seqs[j]
                m2 = degen_masks[j]
                valid_pos = logical_and(m1, m2)  # non-degen positions
                if (s1[valid_pos] == s2[valid_pos]).all():
                    seen.update([i, j])
                    group.update([i, j])

            if group:
                identical_sets.append(set(self.names[i] for i in group))

        return identical_sets


class CodonArrayAlignment(ArrayAlignment):
    """Stores alignment of gapped codons, no degenerate symbols."""
    InputHandlers = {'array': aln_from_array,
                     'seqs': aln_from_array_seqs,
                     'generic': aln_from_generic,
                     'array_aln': aln_from_array_aln,
                     'aln': aln_from_collection,
                     'collection': aln_from_collection,
                     'dict': aln_from_dict,
                     'empty': aln_from_empty,
                     }


def make_gap_filter(template, gap_fraction, gap_run):
    """Returns f(seq) -> True if no gap runs and acceptable gap fraction.

    Calculations relative to template.
    gap_run = number of consecutive gaps allowed in either the template or seq
    gap_fraction = fraction of positions that either have a gap in the template
        but not in the seq or in the seq but not in the template
    NOTE: template and seq must both be ArraySequence objects.
    """
    template_gaps = array(template.gap_vector())

    def result(seq):
        """Returns True if seq adhers to the gap threshold and gap fraction."""
        seq_gaps = array(seq.gap_vector())
        # check if gap amount bad
        if sum(seq_gaps != template_gaps) / float(len(seq)) > gap_fraction:
            return False
        # check if gap runs bad
        if b'\x01' * gap_run in logical_and(seq_gaps,
                                            logical_not(template_gaps)).astype(uint8).tostring():
            return False
        # check if insertion runs bad
        elif b'\x01' * gap_run in logical_and(template_gaps,
                                              logical_not(seq_gaps)).astype(uint8).tostring():
            return False
        return True

    return result


class Alignment(_Annotatable, AlignmentI, SequenceCollection):
    moltype = None  # note: this is reset to ASCII in moltype module

    def __init__(self, *args, **kwargs):
        """Returns new Alignment object: see SequenceCollection."""

        SequenceCollection.__init__(self, *args, **kwargs)

        # need to convert seqs to Aligned objects
        seqs = self.seq_data
        names = self.names

        self._motif_probs = {}
        self._type = self.moltype.gettype()
        lengths = list(map(len, self.seq_data))
        if lengths and (max(lengths) != min(lengths)):
            raise DataError("Not all sequences are the same length:\n" +
                            "max is %s, min is %s" % (max(lengths), min(lengths)))
        aligned_seqs = []
        for s, n in zip(seqs, names):
            if isinstance(s, Aligned):
                s.name = n  # ensure consistency
                aligned_seqs.append(s)
            else:
                aligned_seqs.append(self._seq_to_aligned(s, n))
        self.named_seqs = self.AlignedSeqs = dict(
            list(zip(names, aligned_seqs)))
        self.seq_data = self._seqs = aligned_seqs

    def _coerce_seqs(self, seqs, is_array):
        if not min([isinstance(seq, _Annotatable) or isinstance(seq, Aligned) for seq in seqs]):
            seqs = list(map(self.moltype.make_seq, seqs))
        return seqs

    def _seq_to_aligned(self, seq, key):
        """Converts seq to Aligned object -- override in subclasses"""
        (map, seq) = self.moltype.make_seq(seq, key).parse_out_gaps()
        return Aligned(map, seq)

    def get_tracks(self, policy):
        # drawing code related
        # same as sequence but annotations go below sequence tracks
        return policy.tracks_for_alignment(self)

    def get_child_tracks(self, policy):
        """The only Alignment method required for cogent3.draw"""
        tracks = []
        for label in self.names:
            seq = self.named_seqs[label]
            tracks += seq.get_tracks(policy.copy(seqname=label))
        return tracks

    def __repr__(self):
        seqs = []
        limit = 10
        delimiter = ''
        for (count, name) in enumerate(self.names):
            if count == 3:
                seqs.append('...')
                break
            elts = list(str(self.named_seqs[name])[:limit + 1])
            if len(elts) > limit:
                elts.append('...')
            seqs.append("%s[%s]" % (name, delimiter.join(elts)))
        seqs = ', '.join(seqs)

        return "%s x %s %s alignment: %s" % (len(self.names),
                                             self.seq_len, self._type, seqs)

    def _mapped(self, slicemap):
        align = []
        for name in self.names:
            align.append((name, self.named_seqs[name][slicemap]))
        return self.__class__(moltype=self.moltype, data=align, info=self.info)

    def gapped_by_map(self, keep, **kwargs):
        # keep is a Map
        seqs = []
        for seq_name in self.names:
            aligned = self.named_seqs[seq_name]
            seqmap = aligned.map[keep]
            seq = aligned.data.gapped_by_map(seqmap)
            seqs.append((seq_name, seq))
        return self.__class__(moltype=self.moltype, data=seqs, **kwargs)

    def project_annotation(self, seq_name, annot):
        target_aligned = self.named_seqs[seq_name]
        if annot.parent is not self:
            raise ValueError('Annotation does not belong to this alignment')
        return annot.remapped_to(target_aligned.data, target_aligned.map)

    def get_projected_annotations(self, seq_name, *args):
        aln_annots = self.get_annotations_matching(*args)
        return [self.project_annotation(seq_name, a) for a in aln_annots]

    def get_annotations_from_seq(self, seq_name, *args):
        aligned = self.named_seqs[seq_name]
        return aligned.get_annotations_matching(self, *args)

    def get_annotations_from_any_seq(self, *args):
        result = []
        for seq_name in self.names:
            result.extend(self.get_annotations_from_seq(seq_name, *args))
        return result

    def get_by_seq_annotation(self, seq_name, *args):
        result = []
        for feature in self.get_annotations_from_seq(seq_name, *args):
            segment = self[feature.map.start:feature.map.end]
            segment.name = '%s "%s" %s to %s of %s' % (
                feature.type, feature.name,
                feature.map.start, feature.map.end, self.name or '')
            result.append(segment)
        return result

    def with_masked_annotations(self, annot_types, mask_char=None, shadow=False):
        """returns an alignment with annot_types regions replaced by mask_char
        if shadow is False, otherwise all other regions are masked.

        Arguments:
            - annot_types: annotation type(s)
            - mask_char: must be a character valid for the seq moltype. The
              default value is the most ambiguous character, eg. '?' for DNA
            - shadow: whether to mask the annotated regions, or everything but
              the annotated regions"""
        masked_seqs = []
        for seq in self.seqs:
            # we mask each sequence using these spans
            masked_seqs += [seq._masked_annotations(
                annot_types, mask_char, shadow)]
        new = self.__class__(
            data=masked_seqs, info=self.info, name=self.name)
        return new

    def filtered(self, predicate, motif_length=1, **kwargs):
        """The alignment positions where predicate(column) is true.

        Arguments:
            - predicate: a callback function that takes an tuple of motifs and
              returns True/False
            - motif_length: length of the motifs the sequences should be split
              into, eg. 3 for filtering aligned codons."""
        gv = []
        kept = False
        seqs = [self.get_gapped_seq(n).get_in_motif_size(motif_length,
                                                         **kwargs) for n in self.names]

        positions = list(zip(*seqs))
        for (position, column) in enumerate(positions):
            keep = predicate(column)
            if kept != keep:
                gv.append(position * motif_length)
                kept = keep

        if kept:
            gv.append(len(positions) * motif_length)

        if not gv:
            return None

        locations = [(gv[i], gv[i + 1]) for i in range(0, len(gv), 2)]

        keep = Map(locations, parent_length=len(self))
        return self.gapped_by_map(keep, info=self.info)

    def get_seq(self, seqname):
        """Return a ungapped Sequence object for the specified seqname.

        Note: always returns Sequence object, not ArraySequence.
        """
        return self.named_seqs[seqname].data

    def get_gapped_seq(self, seq_name, recode_gaps=False, moltype=None):
        """Return a gapped Sequence object for the specified seqname.

        Note: always returns Sequence object, not ArraySequence.
        """
        return self.named_seqs[seq_name].get_gapped_seq(recode_gaps,
                                                        moltype=moltype)

    def iter_positions(self, pos_order=None):
        """Iterates over positions in the alignment, in order.

        pos_order refers to a list of indices (ints) specifying the column
        order. This lets you rearrange positions if you want to (e.g. to pull
        out individual codon positions).

        Note that self.iter_positions() always returns new objects, by default
        lists of elements. Use map(f, self.iter_positions) to apply the
        constructor or function f to the resulting lists (f must take a single
        list as a parameter).

        Will raise IndexError if one of the indices in order exceeds the
        sequence length. This will always happen on ragged alignments:
        assign to self.seq_len to set all sequences to the same length.
        """
        get = self.named_seqs.__getitem__
        pos_order = pos_order or range(self.seq_len)
        seq_order = self.names
        aligned_objs = [get(seq) for seq in seq_order]
        seqs = list(map(str, aligned_objs))
        for pos in pos_order:
            yield [seq[pos] for seq in seqs]

    positions = property(iter_positions)

    def with_gaps_from(self, template):
        """Same alignment but overwritten with the gaps from 'template'"""
        if len(self) != len(template):
            raise ValueError("Template alignment must be same length")
        gap = self.alphabet.gap
        tgp = template.alphabet.gap
        result = {}
        for name in self.names:
            seq = self.get_gapped_seq(name)
            if name not in template.names:
                raise ValueError("Template alignment doesn't have a '%s'"
                                 % name)
            gsq = template.get_gapped_seq(name)
            assert len(gsq) == len(seq)
            combo = []
            for (s, g) in zip(seq, gsq):
                if g == tgp:
                    combo.append(gap)
                else:
                    combo.append(s)
            result[name] = combo
        return Alignment(result, alphabet=self.alphabet.with_gap_motif())

    def get_degapped_relative_to(self, name):
        """Remove all columns with gaps in sequence with given name.

        Returns Alignment object of the same class.
        Note that the seqs in the new Alignment are always new objects.

        Arguments:
            - name: sequence name
        """
        if name not in self.names:
            raise ValueError("The alignment doesn't have a sequence named '{0}'"
                             .format(name))

        gap = self.alphabet.gap
        non_gap_cols = [i for i, col in enumerate(self.get_gapped_seq(name))
                        if col != gap]

        return self.take_positions(non_gap_cols)

    def add_from_ref_aln(self, ref_aln, before_name=None, after_name=None):
        """
        Insert sequence(s) to self based on their alignment to a reference
        sequence. Assumes the first sequence in ref_aln.names[0] is the
        reference.

        By default the sequence is appended to the end of the alignment,
        this can be changed by using either before_name or after_name
        arguments.

        Returns Alignment object of the same class.

        Arguments:
            - ref_aln: reference alignment (Alignment object/series) of
              reference sequence and sequences to add.
              New sequences in ref_aln (ref_aln.names[1:] are sequences to add.
              If series is used as ref_aln, it must have the structure
              [['ref_name', SEQ], ['name', SEQ]]
            - before_name: name of the sequence before which
              sequence is added
            - after_name: name of the sequence after which sequence is added
              If both before_name and after_name are specified seqs will be
              inserted using before_name.

        Example:
        Aln1:
        -AC-DEFGHI (name: seq1)
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3)

        Aln2:
        ACDEFGHI   (name: seq1)
        KL--MNPR   (name: seqX)
        KLACMNPR   (name: seqY)
        KL--MNPR   (name: seqZ)

        Out:
        -AC-DEFGHI (name: seq1)
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3)
        -KL---MNPR (name: seqX)
        -KL-ACMNPR (name: seqY)
        -KL---MNPR (name: seqZ)
        """

        if type(ref_aln) != type(self):  # let the seq class try and guess
            ref_aln = self.__class__(ref_aln)

        ref_seq_name = ref_aln.names[0]

        if ref_seq_name not in self.names:
            raise ValueError("The name of reference sequence ({0})"
                             "not found in the alignment \n(names in the alignment:\n{1}\n)"
                             .format(ref_seq_name, "\n".join(self.names)))

        if str(ref_aln.get_gapped_seq(ref_seq_name)) \
                != str(self.get_seq(ref_seq_name)):
            raise ValueError("Reference sequences are unequal."
                             "The reference sequence must not contain gaps")

        temp_aln = None
        for seq_name in ref_aln.names[1:]:
            if seq_name in self.names:
                raise ValueError("The name of a sequence being added ({0})"
                                 "is already present".format(seq_name))

            seq = ref_aln.get_gapped_seq(seq_name)
            new_seq = Aligned(self.named_seqs[ref_seq_name].map, seq)
            if not temp_aln:
                temp_aln = self.__class__({new_seq.name: str(new_seq)})
            else:
                temp_aln = temp_aln.add_seqs(self.__class__({new_seq.name:
                                                             str(new_seq)}))

        aln = self.add_seqs(temp_aln, before_name, after_name)

        return aln

    def replace_seqs(self, seqs, aa_to_codon=True):
        """Returns new alignment with same shape but with data taken from seqs.

        Arguments:
            - aa_to_codon: If True (default) aligns codons from protein
              alignment, or, more generally, substituting in codons from a set
              of protein sequences (not necessarily aligned). For this reason,
              it takes characters from seqs three at a time rather than one at
              a time (i.e. 3 characters in seqs are put in place of 1 character
              in self). If False, seqs must be the same lengths.

        If seqs is an alignment, any gaps in it will be ignored.
        """
        if isinstance(self, ArrayAlignment):
            new = self.to_type(array_align=False)
            result = new.replace_seqs(seqs, aa_to_codon=aa_to_codon)
            return result.to_type(array_align=True)

        if aa_to_codon:
            scale = 3
        else:
            scale = 1

        if hasattr(seqs, 'named_seqs'):
            seqs = seqs.named_seqs
        else:
            seqs = SequenceCollection(seqs).named_seqs

        new_seqs = []
        for label in self.names:
            aligned = self.named_seqs[label]
            seq = seqs[label]

            if isinstance(seq, Aligned):
                seq = seq.data

            if len(seq) != len(aligned.data) * scale:
                raise ValueError("%s has incorrect length" % label)

            new_seqs.append((label, Aligned(aligned.map * scale, seq)))

        return self.__class__(new_seqs, info=self.info)


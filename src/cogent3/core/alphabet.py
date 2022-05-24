"""
alphabet.py

Contains classes for representing alphabets, and more general ordinations that
map between a set of symbols and indices for storing the results in tables.
The provided alphabets are those encountered in biological sequences, but other
alphabets are certainly possible.

WARNING: do access the standard alphabets directly. It is expected that you
will access them through the appropriate moltype. Until the moltype module
has been imported, the alphabets will not know their MolType, which will
cause problems. It is often useful to create alphabets
and/or Enumerations on the fly, however.
"""

import json

from itertools import product

import numpy

from numpy import (
    arange,
    array,
    asarray,
    frombuffer,
    newaxis,
    ravel,
    sum,
    take,
    transpose,
    uint8,
    uint16,
    uint32,
    zeros,
)

from cogent3.util.misc import get_object_provenance


Float = numpy.core.numerictypes.sctype2char(float)
Int = numpy.core.numerictypes.sctype2char(int)

__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", "Andrew Butterfield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class AlphabetError(Exception):
    pass


def get_array_type(num_elements):
    """Returns smallest array type that can contain sequence on num_elements.

    Used to figure out how large a data type is needed for the array in which
    elements are indices from an alphabet. If the data type is too small
    (e.g. you allocated an uint8 array, with 256 possible states (0-255), but
    your data actually have more than 256 states, e.g. tripeptide data with
    20*20*20 = 8000 states), when you assign a state larger than the data type
    can hold you'll get an unexpected result. For example, assigning state
    800 in an array that can only hold 256 different states will actually
    give you the result mod 256:

    >>> a = array(range(10), uint8)
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],'B')
    >>> a[0] = 800
    >>> a
    array([32,  1,  2,  3,  4,  5,  6,  7,  8,  9],'B')
           ^^
           NOTE: element 1 is _not_ 800, but instead 32 -- nasty surprise!

    Getting the size of the necessary array from the alphabet is a good
    solution to this problem.

    WARNING: Will not overflow if somehow you manage to feed it an alphabet
    with more than 2**32 elements, but it seems unlikely that this will
    happen very often in practice...
    """
    if num_elements <= 256:
        return uint8
    elif num_elements <= 2 ** 16:
        return uint16
    return uint32


def _make_translation_tables(a):
    """Makes translation tables between chars and indices.

    Return value is a tuple containing (a) the translation table where
    s.translate(a) -> array data, i.e. mapping characters to numbers, and
    the characters in the original alphabet (e.g. array of 0..4 converted
    to strings of UCAG...).

    This is useful for alphabets where the entries are all single characters
    (e.g.nucleotides or amino acids, but not codons) because we can use
    translate() on the input string to make the array of values we need instead
    of having to convert each character into a Python object and look it up in
    some mapping. Using translate() can be thousands of times faster, so it's
    almost always worth it if you have a choice.
    """
    indices = "".join(map(chr, list(range(len(a)))))
    chars = "".join(a)
    return str.maketrans(indices, chars), str.maketrans(chars, indices)


def _make_complement_array(a, complements):
    """Makes translation array between item indices and their complements."""
    comps = [complements.get(i, i) for i in a]
    return array(list(map(a.index, comps)))


class Enumeration(tuple):
    """An ordered set of objects, e.g. a list of taxon labels or sequence ids.

    An Enumeration maps items to indices, and vice versa.
    Immutable. Must initialize with a sequence of (hashable) objects, in order.

    This is the base class for alphabets. An Alphabet is a special case of
    Enumeration in which all the objects are strings of the same length.

    Stored as a tuple, but remember that if the elements in the tuple are
    mutable you can still mutate them in-place. Don't do this if you want
    your enumeration to work in a predictable fashion.

    Optionally takes a gap parameter that defines the standard gap that will
    be used for output or for operations that act on gaps. Typically, this
    will be '-' or None, depending on the application.
    """

    def __new__(cls, data=None, gap=None, moltype=None):
        """Returns a new Enumeration object.

        data can be any sequence that can be passed to the tuple() constructor.

        Takes gap as an argument but ignores it (handled in __init__).
        """
        data = data or []
        return tuple.__new__(cls, data)

    def __init__(self, data=None, gap=None, moltype=None):
        """Initializes self from data, and optionally a gap.

        An Enumeration object mainly provides the mapping between objects and
        order so that you can convert symbols on an enumeration into numeric
        indices (e.g. recoding UCAG as the numbers 0,1,2,3, or recoding
        the set of species ['Human', 'Mouse', 'Fly'] as indices 0, 1 and 2
        in a matrix.

        Properties:

        _obj_to_index: dict mapping the objects onto indices for fast lookup.

        index: provides the index of an object.

        __getitem__: provides the object at a specified index.

        shape: shape of the data, typically an n x 1 array.

        _allowed_range: stores the range in which the enumeration elements occur
        (used for summing items that match a particular symbol).

        gap: item to be used as a gap -- should typically appear in data too.

        array_type: type of array needed to store all the symbols in the
        enumeration, e.g. if your enumeration has > 256 objects in it you need
        to use uint16, not uint8, because it will wrap around otherwise. Also
        constrains the types to unsigned integer types so you don't
        accidentally use negative numbers as indices (this is very bad when
        doing indexed lookups).
        """
        self.moltype = moltype

        # check if motif lengths are homogeneous -- if so, set length
        try:
            motif_lengths = frozenset(list(map(len, self)))
            if len(motif_lengths) > 1:
                self._motiflen = None
            else:
                self._motiflen = list(motif_lengths)[0]
        except TypeError:  # some motifs don't support __len__, e.g. ints
            self._motiflen = None

        # make the quick_motifset for fast lookups; check for duplicates.
        self._quick_motifset = frozenset(self)
        if len(self._quick_motifset) != len(self):
            # got duplicates: show user what they sent in
            raise TypeError("Alphabet initialized with duplicate values:\n" + str(self))
        self._obj_to_index = dict(list(zip(self, list(range(len(self))))))
        # handle gaps
        self.gap = gap
        if gap and (gap in self):
            gap_index = self.index(gap)
            if gap_index >= 0:
                self.gap_index = gap_index
        try:
            self._gapmotif = self.gap * self._motiflen
        except TypeError:  # self._motiflen was probably None
            self._gapmotif = self.gap

        self.shape = (len(self),)
        # _allowed_range provides for fast sums of matching items
        self._allowed_range = arange(len(self))[:, newaxis]
        self.array_type = get_array_type(len(self))
        self._complement_array = None  # set in moltypes.py for standard types

    def index(self, item):
        """Returns the index of a specified item.

        This goes through an extra object lookup. If you _really_ need speed,
        you can bind self._obj_to_index.__getitem__ directly, but this is not
        recommended because the internal implementation may change."""
        return self._obj_to_index[item]

    def to_indices(self, data):
        """Returns sequence of indices from sequence of elements.

        Raises KeyError if some of the elements were not found.

        Expects data to be a sequence (e.g. list of tuple) of items that
        are in the Enumeration. Returns a list containing the index of each
        element in the input, in order.

        e.g. for the RNA alphabet ('U','C','A','G'), the sequence 'CCAU'
        would produce the result [1,1,2,0], returning the index of each
        element in the input.
        """
        return [self._obj_to_index[e] for e in data]

    def is_valid(self, seq):
        """Returns True if seq contains only items in self."""
        try:
            self.to_indices(seq)
            return True
        except (KeyError, TypeError):
            return False

    def from_indices(self, data):
        """Returns sequence of elements from sequence of indices.

        Specifically, takes as input a sequence of numbers corresponding to
        elements in the Enumeration (i.e. the numbers must all be < len(self).
        Returns a list of the items in the same order as the indices. Inverse
        of to_indices.

        e.g. for the DNA alphabet ('U','C','A','G'), the sequence [1,1,2,0]
        would produce the result 'CCAU', returning the element corresponding
        to each element in the input.

        """
        # if it's a normal Python type, map will work
        return [self[index] for index in data]

    def __pow__(self, num):
        """Returns JointEnumeration with num copies of self.

        A JointEnumeration is an Enumeration of tuples on the original
        enumeration (although these may be mapped to a different data type,
        e.g. a JointAlphabet is still an Alphabet, so its members are
        fixed-length strings).

        For example, a trinucleotide alphabet (or codons) would be a
        JointEnumeration on the nucleotides. So RnaBases**3 is RnaCodons (except
        for some additional logic for converting between tuples of one-letter
        strings and the desired 3-letter strings). All subenumerations of a
        JointEnumeration made by __pow__ are identical.
        """
        return JointEnumeration([self] * num, moltype=self.moltype)

    def __mul__(self, other):
        """Returns JointEnumeration between self and other.

        Specifically, returns a JointEnumeration whose elements are (a,b) where
        a is an element of the first enumeration and b is an element of the
        second enumeration. For example, a JointEnumeration of 'ab' and 'cd'
        would have the four elements ('a','c'),('a','d'),('b','c'),('b',d').

        A JointEnumeration is an enumeration of tuples on more than one
        enumeration, where the first element in each tuple comes from the
        first enumeration, the second from the second enumeration, and so on.
        JointEnumerations are useful as the basis for contingency tables,
        transition matrices, counts of dinucleotides, etc.
        """
        if self.moltype is other.moltype:
            moltype = self.moltype
        else:
            moltype = None
        return JointEnumeration([self, other], moltype=moltype)

    def counts(self, a):
        """Returns array containing counts of each item in a.

        For example, on the enumeration 'UCAG', the sequence 'CCUG' would
        return the array [1,2,0,1] reflecting one count for the first item
        in the enumeration ('U'), two counts for the second item ('C'), no
        counts for the third item ('A'), and one count for the last item ('G').

        The result will always be a vector of Int with length equal to
        the length of the enumeration. We return Int and non an unsigned
        type because it's common to subtract counts, which produces surprising
        results on unit types (i.e. wrapraround to maxint) unless the type
        is explicitly coerced by the user.

        Sliently ignores any unrecognized indices, e.g. if your enumeration
        contains 'TCAG' and you get an 'X', the 'X' will be ignored because
        it has no index in the enumeration.
        """
        try:
            data = ravel(a)
        except ValueError:  # ravel failed; try coercing to array
            try:
                data = ravel(array(a))
            except ValueError:  # try mapping to string
                data = ravel(array(list(map(str, a))))
        return sum(asarray(self._allowed_range == data, Int), axis=-1)

    def _get_pairs(self):
        """Accessor for pairs, lazy evaluation."""
        if not hasattr(self, "_pairs"):
            self._pairs = self ** 2
        return self._pairs

    pairs = property(_get_pairs)

    def _get_triples(self):
        """Accessor for triples, lazy evaluation."""
        if not hasattr(self, "_triples"):
            self._triples = self ** 3
        return self._triples

    Triples = property(_get_triples)


class JointEnumeration(Enumeration):
    """Holds an enumeration composed of subenumerations. Immutable.

    JointEnumeration[i] will return tuple of items on each of the constituent
    alphabets. For example, a JointEnumeration between the enumerations 'ab'
    and 'cd' would have four elements: ('a','c'),('a','d'),('b','c'),('b','d').
    (note that if doing a JointAlphabet, these would be strings, not tuples).
    Note that the two enumerations do not have to be the same, although it is
    often convenient if they are (e.g. pair enumerations that underlie
    substitution matrices).
    """

    def __new__(cls, data=None, gap=None, moltype=None):
        """Fills in the tuple with tuples from the enumerations in data."""
        data = data or []
        sub_enums = cls._coerce_enumerations(data)
        return Enumeration.__new__(cls, product(*sub_enums), moltype=moltype)

    def __init__(self, data=None, gap=None, moltype=None):
        """Returns a new JointEnumeration object. See class docstring for info.

        Expects a list of Enumeration objects, or objects that can be coerced
        into Enumeration objects (basically, anything that can be a tuple).

        Does NOT have an independent concept of a gap -- gets the gaps from the
        constituent subenumerations.
        """
        data = data or []
        self.sub_enumerations = self._coerce_enumerations(data)
        sub_enum_lengths = list(map(len, self.sub_enumerations))
        # build factors for combining symbols.
        curr_factor = 1
        sub_enum_factors = [curr_factor]
        for i in sub_enum_lengths[-1:0:-1]:
            curr_factor *= i
            sub_enum_factors = [curr_factor] + sub_enum_factors
        self._sub_enum_factors = transpose(array([sub_enum_factors]))

        try:
            # figure out the gaps correctly
            gaps = [i.gap for i in self.sub_enumerations]
            self.gap = tuple(gaps)
            gap_indices = array([i.gap_index for i in self.sub_enumerations])
            gap_indices *= sub_enum_factors
            self.gap_index = sum(gap_indices)
        except (TypeError, AttributeError):  # index not settable
            self.gap = None

        super(JointEnumeration, self).__init__(self, self.gap)
        # remember to reset shape after superclass init
        self.shape = tuple(sub_enum_lengths)

    def __getnewargs_ex__(self, *args, **kw):
        data = self.to_rich_dict(for_pickle=True)
        r = tuple([data[k] for k in ("data", "gap", "moltype")])
        return r, {}

    def to_rich_dict(self, for_pickle=False):
        data = {
            "data": ["".join(e) for e in self.sub_enumerations],
            "gap": self.gap,
            "moltype": self.moltype,
        }
        if not for_pickle:
            data["type"] = get_object_provenance(self)
            data["moltype"] = self.moltype.label
            data["version"] = __version__
        return data

    def to_json(self):
        """returns result of json formatted string"""
        data = self.to_rich_dict(for_pickle=False)
        return json.dumps(data)

    def _coerce_enumerations(cls, enums):
        """Coerces putative enumerations into Enumeration objects.

        For each object passed in, if it's an Enumeration object already, use
        that object without translation/conversion. If it isn't, call the
        Enumeration constructor on it and append the new Enumeration to the
        result.

        Note that this means you can construct JointEnumerations where the
        subenumerations have the same data but are different objects --
        in general, you probably don't want to do this (i.e. you should make
        it into an Enumeration beforehand and pass n references to that
        Enumeration in as a list:
            a = Enumeration('abc')
            j = JointEnumeration([a,a,a])
        ... not
            a = JointEnumeration(['abc','abc','abc'])
        """
        result = []
        for a in enums:
            if isinstance(a, Enumeration):
                result.append(a)
            else:
                result.append(Enumeration(a))
        return result

    def pack_arrays(self, arrays):
        """Packs parallel arrays on subenums to single array on joint enum.

        WARNING: must pass arrays as single object.

        This method takes a single array in which each row is an array of
        indices on the appropriate subenumeration. For example, you might have
        arrays for the bases at the first, second, and third positions of each
        codon in a gene, and want to pack them together into a single Codons
        object. pack_arrays() allows you to do this without having to explicitly
        interleave the arrays into a single sequence and then convert it back
        on the JointEnumeration.

        Notes:

        - Expects a single array object, where the rows (first dimension)
        correspond to information about the same set of data on a different
        enumeration (or, as in the case of codons, on the same enumeration but
        at a different position). This means that if you're constructing the
        array on the fly, the number of elements you have in each enumeration
        must be the same.

        - Arrays must already be converted into indices -- for example, you
        can't pass in raw strings as sequences.

        - This method is the inverse of unpack_arrays().

        - Uses self.array_type to figure out the type of array to return (e.g.
        the amino acids may use a character array, but you need a larger
        data type to store indices on a JointEnumeration of pairs or triples of
        amino acids).
        """
        return sum(self._sub_enum_factors * array(arrays, self.array_type), axis=0)

    def unpack_arrays(self, a):
        """Unpacks array on joint enum to individual arrays on subenums.

        Returns result as single numpy array object.

        This method takes a single vector of indices on the appropriate
        JointEnumeration, and returns an array where the rows are, in order,
        vectors of the appropriate indices on each subenumeration. For example,
        you might have a sequence on the Codons enumeration, and want to unpack
        it into the corresponding sequences of the first, second, and third
        position bases in each codon. unpack_arrays() allows you to do this .

        Notes:

        - Will always return a single array object, with number of rows equal
        to the number of subenumerations in self.

        - Will always return a value for each enumeration for each item, e.g.
        a sequence on the codon enumeration will always return three sequences
        on the individual enumerations that all have the same length (packed
        into a single array).

        - Output will aill always use the same typecode as the input array.
        """
        a = array(a)
        num_enums = len(self.sub_enumerations)
        result = zeros((num_enums, len(a)))
        lengths = self.shape
        # count backwards through the enumeration lengths and add to array
        for i in range(num_enums - 1, -1, -1):
            length = lengths[i]
            result[i] = a % length
            a //= array(length, a.dtype.char)
        return result

    # the following, _coerce_enumerations, is a class method because we use
    # it in __new__ before we have an instance to call it on.
    _coerce_enumerations = classmethod(_coerce_enumerations)


class Alphabet(Enumeration):
    """An ordered set of fixed-length strings, e.g. the 61 sense codons.

    ambiguities (e.g. N for any base in DNA) are not considered part of the
    alphabet itself, although a sequence is valid on the alphabet even if
    it contains ambiguities that are known to the alphabet.
    A gap is considered a separate motif and is not part of the alphabet itself.

    The typical use is for the Alphabet to hold nucleic acid bases, amino acids,
    or codons.

    The moltype, if supplied, handles ambiguities, coercion of the sequence
    to the correct data type, and complementation (if appropriate).
    """

    # make this exception avalable to objects calling alphabet methods.
    AlphabetError = AlphabetError

    def __new__(cls, motifset, gap="-", moltype=None):
        """Returns a new Alphabet object."""
        return Enumeration.__new__(cls, data=motifset, gap=gap, moltype=moltype)

    def __init__(self, motifset, gap="-", moltype=None):
        """Returns a new Alphabet object."""
        super(Alphabet, self).__init__(data=motifset, gap=gap, moltype=moltype)

    def __getnewargs_ex__(self, *args, **kw):
        data = self.to_rich_dict(for_pickle=True)
        r = tuple([data[k] for k in ("motifset", "gap", "moltype")])
        return r, {}

    def to_rich_dict(self, for_pickle=False):
        data = {"motifset": tuple(self), "gap": self.gap, "moltype": self.moltype}
        if not for_pickle:
            data["type"] = get_object_provenance(self)
            data["moltype"] = self.moltype.label
            data["version"] = __version__
            if hasattr(self, "_gc"):
                data["genetic_code"] = self._gc.name
        return data

    def to_json(self):
        """returns result of json formatted string"""
        data = self.to_rich_dict(for_pickle=False)
        return json.dumps(data)

    def get_word_alphabet(self, word_length):
        """Returns a new Alphabet object with items as word_length strings.

        Note that the result is not a JointEnumeration object, and cannot
        unpack its indices. However, the items in the result _are_ all strings.
        """
        crossproduct = [""]
        for a in range(word_length):
            n = []
            for c in crossproduct:
                for m in self:
                    n.append(m + c)
            crossproduct = n
        return Alphabet(crossproduct, moltype=self.moltype)

    def from_seq_to_array(self, sequence):
        """Returns an array of indices corresponding to items in sequence.

        Parameters
        ----------
        sequence: Sequence
         A cogent3 sequence object

        Returns
        -------
        ndarray

        Notes
        -----
        Unlike to_indices() in superclass, this method returns a numpy array
        object. It also breaks the seqeunce into items in the current alphabet
        (e.g. breaking a raw DNA sequence into codons), which to_indices() does
        """
        sequence = sequence.get_in_motif_size(self._motiflen)
        return array(list(map(self.index, sequence)))

    def from_ordinals_to_seq(self, data):
        """Returns a Sequence object corresponding to indices in data.

        Parameters
        ----------
        data: series
            series of int

        Returns
        -------
        Sequence with self.moltype

        Notes
        -----
        Unlike from_indices(), this method uses the MolType to
        coerce the result into a sequence of the correct class.

        Raises an AttributeError if MolType is not set.
        """
        return self.moltype.make_seq("".join(self[i] for i in data))

    def get_matched_array(self, motifs, dtype=Float):
        """Returns an array in which rows are motifs, columns are items in self.

        Result is an array of Float in which a[i][j] indicates whether the ith
        motif passed in as motifs is a symbol that matches the jth character
        in self. For example, on the DNA alphabet 'TCAG', the degenerate symbol
        'Y' would correspond to the row [1,1,0,0] because Y is a degenerate
        symbol that encompasses T and C but not A or G.

        This code is similar to code in the Profile class, and should perhaps
        be merged with it (in particular, because there is nothing likelihood-
        specific about the resulting match table).
        """
        result = zeros([len(motifs), len(self)], dtype)
        obj_to_index = self._obj_to_index
        for (u, ambig_motif) in enumerate(motifs):
            for motif in self.resolve_ambiguity(ambig_motif):
                result[u, obj_to_index[motif]] = 1.0
        return result

    def get_motif_len(self):
        """Returns the length of the items in self, or None if they differ."""
        return self._motiflen

    def get_gap_motif(self):
        """Returns the motif that self is using as a gap. Note that this will
        typically be a multiple of self.gap.
        """
        return self._gapmotif

    def includes_gap_motif(self):
        """Returns True if self includes the gap motif, False otherwise."""
        return self._gapmotif in self

    def _with(self, motifset):
        """Returns a new Alphabet object with same class and moltype as self.

        Will always return a new Alphabet object even if the motifset is the
        same.
        """
        return self.__class__(tuple(motifset), moltype=self.moltype)

    def with_gap_motif(self):
        """Returns an Alphabet object resembling self but including the gap.

        Always returns the same object.
        """
        if self.includes_gap_motif():
            return self
        if not hasattr(self, "gapped"):
            self.gapped = self._with(list(self) + [self.get_gap_motif()])
        return self.gapped

    def get_subset(self, motif_subset, excluded=False):
        """Returns a new Alphabet object containing a subset of motifs in self.

        Raises an exception if any of the items in the subset are not already
        in self. Always returns a new object.
        """
        if isinstance(motif_subset, dict):
            motif_subset = [m for m in motif_subset if motif_subset[m]]
        for m in motif_subset:
            if m not in self:
                raise AlphabetError(m)
        if excluded:
            motif_subset = [m for m in self if m not in motif_subset]
        return self._with(motif_subset)

    def resolve_ambiguity(self, ambig_motif):
        """Returns set of symbols corresponding to ambig_motif.

        Handles multi-character symbols and screens against the set of
        valid motifs, unlike the MolType version.
        """
        # shortcut easy case
        if ambig_motif in self._quick_motifset:
            return (ambig_motif,)

        # resolve each letter, and build the possible sub motifs
        ambiguities = self.moltype.ambiguities
        motif_set = [""]
        ALL = self.moltype.alphabet.with_gap_motif()
        for character in ambig_motif:
            new_motifs = []
            if character == "?":
                resolved = ALL
            elif character == "-":
                resolved = ["-"]
            else:
                try:
                    resolved = ambiguities[character]
                except KeyError:
                    raise AlphabetError(ambig_motif)
            for character2 in resolved:
                for motif in motif_set:
                    new_motifs.append("".join([motif, character2]))

            motif_set = new_motifs

        # delete sub motifs that are not to be included
        motif_set = [motif for motif in motif_set if motif in self._quick_motifset]

        if not motif_set:
            raise AlphabetError(ambig_motif)

        return tuple(motif_set)

    def adapt_motif_probs(self, motif_probs):
        """Prepare an array or dictionary of probabilities for use with
        this alphabet by checking size and order"""
        if hasattr(motif_probs, "keys"):
            sample = list(motif_probs.keys())[0]
            if sample not in self:
                raise ValueError(f"Can't find motif {sample} in alphabet")
            motif_probs = numpy.array([motif_probs[motif] for motif in self])
        else:
            if len(motif_probs) != len(self):
                if len(motif_probs) != len(self):
                    raise ValueError(
                        f"Can't match {len(motif_probs)} probs to {len(self)} alphabet"
                    )
            motif_probs = numpy.asarray(motif_probs)
        assert abs(sum(motif_probs) - 1.0) < 0.0001, motif_probs
        return motif_probs


class CharAlphabet(Alphabet):
    """Holds an alphabet whose items are single chars.

    The general Alphabet can hold items of any type, but this is inconvenient
    if your Alphabet is characters-only because you get back operations on
    the alphabet as tuples of single-character strings instead of the strings
    you probably want. Having a separate represntation for CharAlphabets
    also allows certain efficiencies, such as using translation tables to
    map characters and indices instead of having to extract each element
    searately for remapping.
    """

    def __init__(self, data=None, gap="-", moltype=None):
        """Initializes self from items.

        data should be a sequence (string, list, etc.) of characters that
        are in the alphabet, e.g. 'UCAG' for RNA.

        gap should be a single character that represents the gap, e.g. '-'.
        """
        data = data or []
        super(CharAlphabet, self).__init__(data, gap, moltype=moltype)
        self._indices_to_chars, self._chars_to_indices = _make_translation_tables(data)
        self._char_nums_to_indices = array(range(256), uint8)
        for c, i in self._chars_to_indices.items():
            self._char_nums_to_indices[c] = i

        chars = bytearray(range(256))
        for i, c in self._indices_to_chars.items():
            chars[i] = c
        self._indices_nums_to_chars = array(list(chars), "B").view("c")

    def from_string(self, data):
        """Returns array of indices from string containing elements.

        data should be a string on the alphabet, e.g. 'ACC' for the RNA
        alhabet 'UCAG' would return the array [2,1,1]. This is useful for
        converting strings into arrays of small integers on the alphabet,
        e.g. for reading a Sequence from a string.

        This is on the Alphabet, not the Sequence, because lots of objects
        (e.g. Profile, Alignment) also need to use it.
        """
        vals = str.translate(data, self._chars_to_indices)
        vals = frombuffer(memoryview(vals.encode("utf8")), dtype=uint8)
        return vals

    def is_valid(self, seq):
        """Returns True if seq contains only items in self."""
        try:
            if len(seq) == 0:  # can't be invalid if empty
                return True
            ind = self.to_indices(seq)
            return max(ind) < len(self) and min(ind) >= 0
        except (TypeError, KeyError):
            return False

    def from_array(self, data):
        """Returns array of indices from array containing elements.

        This is useful if, instead of a string, you have an array of
        characters that's been converted into a numpy array. See
        from_string docstring for general behavior.
        """
        return take(self._char_nums_to_indices, data.view("B"))

    def to_chars(self, data):
        """Converts array of indices into array of elements.

        For example, on the 'UCAG' RNA alphabet, an array with the data
        [0,1,1] would return the characters [U,C,C] in a byte array.
        """
        data = array(data)
        return take(self._indices_nums_to_chars, data.astype("B"))

    def to_string(self, data, delimiter="\n"):
        """Converts array of data into string.

        For example, on the 'UCAG' RNA alphabet, an array with the data
        [0,1,1] would return the string 'UCC'. This is the most useful
        conversion mechanism, and is used by e.g. the Sequence object to
        convert the internal data back into strings for output.
        """
        s = data.shape
        if not s:
            return ""
        elif len(s) == 1:
            val = self.to_chars(data)
            val = val.tobytes().decode("utf-8")
            return val
        else:
            return delimiter.join(
                [i.tobytes().decode("utf-8") for i in self.to_chars(data)]
            )

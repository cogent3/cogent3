#!/usr/bin/env python
"""
alphabet.py

Contains classes for representing alphabets, and more general ordinations that
map between a set of symbols and indices for storing the results in tables.
The provided alphabets are those encountered in biological sequences, but other
alphabets are certainly possible.

WARNING: do access the standard Alphabets directly. It is expected that you
will access them through the appropriate MolType. Until the moltype module
has been imported, the Alphabets will not know their MolType, which will
cause problems. It is often useful to create Alphabets
and/or Enumerations on the fly, however.

MolType provides services for resolving ambiguities, or providing the
correct ambiguity for recoding -- will move to its own module.
"""

from cogent.util.array import cartesian_product

import re
import string
from numpy import array, sum, transpose, remainder, zeros, arange, newaxis, \
    ravel, asarray, fromstring, take, uint8, uint16, uint32, take
from string import maketrans, translate
import numpy
Float = numpy.core.numerictypes.sctype2char(float)
Int = numpy.core.numerictypes.sctype2char(int)

__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Andrew Butterfield"]
__license__ = "GPL"
__version__ = "1.5.3"
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
    
    Getting the size of the necessary array from the Alphabet is a good
    solution to this problem.
    
    WARNING: Will not overflow if somehow you manage to feed it an alphabet
    with more than 2**32 elements, but it seems unlikely that this will
    happen very often in practice...
    """
    if num_elements <= 256:
        return uint8
    elif num_elements <= 2**16:
        return uint16
    return uint32

def _make_translation_tables(a):
    """Makes translation tables between chars and indices.
    
    Return value is a tuple containing (a) the translation table where
    s.translate(a) -> array data, i.e. mapping characters to numbers, and
    (b) the translation table where a.tostring().translate(s) -> string of
    the characters in the original alphabet (e.g. array of 0..4 converted
    to strings of UCAG...).
    
    This is useful for alphabets where the entries are all single characters
    (e.g.nucleotides or amino acids, but not codons) because we can use
    translate() on the input string to make the array of values we need instead
    of having to convert each character into a Python object and look it up in
    some mapping. Using translate() can be thousands of times faster, so it's
    almost always worth it if you have a choice.
    """
    indices = ''.join(map(chr, range(len(a))))
    chars = ''.join(a)
    return maketrans(indices, chars), maketrans(chars, indices)

def _make_complement_array(a, complements):
    """Makes translation array between item indices and their complements."""
    comps = [complements.get(i, i) for i in a]
    return array(map(a.index, comps))

class Enumeration(tuple):
    """An ordered set of objects, e.g. a list of taxon labels or sequence ids.
    
    An Enumeration maps items to indices, and vice versa.
    Immutable. Must initialize with a sequence of (hashable) objects, in order.
    
    This is the base class for Alphabets. An Alphabet is a special case of
    Enumeration in which all the objects are strings of the same length.
    
    Stored as a tuple, but remember that if the elements in the tuple are
    mutable you can still mutate them in-place. Don't do this if you want
    your enumeration to work in a predictable fashion.
    
    Optionally takes a Gap parameter that defines the standard gap that will
    be used for output or for operations that act on gaps. Typically, this
    will be '-' or None, depending on the application.
    """
    
    def __new__(cls, data=[], Gap=None, MolType=None):
        """Returns a new Enumeration object.
        
        data can be any sequence that can be passed to the tuple() constructor.
        
        Takes Gap as an argument but ignores it (handled in __init__).
        """
        return tuple.__new__(cls,data)
    
    def __init__(self, data=[], Gap=None, MolType=None):
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
        
        Shape: shape of the data, typically an n x 1 array.
        
        _allowed_range: stores the range in which the enumeration elements occur
        (used for summing items that match a particular symbol).
        
        Gap: item to be used as a gap -- should typically appear in data too.
        
        ArrayType: type of array needed to store all the symbols in the
        enumeration, e.g. if your enumeration has > 256 objects in it you need
        to use uint16, not uint8, because it will wrap around otherwise. Also
        constrains the types to unsigned integer types so you don't
        accidentally use negative numbers as indices (this is very bad when
        doing indexed lookups).
        """
        self.MolType = MolType
        
        #check if motif lengths are homogeneous -- if so, set length
        try:
            motif_lengths = frozenset(map(len, self))
            if len(motif_lengths) > 1:
                self._motiflen = None
            else:
                self._motiflen = list(motif_lengths)[0]
        except TypeError:       #some motifs don't support __len__, e.g. ints
            self._motiflen = None
        
        #make the quick_motifset for fast lookups; check for duplicates.
        self._quick_motifset = frozenset(self)
        if len(self._quick_motifset) != len(self):
            #got duplicates: show user what they sent in
            raise TypeError, 'Alphabet initialized with duplicate values:\n' +\
                str(self)
        self._obj_to_index = dict(zip(self, range(len(self))))
        #handle gaps
        self.Gap = Gap
        if Gap and (Gap in self):
            gap_index = self.index(Gap)
            if gap_index >= 0:
                self.GapIndex = gap_index
        try:
            self._gapmotif = self.Gap * self._motiflen
        except TypeError:       #self._motiflen was probably None
            self._gapmotif = self.Gap
        
        self.Shape = (len(self),)
        #_allowed_range provides for fast sums of matching items
        self._allowed_range = arange(len(self))[:,newaxis]
        self.ArrayType = get_array_type(len(self))
        self._complement_array = None   #set in moltypes.py for standard types
    
    def index(self, item):
        """Returns the index of a specified item.
        
        This goes through an extra object lookup. If you _really_ need speed,
        you can bind self._obj_to_index.__getitem__ directly, but this is not
        recommended because the internal implementation may change."""
        return self._obj_to_index[item]
    
    def toIndices(self, data):
        """Returns sequence of indices from sequence of elements.
        
        Raises KeyError if some of the elements were not found.
        
        Expects data to be a sequence (e.g. list of tuple) of items that
        are in the Enumeration. Returns a list containing the index of each
        element in the input, in order.
        
        e.g. for the RNA alphabet ('U','C','A','G'), the sequence 'CCAU'
        would produce the result [1,1,2,0], returning the index of each
        element in the input.
        """
        return map(self._obj_to_index.__getitem__, data)
    
    def isValid(self, seq):
        """Returns True if seq contains only items in self."""
        try:
            self.toIndices(seq)
            return True
        except (KeyError, TypeError):
            return False
    
    def fromIndices(self, data):
        """Returns sequence of elements from sequence of indices.
        
        Specifically, takes as input a sequence of numbers corresponding to
        elements in the Enumeration (i.e. the numbers must all be < len(self).
        Returns a list of the items in the same order as the indices. Inverse
        of toIndices.
        
        e.g. for the DNA alphabet ('U','C','A','G'), the sequence [1,1,2,0]
        would produce the result 'CCAU', returning the element corresponding
        to each element in the input.
        
        """
        #if it's a normal Python type, map will work
        try:
            return map(self.__getitem__, data)
        #otherwise, it's probably an array object.
        except TypeError:
            try:
                data = map(int, data)
            except (TypeError, ValueError): #might be char array?
                print "DATA", data
                print "FIRST MAP:", map(str, data)
                print "SECOND MAP:", map(ord, map(str, data))
                data = map(ord, map(str, data))
            return(map(self.__getitem__, data))
    
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
        return JointEnumeration([self]*num, MolType=self.MolType)
    
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
        if self.MolType is other.MolType:
            MolType = self.MolType
        else:
            MolType = None
        return JointEnumeration([self, other], MolType=MolType)
    
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
        except ValueError:  #ravel failed; try coercing to array
            try:
                data = ravel(array(a))
            except ValueError: #try mapping to string
                data = ravel(array(map(str, a)))
        return sum(asarray(self._allowed_range == data, Int), axis=-1)
    
    def _get_pairs(self):
        """Accessor for pairs, lazy evaluation."""
        if not hasattr(self, '_pairs'):
            self._pairs = self**2
        return self._pairs
    
    Pairs = property(_get_pairs)
    
    def _get_triples(self):
        """Accessor for triples, lazy evaluation."""
        if not hasattr(self, '_triples'):
            self._triples = self**3
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
    def __new__(cls, data=[], Gap=None, MolType=None):
        """Fills in the tuple with tuples from the enumerations in data."""
        sub_enums = cls._coerce_enumerations(data)
        return Enumeration.__new__(cls, cartesian_product(sub_enums), \
            MolType=MolType)
    
    def __init__(self, data=[], Gap=None, MolType=None):
        """Returns a new JointEnumeration object. See class docstring for info.
        
        Expects a list of Enumeration objects, or objects that can be coerced
        into Enumeration objects (basically, anything that can be a tuple).
        
        Does NOT have an independent concept of a gap -- gets the gaps from the
        constituent subenumerations.
        """
        self.SubEnumerations = self._coerce_enumerations(data)
        sub_enum_lengths = map(len, self.SubEnumerations)
        #build factors for combining symbols.
        curr_factor = 1
        sub_enum_factors = [curr_factor]
        for i in sub_enum_lengths[-1:0:-1]:
            curr_factor *= i
            sub_enum_factors = [curr_factor] + sub_enum_factors
        self._sub_enum_factors = transpose(array([sub_enum_factors]))
        
        try:
            #figure out the gaps correctly
            gaps = [i.Gap for i in self.SubEnumerations]
            self.Gap = tuple(gaps)
            gap_indices = array([i.GapIndex for i in self.SubEnumerations])
            gap_indices *= sub_enum_factors
            self.GapIndex = sum(gap_indices)
        except (TypeError, AttributeError): #index not settable
            self.Gap = None
        
        super(JointEnumeration, self).__init__(self, self.Gap)
        #remember to reset shape after superclass init
        self.Shape = tuple(sub_enum_lengths)
    
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
    
    def packArrays(self, arrays):
        """Packs parallel arrays on subenums to single array on joint enum.
        
        WARNING: must pass arrays as single object.
        
        This method takes a single array in which each row is an array of
        indices on the appropriate subenumeration. For example, you might have
        arrays for the bases at the first, second, and third positions of each
        codon in a gene, and want to pack them together into a single Codons
        object. packArrays() allows you to do this without having to explicitly
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
        can't pass in raw strings as sequences, but you can pass in the data
        from cogent.seqsim.Sequence objects.
        
        - This method is the inverse of unpackArrays().
        
        - Uses self.ArrayType to figure out the type of array to return (e.g.
        the amino acids may use a character array, but you need a larger
        data type to store indices on a JointEnumeration of pairs or triples of
        amino acids).
        """
        return sum(self._sub_enum_factors * array(arrays, self.ArrayType), axis=0)
    
    def unpackArrays(self, a):
        """Unpacks array on joint enum to individual arrays on subenums.
        
        Returns result as single numpy array object.
        
        This method takes a single vector of indices on the appropriate
        JointEnumeration, and returns an array where the rows are, in order,
        vectors of the appropriate indices on each subenumeration. For example,
        you might have a sequence on the Codons enumeration, and want to unpack
        it into the corresponding sequences of the first, second, and third
        position bases in each codon. unpackArrays() allows you to do this .
        
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
        num_enums = len(self.SubEnumerations)
        result = zeros((num_enums, len(a)))
        lengths = self.Shape
        # count backwards through the enumeration lengths and add to array
        for i in range(num_enums-1, -1, -1):
            length = lengths[i]
            result[i] = a % length
            a /= array(length,a.dtype.char)
        return result
    
    # the following, _coerce_enumerations, is a class method because we use
    # it in __new__ before we have an instance to call it on.
    _coerce_enumerations = classmethod(_coerce_enumerations)

class Alphabet(Enumeration):
    """An ordered set of fixed-length strings, e.g. the 61 sense codons.
    
    Ambiguities (e.g. N for any base in DNA) are not considered part of the
    alphabet itself, although a sequence is valid on the alphabet even if
    it contains ambiguities that are known to the alphabet.
    A gap is considered a separate motif and is not part of the alphabet itself.
    
    The typical use is for the Alphabet to hold nucleic acid bases, amino acids,
    or codons.
    
    The MolType, if supplied, handles ambiguities, coercion of the sequence
    to the correct data type, and complementation (if appropriate).
    """
    
    # make this exception avalable to objects calling alphabet methods.
    AlphabetError = AlphabetError
    
    def __new__(cls, motifset, Gap='-', MolType=None):
        """Returns a new Alphabet object."""
        return Enumeration.__new__(cls, data=motifset, Gap=Gap, \
            MolType=MolType)
    
    def __init__(self, motifset, Gap='-', MolType=None):
        """Returns a new Alphabet object."""
        super(Alphabet, self).__init__(data=motifset, Gap=Gap, \
            MolType=MolType)
    
    def getWordAlphabet(self, length):
        """Returns a new Alphabet object with items as length-n strings.
        
        Note that the result is not a JointEnumeration object, and cannot
        unpack its indices. However, the items in the result _are_ all strings.
        """
        crossproduct = ['']
        for a in range(length):
            n = []
            for c in crossproduct:
                for m in self:
                    n.append(m+c)
            crossproduct = n
        return Alphabet(crossproduct, MolType=self.MolType)
    
    def fromSequenceToArray(self, sequence):
        """Returns an array of indices corresponding to items in sequence.
        
        Unlike toIndices() in superclass, this method returns a numpy array
        object. It also breaks the seqeunce into items in the current alphabet
        (e.g. breaking a raw DNA sequence into codons), which toIndices() does
        not do. It also requires the sequence to be a Sequence object rather
        than an arbitrary string, tuple, etc.
        """
        sequence = sequence.getInMotifSize(self._motiflen)
        return array(map(self.index, sequence))
    
    def fromOrdinalsToSequence(self, data):
        """Returns a Sequence object corresponding to indices in data.
        
        Unlike fromIndices() in superclass, this method uses the MolType to
        coerce the result into a sequence of the correct class. Note that if
        the MolType is not set, this method will raise an AttributeError.
        """
        result = ''
        return self.MolType.makeSequence(''.join(self[i] for i in data))
    
    def fromAmbigToLikelihoods(self, motifs, dtype=Float):
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
            for motif in self.resolveAmbiguity(ambig_motif):
                result[u, obj_to_index[motif]] = 1.0
        return result
    
    def getMotifLen(self):
        """Returns the length of the items in self, or None if they differ."""
        return self._motiflen
    
    def getGapMotif(self):
        """Returns the motif that self is using as a gap. Note that this will
        typically be a multiple of self.Gap.
        """
        return self._gapmotif
    
    def includesGapMotif(self):
        """Returns True if self includes the gap motif, False otherwise."""
        return self._gapmotif in self
    
    def _with(self, motifset):
        """Returns a new Alphabet object with same class and moltype as self.
        
        Will always return a new Alphabet object even if the motifset is the
        same.
        """
        return self.__class__(tuple(motifset), MolType=self.MolType)
    
    def withGapMotif(self):
        """Returns an Alphabet object resembling self but including the gap.
        
        Always returns the same object.
        """
        if self.includesGapMotif():
            return self
        if not hasattr(self, 'Gapped'):
            self.Gapped = self._with(list(self) + [self.getGapMotif()])
        return self.Gapped
    
    def getSubset(self, motif_subset, excluded=False):
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
    
    def resolveAmbiguity(self, ambig_motif):
        """Returns set of symbols corresponding to ambig_motif.
        
        Handles multi-character symbols and screens against the set of
        valid motifs, unlike the MolType version.
        """
        # shortcut easy case
        if ambig_motif in self._quick_motifset:
            return (ambig_motif,)
        
        # resolve each letter, and build the possible sub motifs
        ambiguities = self.MolType.Ambiguities
        motif_set = ['']
        ALL = self.MolType.Alphabet.withGapMotif()
        for character in ambig_motif:
            new_motifs = []
            if character == '?':
                resolved = ALL
            elif character == '-':
                resolved = ['-']
            else:
                try:
                    resolved = ambiguities[character]
                except KeyError:
                    raise AlphabetError(ambig_motif)
            for character2 in resolved:
                for motif in motif_set:
                    new_motifs.append(''.join([motif, character2]))
            
            motif_set = new_motifs
        
        # delete sub motifs that are not to be included
        motif_set = [motif for motif in motif_set if motif in self._quick_motifset]
        
        if not motif_set:
            raise AlphabetError(ambig_motif)
        
        return tuple(motif_set)
    
    def adaptMotifProbs(self, motif_probs):
        """Prepare an array or dictionary of probabilities for use with
        this alphabet by checking size and order"""
        if hasattr(motif_probs, 'keys'):
            sample = motif_probs.keys()[0]
            if sample not in self:
                raise ValueError("Can't find motif %s in alphabet" %
                                sample)
            motif_probs = numpy.array(
                    [motif_probs[motif] for motif in self])
        else:
            if len(motif_probs) != len(self):
                if len(motif_probs) != len(self):
                    raise ValueError("Can't match %s probs to %s alphabet" %
                            (len(motif_probs), len(self)))
            motif_probs = numpy.asarray(motif_probs)
        assert abs(sum(motif_probs)-1.0) < 0.0001, motif_probs
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
    def __init__(self, data=[], Gap='-', MolType=None):
        """Initializes self from items.
        
        data should be a sequence (string, list, etc.) of characters that
        are in the alphabet, e.g. 'UCAG' for RNA.
        
        Gap should be a single character that represents the gap, e.g. '-'.
        """
        super(CharAlphabet, self).__init__(data, Gap, MolType=MolType)
        self._indices_to_chars, self._chars_to_indices = \
            _make_translation_tables(data)
        self._char_nums_to_indices = array(self._chars_to_indices,'c').view('B')
        self._indices_nums_to_chars = array(self._indices_to_chars, 'c')
    
    def fromString(self, data):
        """Returns array of indices from string containing elements.
        
        data should be a string on the alphabet, e.g. 'ACC' for the RNA
        alhabet 'UCAG' would return the array [2,1,1]. This is useful for
        converting strings into arrays of small integers on the alphabet,
        e.g. for reading a Sequence from a string.
        
        This is on the Alphabet, not the Sequence, because lots of objects
        (e.g. Profile, Alignment) also need to use it.
        """
        return fromstring(translate(data, self._chars_to_indices), uint8)
    
    def isValid(self, seq):
        """Returns True if seq contains only items in self."""
        try:
            if len(seq) == 0:   #can't be invalid if empty
                return True
            ind = self.toIndices(seq)
            return max(ind) < len(self) and min(ind) >= 0
        except (TypeError, KeyError):
            return False
    
    def fromArray(self, data):
        """Returns array of indices from array containing elements.
        
        This is useful if, instead of a string, you have an array of
        characters that's been converted into a numpy array. See
        fromString docstring for general behavior.
        """
        return take(self._char_nums_to_indices, data.view('B'))
    
    def toChars(self, data):
        """Converts array of indices into array of elements.
        
        For example, on the 'UCAG' RNA alphabet, an array with the data
        [0,1,1] would return the characters [U,C,C] in a byte array.
        """
        return take(self._indices_nums_to_chars, data.astype('B'))
    
    def toString(self, data, delimiter='\n'):
        """Converts array of data into string.
        
        For example, on the 'UCAG' RNA alphabet, an array with the data
        [0,1,1] would return the string 'UCC'. This is the most useful
        conversion mechanism, and is used by e.g. the Sequence object to
        convert the internal data back into strings for output.
        """
        s = data.shape
        if not s:
            return ''
        elif len(s) == 1:
            return self.toChars(data).tostring()
        else:
            return delimiter.join([i.tostring() for i in self.toChars(data)])
    

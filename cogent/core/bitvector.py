#!/usr/bin/env python
"""Classes for dealing with bitvectors: arrays of 1 and 0.

bitvectors are often much faster than numpy arrays for bit operations,
especially if several characters (e.g. of DNA) can be packed into a
single byte.

Usage: v = Bitvector('11000101')

Provides the following major classes and factory functions:

    ShortBitvector and LongBitvector: subclass ImmutableBitvector and are
    produced by the factory function Bitvector(). Short uses an int, while
    Long uses a long. ShortBitvectors in particular are _very_ fast, but are
    limited to 31 bits. Bitvector() will return the appropriate kind of 
    vector for your sequence.

    MutableBitvector: always a LongBitvector. Can be changed through 
    __getitem__, e.g. vec[3] = 0. Use freeze() and thaw() methods to
    convert between mutable and immutable Bitvectors (both types define both
    methods).

    VectorFromCases: constructs a Bitvector from letters in a string.

    VectorFromMatches: constructs a Bitvector from a text and a pattern, with
    1 at the positions where there's a match. Pattern can be a string or a
    regex.

    PackedBases: subclasses LongBitvector: provides a way of storing nucleic
    acid bases compactly and in a format convenient for assessing sequence
    similarity.

    Note: there isn't a general VectorFromX factory function because if you
    have a function, you can always do Bitvector(map(func, items)), making
    VectorFromX trivial.
"""

from cogent.util.misc import Delegator
import re
from string import maketrans
from operator import and_, or_, xor
from numpy import log2
from sys import maxint

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

_bits_in_int = int(round(log2(maxint)))


def is_nonzero_string_char(char):
    """Tests whether input character is not '0', returning '1' or '0'."""
    if char == '0' or char == '':
        return '0'
    else:
        return '1'

def is_nonzero_char(item):
    """Tests whether input item is nonzero, returning '1' or '0'."""
    if isinstance(item, str):
        return is_nonzero_string_char(item)
    else:
        if item:
            return '1'
        else:
            return '0'

def seq_to_bitstring(seq):
    """Converts sequence to string of 1 and 0, returning string of '1'/'0'."""
    if not seq:
        return ''
    if isinstance(seq, str):
        return ''.join(map(is_nonzero_string_char, seq))
    else:
        return ''.join(map(is_nonzero_char, seq))
        
def is_nonzero_string_int(char):
    """Tests whether input character is not '0', returning 1 or 0."""
    if char == '0' or char == '':
        return 0
    else:
        return 1

def is_nonzero_int(item):
    """Tests whether input item is nonzero, returning 1 or 0."""
    if isinstance(item, str):
        return is_nonzero_string_int(item)
    else:
        if item:
            return 1
        else:
            return 0

def seq_to_bitlist(seq):
    """Converts sequence to list of 1 and 0."""
    if not seq:
        return []
    if isinstance(seq, str):
        return map(is_nonzero_string_int, seq)
    else:
        return map(is_nonzero_int, seq)

def num_to_bitstring(num, length):
    """Returns string of bits from number, truncated at length.

    Algorithm: While the number storing the bitvector is greater than
    zero, find the last digit as (num mod 2) and then bit-shift the string to 
    the right to delete the last digit. Add the digits into a string from 
    right to left (in reverse order).

    In other words, when length is too small for the number, takes the last
    n bits of the number in order.

    Warning: this is not terribly efficient, since the entire number
    must be rewritten at each shift step. In mutable bitvectors,
    use caching to reduce the number of conversions required.
               
    """
    bits = [0] * length #start with all zeroes
    #successively find the last bit of curr, replacing the
    #appropriate element of bits as each bit is processed and
    #removed
    while num and length:
        bits[length-1] = 1 & num
        num >>= 1
        length -= 1
    return ''.join(map(str, bits))
    
def bitcount(num, length, value=1):
    """Counts the bits in num that are set to value (1 or 0, default 1)."""
    one_count = 0
    curr_length = length
    while num and curr_length:
        one_count += 1 & num
        num >>= 1
        curr_length -= 1
    if value:
        return one_count
    else:
        return length - one_count

class ImmutableBitvector(object):
    """Generic interface for immutable bitvectors."""

    def __init__(self, items='', length=None):
        """Returns new Bitvector."""
        if length is not None:
            self._length = length
        else:
            self._length = len(items)

    def __str__(self):
        """Returns string representation of bitvector."""
        return num_to_bitstring(self, self._length)

    def __len__(self):
        """Returns length of the bitvector."""
        return self._length

    def op(self, other, func):
        """Performs bitwise op on self and other, returning new Bitvector."""
        self_length = self._length
        if not isinstance(other, ImmutableBitvector):
            #coerce to comparable type
            try:
                other = other.freeze()
                other_length = other._length
            except:
                other = Bitvector(other, self_length)
                other_length = self_length
        else:
            other_length = other._length
            
        if self_length == other_length:
        #if they're the same length, just combine the vectors
            return Bitvector(func(self, other), self_length)
        else:
        #need to find which is larger, right-shift, and combine
            diff = self_length - other_length
            if diff < 0:    #re-bind self and other so self is now longest
                #tuple unpacking swaps self and other
                self_length, other_length = other_length, self_length
                self, other = other, self
            #shift the shorter vector and combine with bitwise function
            return Bitvector(func(self >> abs(diff), other), other_length)

    def __getitem__(self, item):
        """Returns self[item] as 1 or 0, i.e. as integers and not strings.

        key should be an index less than the length of the bitvector; negative
        indexes are handled in the usual manner.

        Algorithm:  uses bit masks to figure out whether the specific
                    position is occupied. Caches masks in Bitvector._masks, so
                    should be very fast when only a few different lengths of
                    bitvector are used in a program. Will be somewhat
                    inefficient if many large bitvectors of different sizes
                    are each used occasionally.
        """
        #check that there actually are items in the vector
        length = self._length
        if not length:
            raise IndexError, "Can't find items in empty vector!"
        
        if isinstance(item, slice):
            return Bitvector(''.join(map(str,[self[i] \
                for i in range(*item.indices(length))])))
        else:
            #transform keys to allow negative index
            if item < 0:
                item += length
            #check that key is in bounds
            if (item < 0) or (item >= length):
                raise IndexError, "Index %s out of range." % (item,)
            move_to = length - item - 1
            if move_to >= _bits_in_int:
                result =  self & (1L << move_to)
            else:
                result = self & (1 << move_to)
            if result:
                return 1
            else:
                return 0

    def bitcount(self, value=1):
        """Counts the bits in self that match value."""
        return bitcount(self, self._length, value)

    def __repr__(self):
        """Produces standard object representation instead of int."""
        c = self.__class__
        return "<%s.%s object at %s>" % (c.__module__,c.__name__,hex(id(self)))

    def freeze(self):
        """Returns ImmutableBitvector containing data from self."""
        return self

    def thaw(self):
        """Returns MutableBitvector containing data from self."""
        return MutableBitvector(self)

    def stateChanges(self):
        """Returns list of indices where state changes from 0->1 or 1->0."""
        #bail out if no elements rather than raising IndexError
        length = len(self)
        if not length:
            return []
        bits = list(self)
        changes = [i for i in range(1, length) if self[i] != self[i-1]]
        return [0] + changes + [length] #always implicitly add start and end

    def divideSequence(self, seq, state=None):
        """Divides sequence into a list of subsequences, cut using stateChanges.

        The list will contain the slices for the indices containing the given
        state.  If no state given (state = None) the list will consist of all
        of the subsequences, cut at the indices.
        
        The list and self[0] are returned as a tuple.

        Truncates whichever sequence is shorter.
        """
        #bail out rather than raising IndexError if vector or sequence is empty
        if not (len(self) and seq):
            return ([], 0)

        cut_list = self.stateChanges()
        cut_seq = []
        first = 0
        cut_index = len(cut_list)-1
        #If user supplied a specific state, return only pieces in that state
        if state is not None:
            #test whether we want to include the first segment or not
            exclude_start = self[0] != state
            for i in range(cut_index):
                if i % 2 == exclude_start:
                    cut_seq.append(seq[cut_list[i]:cut_list[i+1]])
        else:
            #No state supplied: return pieces in all states
            for i in range(cut_index):
                cut_seq.append(seq[cut_list[i]:cut_list[i+1]])
        first = self[0]

        return filter(None, cut_seq), first

class ShortBitvector(ImmutableBitvector, int):
    """Short bitvector, stored as an int."""
    __slots__ = ['_length']
    
    def __new__(cls, items='', length='ignored'):
        """Creates new bitvector from sequence of items."""
        if isinstance(items, str):
            if items:   #guard against empty string
                return int.__new__(cls, items, 2)
            else:
                return int.__new__(cls, 0)
        else:
            return int.__new__(cls, items)

    def __or__(self, other):
        """Returns position-wise OR (true if either true)."""
        if isinstance(other, long):
            return LongBitvector(self).__or__(other)
        else:
            return self.op(other, int.__or__)

    def __and__(self, other):
        """Returns position-wise AND (true if both true)."""
        if isinstance(other, long):
            return LongBitvector(self).__and__(other)
        else:
            return self.op(other, int.__and__)

    def __xor__(self, other):
        """Returns position-wise XOR of self and other (true if states differ).
        """
        if isinstance(other, long):
            return LongBitvector(self).__xor__(other)
        else:
            return self.op(other, int.__xor__)

    def __invert__(self):
        """Returns complement (replace 1 with 0 and vice versa)."""
        length = self._length
        if length >= _bits_in_int:
            mask = (1L << length) - 1
        else:
            mask = (1 << length) - 1
        return Bitvector(mask & ~int(self), length)

class LongBitvector(ImmutableBitvector, long):
    """Long bitvector, stored as a long."""
    
    def __new__(cls, items='', length='ignored'):
        """Creates new bitvector from sequence of items."""
        if isinstance(items, str):
            if items:   #guard against empty string
                return long.__new__(cls, items, 2)
            else:
                return long.__new__(cls, 0)
        else:
            return long.__new__(cls, items)

    def __or__(self, other):
        """Returns position-wise OR (true if either true)."""
        return self.op(other, long.__or__)

    def __and__(self, other):
        """Returns position-wise AND (true if both true)."""
        return self.op(other, long.__and__)

    def __xor__(self, other):
        """Returns position-wise XOR of self and other (true if states differ).
        """
        return self.op(other, long.__xor__)

    def __invert__(self):
        """Returns complement (replace 1 with 0 and vice versa)."""
        length = self._length
        return Bitvector(((1L << length) - 1) & ~long(self), length)

def Bitvector(items='', length=None, constructor=None):
    """Factory function returning short or long Bitvector depending on length.
    
    Note: uses explict test rather than try/except to fix memory leak.

    Now is the only way to convert an arbitrary sequence of true/false values
    into a bitvector.
    """
    #convert whatever was passed as 'items' into a number
    if isinstance(items, ImmutableBitvector):
        num = items
        if length is None:
            length = len(items)
    elif isinstance(items, int) or isinstance(items, long):
        num = items
        if length is None:
            raise TypeError, "Must specify length if initializing with number."
    else:
        bitstring = seq_to_bitstring(items)
        if bitstring:
            num = long(bitstring, 2)
            if length is None:
                length = len(items)
        else:
            num = 0
            length = 0
    #if the constructor was not passed explicitly, guess it
    if not constructor:
        if 0 <= num <= maxint:
            constructor = ShortBitvector
        else:
            constructor = LongBitvector
                
    return constructor(num, length)

# Original version with memory leak under Python 2.3
#    try:
#        return ShortBitvector(items, length)
#    except (OverflowError, TypeError):
#        return LongBitvector(items, length)


class MutableBitvector(Delegator):
    """Array of bits (0 or 1) supporting set operations. Supports __setitem__.
    """
    def __init__(self, items='', length=None):
        """Initializes the Bitvector class, storing data in long _vector.

        Items can be any sequence with members that evaluate to True or False.

        Private data:
            _vector:     long representing the bitvector
            _string:     string representation of the vector
            _is_current: boolean indicating whether the string is known to
                         match the long data
        """
        Delegator.__init__(self, None)
        self.replace(items, length)

    def replace(self, items='', length=None):
        """Replaces the contents of self with data in items.

        Items should be a sequence with elements that evaluate to True or False.

        If items is another BitVector, uses an efficient method based on the
        internal data. Otherwise, will loop through items evaluating each to
        True or False.

        Primarily used internally, but exposed publically because it's often
        convenient to replace the contents of an existing bitvector.

        Usage:  vec.update(foo) where foo is a sequence or bitvector.

        """
        vec = Bitvector(items, length)
        self._handler = vec     #make sure attributes go to the new vec
        self._is_current = False
        self._string = ''

    def __str__(self):
        """Prints the bits in the vector as a string of '1' and '0'."""
        #check that the string is up to date; when it is, return it.
        if not self._is_current:
            self._string = str(self._handler)
            self._is_current = True
        return self._string

    def __len__(self):
        """Need to pass directly to handler."""
        return len(self._handler)

    def __getitem__(self, *args):
        """Gets item from vector: needs to be passed through explicitly."""
        return self._handler.__getitem__(*args)

    def __setitem__(self, item, value):
        """Sets self[item] to value: value must be 0 or 1, or '0' or '1'.

        Uses masks to alter the bit at a specific position. To set a bit
        to zero, make a string that is 1 everywhere except at that position
        and then use bitwise &. To set a bit to 1, make a string that is 0
        everywhere except at that position and then use bitwise |.

        Warning: this is quite a slow operation, especially on long vectors.
        """
        length = self._length
        #handle slice assignment. Warning: does not allow change in length!
        if isinstance(item, slice):
            for i, val in zip(range(*item.indices(length)), value):
                self[i] = is_nonzero_int(val)
            return
        
        #otherwise, check that the key is in range
        vec = self._handler
        curr_val = vec[item]
        val = is_nonzero_int(value)

        if val != curr_val: #need only do anything if value changed!
            #find index of negative values
            if item < 0:
                item += length
            #check that key is in bounds
            if (item < 0) or (item >= length):
                raise IndexError, "Index %s out of range." % (item,)
            
            #figure out offset
            move_to = length - item - 1
            if value == 0:
                #make mask of '1', and combine with &
                result = vec &((1L << length) - 1) ^ (1L << move_to)

            elif value == 1:
                #make mask of '0', and combine with |
                result = vec | (1L << move_to)

            else:
                #complain if we got anything else
                raise ValueError, "Item not 0 or 1: " + str(value)
            
            self.replace(result)

    def __cmp__(self, other):
        """Comparison needs to be passed through explicitly."""
        return cmp(self._handler, other)

    def op(self, other, func):
        """Returns position-wise op for self and other.

        self[i] == 1 if op(self[i], other[i]) == 1

        Internal method: used to implement bitwise operations on pairs of
        bitvectors. Truncates the longer vector to fit the shorter one by
        right-shifting, i.e. the leftmost bits of the longer vector will
        be combined with all the bits of the shorter vector.

        Handles the case where both the operands are bitvectors specially for
        speed, to greatest advantage when they are equal length.

        Uses bool() to convert elements of sequences, in preference to checking
        whether they directly represent 0 or 1. This makes the operations more
        flexible, but beware of bad data!
        """
        return self.__class__(self._handler.op(other, func))

    def __or__(self, other):
        """Returns position-wise OR (true if either true).
        """
        return self.op(other, or_)

    def __and__(self, other):
        """Returns position-wise AND (true if both true).
        """
        return self.op(other, and_)

    def __xor__(self, other):
        """Returns position-wise XOR of self and other (true if states differ).
        """
        return self.op(other, xor)

    def __invert__(self):
        """Returns complement (replace 1 with 0 and vice versa).
        """
        return self.__class__(~self._handler)

    def thaw(self):
        """Returns mutable version of self."""
        return self

    def freeze(self):
        """Returns immutable version of self."""
        return self._handler

def VectorFromCases(seq, constructor=LongBitvector):
    """Returns vector with 0(1) for lower(upper)case positions in sequence.

    Primarily used for identifying binding sites denoted by uppercase bases.
    """
    return constructor(''.join(map(str, map(int, [i.isupper() for i in seq]))))

def VectorFromMatches(seq, pattern, overlapping=1, constructor=LongBitvector):
    """Replaces self with 1(0) for each position in sequence (not) matching.

    Usage:  vec.toMatches('agagac', 'aga', 1) == vec('111110')
            vec.toMatches('agagac', 'aga', 0) == vec('111000')

    The first argument must be a type that is searchable by a re object,
    typically a string. The second argument must be a pattern to search:
    this can be a string, a list of strings, or an re object. The third,
    optional argument specifies whether matches may overlap or not.

    If a list is passed, its elements will be converted into strings, and
    then into a regular expression using alternation. Lists of regular
    expressions are not supported directly.

    Zero-length patterns will typically never be counted as matching,
    rather than being counted as matching everywhere.
    """
    #allocate space for the results
    seqlength = len(seq)
    result = [0] * seqlength

    if seqlength:   #will return empty string if seq empty
        if isinstance(pattern, str):
            #handle strings by using the index string method
            patlength = len(pattern)
            #will skip zero-length patterns
            if patlength:
                #until we reach the end of the sequence, find the next
                #index that contains a match. When a match has been
                #processed by inserting '1' into the list at each position
                #that matches, either move to the next position (if matches
                #can overlap) or to the position after the end of the match
                #(if matches cannot overlap).
                start = 0
                while start + patlength - 1 <= seqlength:
                    try:
                        match = seq.index(pattern, start)
                    except ValueError:
                        break   #no more matches
                    #insert '1' into the list in the appropriate slice
                    result[match:match + patlength] = [1] * patlength
                    #move to the first index where the next match could be
                    if overlapping:
                        start = match + 1
                    else:
                        start = match + patlength
        else:
            #if not a string, assume regex
            #first, attempt to alternate the items in pattern as though it
            #were a sequence type
            try:
                pattern = re.compile('|'.join(map(str, pattern)))
            except:
                pass
            #use same strategy as above to locate matches of pattern, find
            #the lengths, insert '1' into the list at the appropriate
            #places, and move on to the next position where a match might
            #start.
            start = 0
            while start < seqlength:
                match = pattern.search(seq, start)
                try:
                    #find the length of the match, and the start and end
                    #indexes. Use these to define the positions where '1'
                    #will appear in the result.
                    match_start, match_end = match.span()
                    match_length = match_end - match_start
                    result[match_start:match_end] = [1] * match_length
                    #figure out how many positions to advance for the next
                    #possible match start
                    if overlapping or (match_length == 0):
                        start = match_start + 1
                    else:
                        start = match_end
                except AttributeError, TypeError:
                    #match is None if no more in string
                    break
    else:   #sequence was zero-length
        result = []
    return constructor(''.join(map(str, result)))

def VectorFromRuns(seq, length, constructor=LongBitvector):
    """Returns vector with 1 at positions specified by sequence of (start,len).

    seq should be a sequence of 2-item sequences, where the first element is
    the index where the run of 1's starts and the second element is the number
    of 1's in the run.

    length should be an int giving the total length of the sequence.
    """
    if not length:
        return constructor('')
    
    bits = ['0']*length
    for start, run in seq:
        if start + run > length:
            raise IndexError, "start %s + run %s exceeds seq length %s" % \
                (start, run, length)
        bits[start:start+run] = ['1']*run
    return constructor(''.join(bits))

def VectorFromSpans(seq, length, constructor=LongBitvector):
    """Returns vector with 1 at positions specified by sequence of (start,end).

    seq should be a sequence of 2-item sequences, where the first element is
    the index where the run of 1's starts and the second element is the index
    where the run stops (standard Python slice notation, i.e. the start and
    stop you would normally use in a slice, where the slice does not include
    the element at i[stop]).

    length should be an int giving the total length of the sequence.
    """
    if not length:
        return constructor('')
    
    bits = ['0']*length
    for start, stop in seq:
        if stop > length:
            raise IndexError, "stop %s exceeds seq length %s" % (stop, length)
        run = stop - start
        bits[start:stop] = ['1']*run
    
    return constructor(''.join(bits))

def VectorFromPositions(seq, length, constructor=LongBitvector):
    """Returns vector with 1 at positions specified by sequence of positions.

    seq should be a sequence of ints (position).
    
    length should be an int giving the total length of the sequence.
    """
    if not length:
        return constructor('')
    
    bits = ['0']*length
    for i in seq:
        bits[i] = '1'
    return constructor(''.join(bits))



class PackedBases(LongBitvector):
    """Stores unambiguous nucleotide sequences as bitvectors, 4 per byte.

    Rna controls whether __str__ produces RNA or DNA.

    Each base is stored as 2 bits. The first bit encodes purine/pyrmidine,
    while the second encodes weak/strong.

    NOTE: len() returns the length in _bits_, not the length in bases.
    """
    #translation table to turn bases into sequences. First bit encodes purine/
    #pyrimindine; second bit encodes weak/strong.
    _sequence_to_bits = maketrans('AGUTCagutc', '0122301223')
    _rna_bases = {0:'A', 1:'G', 2:'U', 3:'C'}
    _dna_bases = {0:'A', 1:'G', 2:'T', 3:'C'}

    def __new__(cls, sequence='', Rna='ignored'):
        """Packs bases in sequence into a 2-bit per base bitvector.

        Usage:  vec.toBases('AcgcaAc') == vec('00110111000011')

        Encoding is to use the first bit as purine vs pyrimidine (purine = 0),
        and to use the second bit as weak vs. strong (weak = 0). This means
        that A == 00, G == 01, U == 10, and C == 11.
        """
        seq = sequence.translate(cls._sequence_to_bits)
        if not seq:
            return long.__new__(cls, 0L)
        else:
            return long.__new__(cls, long(seq, 4))  #note base 4 conversion

    def __init__(self, sequence='', Rna=True):
        """Returns new PackedBases object."""
        self._length = len(sequence)*2
        self.Rna = Rna

    def __str__(self):
        """Unpacks bases from sequence using encoding scheme in toBases."""
        length = self._length/2
        bits = [0] * length #allocate space using '00' pattern
        curr = self
        #successively find the last two bit of curr, replacing the
        #appropriate element of bits as each bit is processed and
        #removed
        while curr and length:
            bits[length-1] = curr & 3
            curr >>= 2
            length -= 1
        #convert numbers back into string
        if self.Rna:
            num2base = self._rna_bases
        else:
            num2base = self._dna_bases
        return ''.join([num2base[i] for i in bits])


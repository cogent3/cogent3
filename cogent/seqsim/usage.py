#!/usr/bin/env python
"""usage.py: usage of symbols, including substitutions on pairwise alphabets.

Revision History

Created 10/12/04 by Rob Knight.

9/14/05 Rob Knight: Changed Usage constructor to allow Alphabet on the instance
level, and to eliminate the precalculated flag which was not used. Added
entropy method.

7/20/07 Mike Robeson: Under PairMatrix.__init__ changed 'if data:' to
        'if data != None:
8/3/07 Daniel McDonald: Code now relies on numpy and cogent with the exception
        of the one scipy function that still needs to be removed
"""
from cogent.maths.scipy_optimize import fmin, brent
from cogent.util.array import scale_trace, norm_diff, \
    has_neg_off_diags, sum_neg_off_diags, with_diag, without_diag

from cogent.core.alphabet import get_array_type
from cogent.core.usage import RnaBases, DnaBases, DnaPairs, RnaPairs, Codons
from cogent.core.sequence import ModelSequence, ModelDnaSequence, \
    ModelRnaSequence
from operator import add, sub, mul, div
from cogent.maths.matrix_logarithm import logm
from cogent.maths.stats.util import FreqsI
from cogent.maths.matrix_exponentiation import FastExponentiator as expm
from numpy import zeros, array, max, diag, log, nonzero, product, cumsum, \
                  searchsorted, exp, diagonal, choose, less, repeat, average,\
                  logical_and, logical_or, logical_not, transpose, compress,\
                  ravel, concatenate, equal, log, dot, identity, \
                  newaxis as NewAxis, sum, take, reshape, any, all, asarray
from numpy.linalg import eig
from numpy.linalg import inv as inverse
from numpy.random import random as randarray

ARRAY_TYPE = type(array([0]))

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Mike Robeson", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class Usage(FreqsI):
    """Stores usage on a particular alphabet. Abstract class.
    
    Note: Usage is abstract because most subclasses (e.g. CodonUsage,
    AminoAcidUsage) have specific methods that depend on their alphabets.
    Allowing generic Usage objects is disallowed to enforce use of the
    appropriate Usage object for specific situations.

    Supports most of the Cogent FreqsI interface.
    """
    Alphabet = None # concrete subclasses have specific alphabets
    
    def __init__(self, data=None, Alphabet=None):
        """Returns a new Usage object from array of symbol freqs.
        
        Will interpret many different kinds of data, including precalculated
        frequencies, arrays of symbols, and cogent.core.sequence.ModelSequence 
        objects.

        Warning: it guesses whether you passed in frequencies or symbols based
        on the length of the array, so for example Usage(DnaSequence('ATCG'))
        will _not_ give the result you expect. If you know the data type,
        use the alternative class method constructors.
        """
        if Alphabet is not None:
            self.Alphabet = Alphabet
        if not self.Alphabet:
            raise TypeError, "Usage subclasses must define alphabet."""
        if isinstance(data, Usage):
            self._data = data._data
        else:
            self._data = zeros(len(self), 'float64')
            if any(data):
                self += data

    def __getitem__(self, i):
        """Returns item based on alphabet."""
        return self._data[self.Alphabet.index(i)]

    def __setitem__(self, key, val):
        """Sets item based on alphabet."""
        self._data[self.Alphabet.index(key)] = val

    def __str__(self):
        """Prints as though it were a tuple of key,value pairs."""
        return str(self.items())

    def __repr__(self):
        """String representation of self."""
        return ''.join([self.__class__.__name__, '(', repr(self._data), ')'])

    def __iter__(self):
        """Iterates over keys, like a dict."""
        return iter(self.Alphabet)

    def __eq__(self, other):
        """Tests whether two Usage objects have the same data."""
        if hasattr(other, '_data'):
            return all(self._data == other._data)
        #if we get here, didn't compare equal
        try:
            return all(self._data == self.__class__(other)._data)
        except:
            return False

    def __ne__(self, other):
        """Returns True if self and other are not equal."""
        if hasattr(other, '_data'):
            return any(self._data != other._data)
        #if we get here, didn't compare equal
        try:
            return any(self._data != self.__class__(other)._data)
        except:
            return True

    def __iadd__(self, other):
        """Adds data to self in-place."""
        #check if other is nonzero; skip if it isn't
        try:
            if not other:
                return self
        except ValueError:
            if not any(other):
                return self
        #first, check if it's a Usage object
        if isinstance(other, Usage):
            self._data += other._data
            return self
        #then, check if it's one of our ModelSequence objects
        ac = self.Alphabet.counts
        if isinstance(other, ModelSequence):
            self._data += ac(other._data)
            return self
        #if it's the same length as self, try to add it as frequencies
        try:
            if len(other) == len(self):
                self._data += other
                return self
        except TypeError:
            pass
        #then try to convert it using the alphabet
        #WARNING: this will silently ignore unknown keys!
        #since we know other wasn't nonzero, we won't accept
        #the result if we can't convert anything.
        try:
            other_freqs = ac(other)
            #check if we actually converted anything...
            if any(other_freqs):
                self._data += other_freqs
                return self
        except (IndexError, KeyError, TypeError):
            pass
        #then use the generic conversion function
        f = self._find_conversion_function(other)
        if f:
            f(other, op=add)
            return self
        else:
            raise TypeError, "Could not convert this to freqs: %s" % other

    def __isub__(self, other):
        """Subtracts data from self in-place."""
        #check if other is nonzero; skip if it isn't
        try:
            if not other:
                return self
        except ValueError:
            if not any(other):
                return self
        #first, check if it's a Usage object
        if isinstance(other, Usage):
            self._data -= other._data
            return self
        #then, check if it's one of our ModelSequence objects
        ac = self.Alphabet.counts
        if isinstance(other, ModelSequence):
            self._data -= ac(other._data)
            return self
        #if it's the same length as self, try to add it as frequencies
        try:
            if len(other) == len(self):
                self._data -= other
                return self
        except TypeError:
            pass
        #then try to convert it using the alphabet
        #WARNING: this will silently ignore unknown keys!
        #since we know other wasn't nonzero, we won't accept
        #the result if we can't convert anything.
        try:
            other_freqs = ac(other)
            #check if we actually converted anything...
            if other_freqs.any():
                self._data -= other_freqs
                return self
        except (IndexError, KeyError, TypeError):
            pass
        #then use the generic conversion function
        f = self._find_conversion_function(other)
        if f:
            f(other, op=sub)
            return self
        else:
            raise TypeError, "Could not convert this to freqs: %s" % other
  
    def __mul__(self, other):
        """Multiplies self by other (assumed scalar)."""
        return self.__class__(self._data * other)

    def __imul__(self, other):
        """Multiplies self by other in-place (assumed scalar)."""
        self._data *= other
        
    def __div__(self, other):
        """Divides self by other (assumed scalar). Always true division."""
        return self.__class__(self._data / (other))

    def __idiv__(self, other):
        """Divides self by other (assumed scalar) inplace. Maybe int division."""
        self._data /= other

    def scale_sum(self, sum_=1.0):
        """Returns copy of self scaled to specified sum."""
        return self.__class__(self._data * (sum_/sum(self._data)))

    def scale_max(self, max_=1.0):
        """Returns copy of self scaled to specified maximum (default 1)."""
        return self.__class__(self._data * (max_/max(self._data)))

    def probs(self):
        """Returns copy of self scaled so that the sum is 1."""
        return self.__class__(self._data / (sum(self._data)))

    def randomIndices(self, length, random_vector=None):
        """Produces random indices according to symbol freqs."""
        freqs = cumsum(self._data/sum(self._data))[:-1]
        if random_vector is None:
            random_vector=randarray(length)
        return searchsorted(freqs, random_vector)

    def fromSeqData(cls, seq, Alphabet=None):
        """Returns new Usage object from Sequence object."""
        return cls.fromArray(seq._data, Alphabet=Alphabet)
    
    def fromArray(cls, a, Alphabet=None):
        """Returns new Usage object from array."""
        return cls(cls.Alphabet.counts(a), Alphabet=Alphabet)

    fromSeqData = classmethod(fromSeqData)
    fromArray = classmethod(fromArray)
    
    #following code is to support FreqsI
    def get(self, key, default):
        """Returns self._data[self.Alphabet.index(key) if present, or default."""
        try:
            return self._data[self.Alphabet.index(key)]
        except (KeyError, IndexError, TypeError):
            return default

    def values(self):
        """Returns list of keys in self (i.e. the alphabet)."""
        return list(self._data)

    def keys(self):
        """Returns list of values in self (i.e. the data)."""
        return list(self.Alphabet)

    def items(self):
        """Returns list of (key, value) pairs in self."""
        return zip(self.Alphabet, self._data)

    def isValid(self):
        """Always valid (except for negative numbers), so override."""
        return min(self._data) >= 0

    def copy(self):
        """Return copy of self with same alphabet, not sharing data."""
        return self.__class__(self._data.copy())

    def __delitem__(self, key):
        """Can't really delete items, but raise error if in alphabet."""
        if key in self.Alphabet:
            raise KeyError, "May not delete required key %s" % key

    def purge(self):
        """Can't contain anything not in alphabet, so do nothing."""
        pass

    def normalize(self, total=1.0, purge=True):
        """Converts counts into probabilities, normalized to 1 in-place.
        
        Changes result to Float64. Purge is always treated as True. 
        """
        if self._data is not None and self._data.any():
            self._data = self._data / (total * sum(self._data))

    def choice(self, prob):
        """Returns item corresponding to Pr(prob)."""
        if prob > 1:
            return self.Alphabet[-1]
        summed = cumsum(self._data/sum(self._data))
        return self.Alphabet[searchsorted(summed, prob)]

    def randomSequence(self, n):
        """Returns list of n random choices, with replacement."""
        if not self:
            raise IndexError, "All frequencies are zero."
        return list(choose(self.randomIndices(n), self.Alphabet))

    def subset(self, items, keep=True):
        """Sets all frequencies not in items to 0.

        If keep is False, sets all frequencies in items to 0.
        """
        if keep:
            for i in self.Alphabet:
                if i not in items:
                    self[i] = 0
        else:
            for i in items:
                try:
                    self[i] = 0
                except KeyError:
                    pass

    def scale(self, factor=1, offset=0):
        """Linear transform of values in freqs where val= factor*val + offset."""
        self._data = factor * self._data + offset

    def __len__(self):
        """Returns length of alphabet."""
        return len(self.Alphabet)

    def setdefault(self, key, default):
        """Returns self[key] or sets self[key] to default."""
        if self[key]:
            return self[key]
        else:
            self[key] = default
            return default
        

    def __contains__(self, key):
        """Returns True if key in self."""
        try:
            return key in self.Alphabet
        except TypeError:
            return False
    
    def __nonzero__(self):
        """Returns True if self is nonzero."""
        return bool(sum(self._data) != 0)
       
    def rekey(self, key_map, default=None, constructor=None):
        """Returns new Freqs with keys remapped using key_map.

        key_map should be a dict of {old_key:new_key}.
        
        Values are summed across all keys that map to the same new value.
        Keys that are not in the key_map are omitted (if default is None),
        or set to the default.

        constructor defaults to self.__class__. However, if you're doing
        something like mapping amino acid frequencies onto charge frequencies,
        you probably want to specify the constructor since the result won't
        be valid on the alphabet of the current class.

        Note that the resulting Freqs object is not required to contain
        values for all the possible keys.
        """
        if constructor is None:
            constructor = self.__class__
        result = constructor()
        for key, val in self.items():
            new_key = key_map.get(key, default)
            curr = result.get(new_key, 0)
            try:
                result[new_key] = curr + val
            except KeyError:
                pass
        return result 

    def entropy(self, base=2):
        """Returns Shannon entropy of usage: sum of p log p."""
        ln_base = log(base)
        flat = ravel(self._data)
        total = sum(flat)
        if not total:
            return 0
        flat /= total
        ok_indices = nonzero(flat)[0]
        ok_vals = take(flat, ok_indices, axis=0)
        return -sum(ok_vals * log(ok_vals))/ln_base
   
class DnaUsage(Usage):
    """Stores usage on the DNA alphabet."""
    Alphabet = DnaBases
        
class RnaUsage(Usage):
    """Stores usage on the RNA alphabet."""
    Alphabet = RnaBases

class CodonUsage(Usage):
    """Stores usage on the Codon alphabet."""
    Alphabet = Codons

class DnaPairUsage(Usage):
    """Stores usage on the DnaPairs alphabet."""
    Alphabet = DnaPairs

class RnaPairUsage(Usage):
    """Stores usage on the RnaPairs alphabet."""
    Alphabet = RnaPairs

class PairMatrix(object):
    """Base class for Counts, Probs, and Rates matrices. Immutable.
    
    Holds any numeric relationship between pairs of objects on a JointAlphabet.
    Note that the two SubEnumerations of the JointAlphabet need not be the same,
    although many subclasses of PairMatrix will require that the two
    SubEnumerations _are_ the same because their methods assume square matrices.
    """
    def __init__(self, data, Alphabet, Name=None):
        """Returns new PairMatrix object containing data.
        
        WARNING: Alphabet must be a JointAlphabet where the two SubEnumerations
        are the same.
        """
        self.Alphabet = Alphabet
        if any(data):
            self._data = reshape(array(data, 'd'), Alphabet.Shape)
        else:
            self._data = zeros(Alphabet.Shape, 'd')
        self.Name = Name
 
    def toMatlab(self):
        """Returns Matlab-formatted string representation."""
        if self.Name is None:
            name = 'm'
        else:
            name = str(self.Name)
        return ''.join([name, '=', '[', \
            ';\n'.join([' '.join(map(str, r)) for r in self._data]), '];\n'])

    def __str__(self):
        """Returns string representation of array held in self."""
        return str(self._data)

    def __repr__(self):
        """Returns string representation of self."""
        return ''.join([self.__class__.__name__, '(', repr(self._data), \
            ',', repr(self.Alphabet), ',', repr(self.Name), ')'])

    def __getitem__(self, args):
        """__getitem__ passes everything to internal array.
        
        WARNING: m[a,b] will work where a and b are symbols in the alphabet,
        but m[a][b] will fail. This is because m[a] produces an array object
        with the corresponding row, which is then passed b as an index. Because
        the array object doesn't have the alphabet, it can't map the index into
        a number.

        Slicing is not supported.
        """
        # First, test whether args are in the JointAlphabet. Will always be tuple.
        if isinstance(args, tuple):
            try:
                return ravel(self._data)[self.Alphabet.index(tuple(args))]
            except (KeyError, TypeError):
                pass
        return self._data[self.Alphabet.SubEnumerations[0].index(args)]

    def __len__(self):
        """Returns number of rows."""
        return len(self._data)

    def empty(cls, Alphabet):
        """Class method: returns empty matrix sized for alphabet."""
        return cls(zeros(Alphabet.Shape), Alphabet)

    empty = classmethod(empty)

    def __eq__(self, other):
        """Tests whether two Usage objects have the same data."""
        try:
            return all(self._data == other._data)
            #return not bool(all(self._data != other._data))
        except:
            return False

    def __ne__(self, other):
        """Returns True if self and other are not equal."""
        try:
            return any(self._data != other._data)
            #return bool(all(self._data != other._data))
        except:
            return False

    def __iter__(self):
        """Iterates over rows in data."""
        return iter(self._data)

class Counts(PairMatrix):
    """Holds the data for a matrix of counts. Immutable.
    """
    
    def toProbs(self):
        """Returns copy of self where rows sum to 1."""
        return Probs(self._data/ (sum(self._data, 1)[:,NewAxis]), \
            self.Alphabet)

    def fromPair(cls, first, second, Alphabet, average=True):
        """Class method: returns new Counts from two sequences.
        """
        size = len(Alphabet.SubEnumerations[-1])
        #if they're ModelSequence objects, use the _data attribute
        if hasattr(first, '_data'):
            first, second = first._data, second._data

        #figure out what size we need the result to go in: note that the
        #result is on a pair alphabet, so the data type of the single
        #alphabet (that the sequence starts off in) might not work.
        data_type = get_array_type(product(map(len, Alphabet.SubEnumerations)))
        first = asarray(first, data_type)
        second = asarray(second, data_type)
        items = first * size + second
        
        counts = reshape(Alphabet.counts(items), Alphabet.Shape)
        if average:
            return cls((counts + transpose(counts))/2.0, Alphabet)
        else:
            return cls(counts, Alphabet)

    fromPair = classmethod(fromPair)

    def _from_triple_small(cls, first, second, outgroup, Alphabet):
        """Class method: returns new Counts for first from three sequences.

        Sequence order is first, second, outgroup.

        Use this method when the sequences are short and/or the alphabet is
        small: relatively memory intensive because it makes an array the size
        of the seq x the alphabet for each sequence. Fast on short sequences,
        though.

        NOTE: requires input to either all be ModelSequence objects, or all not
        be ModelSequence objects. Could change this if desirable.
        """
        #if they've got data, assume ModelSequence objects. Otherwise, arrays.
        if hasattr(first, '_data'):
            first, second, outgroup = first._data, second._data, outgroup._data

        size = len(Alphabet.SubEnumerations[-1])
        a_eq_b = equal(first, second)
        a_ne_b = logical_not(a_eq_b)
        a_eq_x = equal(first, outgroup)
        b_eq_x = equal(second, outgroup)

        #figure out what size we need the result to go in: note that the
        #result is on a pair alphabet, so the data type of the single
        #alphabet (that the sequence starts off in) might not work.
        data_type = get_array_type(product(map(len, Alphabet.SubEnumerations)))
        first = asarray(first, data_type)
        second = asarray(second, data_type)
       
        b_to_a = second*size + first
        a_to_a = first*size + first

        b_to_a_items = compress(logical_and(b_eq_x, a_ne_b), b_to_a)
        a_to_a_items = compress(logical_or(a_eq_b, a_eq_x), a_to_a)
        items = concatenate((b_to_a_items, a_to_a_items))
        counts = reshape(Alphabet.counts(items), Alphabet.Shape)

        return cls(counts, Alphabet)

    def _from_triple_large(cls, first, second, outgroup, Alphabet):
        """Same as _from_triple except copes with very long sequences.
        
        Specifically, allocates an array for the frequencies of each type,
        walks through the triple one base at a time, and updates the
        appropriate cell. Faster when alphabet and/or sequences are large;
        also avoids memory issues because it doesn't allocate the seq x
        alphabet array.

        NOTE: requires input to either all be ModelSequence objects, or all not
        be ModelSequence objects. Could change this if desirable.

        WARNING: uses float, not int, as datatype in return value.
        """
        #figure out if we already have the data in terms of alphabet indices.
        #if not, we need to convert it.
        if hasattr(first, '_data'):
            first, second, outgroup = first._data, second._data, outgroup._data
        else:
            if hasattr(Alphabet, 'toIndices'):
                converter = Alphabet.toIndices
            else:
                converter = Alphabet.fromSequenceToArray

            # convert to alphabet indices
            first, second, outgroup = map(asarray, map(converter,
                                        [first, second, outgroup]))
        # only include positions where all three not different
        valid_posn = logical_not(logical_and(logical_and(first != outgroup,
                                                        second != outgroup),
                                                        first != second))
        valid_pos = [index for index, val in enumerate(valid_posn) if val]
        first = first.take(valid_pos)
        second = second.take(valid_pos)
        outgroup = outgroup.take(valid_pos)
        out_diffs = logical_and(first == second, first != outgroup)
        counts = zeros((len(Alphabet.SubEnumerations[0]), \
            len(Alphabet.SubEnumerations[0])))
        for x, y, out_diff in zip(outgroup, first,
                                       out_diffs):
            if out_diff:
                counts[y,y] += 1
            else:
                counts[x,y] += 1
        return cls(counts, Alphabet)

    def fromTriple(cls, first, second, outgroup, Alphabet, threshold=1e6):
       """Reads counts from triple of sequences, method chosen by data size."""
       if len(first) * len(Alphabet) > threshold:
           return cls._from_triple_large(first, second, outgroup, Alphabet)
       else:
           return cls._from_triple_small(first, second, outgroup, Alphabet)

    fromTriple = classmethod(fromTriple)
    _from_triple_small = classmethod(_from_triple_small)
    _from_triple_large = classmethod(_from_triple_large)
       
class Probs(PairMatrix):
    """Holds the data for a probability matrix. Immutable."""
    
    def isValid(self):
        """Returns True if all values positive and each row sums to 1."""
        for row in self:
            if sum(row) != 1.0 or min(row) < 0.0:
                return False
        return True

    def makeModel(self, seq):
        """Returns substitution model for seq based on self's rows."""
        return take(self._data, seq, axis=0)

    def mutate(self, seq, random_vector=None):
        """Returns mutated version of seq, according to self.

        seq should behave like a Numeric array.
        
        random_vector should be vector of 0 and 1 of same length as sequence,
        if supplied.

        Result is always an array, not coerced into seq's class.
        """
        sums = cumsum(self._data, 1)
        model = take(sums, seq, axis=0)
        if random_vector is None:
            random_vector = randarray(seq.shape)
        return sum(transpose(model)[:-1] < random_vector, axis=0)
        #transpose needed to align frames
        

    def toCounts(self, num):
        """Returns count matrix with approximately num counts.

        Rounding error may prevent counts from summing exactly to num.
        """
        num_rows = len(self)
        return Counts(self._data * (num/num_rows), self.Alphabet)

    def toRates(self, normalize=False):
        """Returns rate matrix. Does not normalize by default."""
        return Rates(logm(self._data), self.Alphabet, self.Name, normalize)

    def random(cls, Alphabet, diags=None):
        """Makes random P-matrix with specified diag elements and size.

        diags can be a single float, or vector of values with same number
        of chars as individual alphabet (e.g. list of 4 elements will act
        as elements for the 4 bases).
        """
        shape = Alphabet.Shape
        if diags is None:
            result = randarray(shape)
            return cls(result/sum(result, 1)[:,NewAxis], Alphabet)
        else:
            single_size = shape[0]
            diags = array(diags, 'd')
            #handle scalar case
            if not diags.shape:
                diags = reshape(diags, (1,))
            if len(diags) == 1:
                diags = repeat(diags, single_size)
            temp = randarray((single_size, single_size-1))
            temp *= ((1.0-diags)/sum(temp, 1))[:,NewAxis]
            result = diag(diags)
            for r, row in enumerate(temp):
                result[r][:r] = row[:r]
                result[r][r+1:] = row[r:]
            return cls(result, Alphabet)

    random = classmethod(random)

class Rates(PairMatrix):
    """Holds the data for a rate matrix. Immutable."""

    def __init__(self, data, Alphabet, name=None, normalize=False):
        """Returns new Rates matrix, normalizing trace to -1 if necessary."""
        data = array(data)
        #check for complex input array
        if data.dtype == 'complex128':
            self.imag = data.imag
            data = data.real
        super(Rates, self).__init__(data, Alphabet)
        if normalize:
            self._normalize_inplace()

    def isComplex(self):
        """Returns True if self has a complex component."""
        return hasattr(self, 'imag')

    def isSignificantlyComplex(self, threshold=0.1):
        """Returns True if complex component is above threshold."""
        if hasattr(self, 'imag'):
            return sum(ravel(self.imag)) > threshold
        else:
            return False

    def isValid(self, threshold=1e-7):
        """Rate matrix is valid if rows sum to 0 and no negative off-diags.
        
        threshold gives maximum error allowed in row sums.
        """
        if max(abs(sum(self._data, -1)) > threshold):
            return False
        return not has_neg_off_diags(self._data)
    
    def _normalize_inplace(self):
        """Normalizes trace to -1, in-place.
        
        Should only call during __init__, since it mutates the object.
        WARNING: Only normalizes real component.
        """
        scale_trace(self._data)

    def normalize(self):
        """Returns normalized copy of self where trace is -1.
        
        WARNING: Only normalizes real component.
        """
        return Rates(self._data, self.Alphabet, normalize=True)

    def _get_diagonalized(self):
        """Gets diagonalization of self as u, v, w; caches values."""
        if not hasattr(self, '_diag_cache'):
            error_tolerance = 1e-4  #amount of error allowed in product
            eigenvalues, eigenvectors = eig(self._data)
            u = transpose(eigenvectors)
            v = eigenvalues
            w = inverse(u)
            #check that the diagonalization actually worked by multiplying
            #the results back together
            result = dot(dot(u,v),w)
            if abs(sum(ravel(result))) > error_tolerance:
                raise ValueError, "Diagonalization failed with erroneous result."
            self._diag_cache = u, v, w
        return self._diag_cache

    _diagonalized = property(_get_diagonalized)
        
    def toProbs(self, time=1.0):
        """Returns probs at exp(self*scale_factor).
        
        The way this works is by diagonalizing the rate matrix so that u is
        the matrix with eigenvectors as columns, v is a vector of eigenvalues,
        and w is the inverse of u. u * diag(v) * w reconstructs the original
        rate matrix. u * diag(exp(v*t)) * w exponentiates the rate matrix to
        time t.

        This is more expensive than a single exponentiation if the rate matrix
        is going to be sxponentiated only once, but faster if it is to be
        exponentiated to many different time points.

        Note that the diagonalization is not the same as the svd.

        If the diagonalization fails, we use the naive version of just
        multiplying the rate matrix by the time and exponentiating.
        """
        try:
            u, v, w = self._diagonalized
            #scale v to the right time by exp(v_0*t)
            v = diag(exp(v * time))
            return Probs(dot(dot(u,v), w), self.Alphabet)
        except:
            return Probs(expm(self._data)(time), self.Alphabet)

    def _timeForSimilarity_naive(self, similarity, freqs=None):
        """Returns time exponent so that exp(q*time) diverges to right distance.

        Takes symbol freqs into account if specified; otherwise assumes equal.

        freqs: vector of frequencies, applied to each row successively.

        WARNING: Factor of 5 slower than timeForSimilarity. Included for 
        testing that results are identical.
        """
        q = self._data
        if freqs is None:
            def similarity_f(t):
                return abs(average(diagonal(expm(q)(t)))-similarity)
        else:
            def similarity_f(t):
                return abs(sum(diagonal(expm(q)(t)*freqs)) - similarity)
        initial_guess = array([1.0])
        result = fmin(similarity_f, initial_guess, disp=0)
        #disp=0 turns off fmin messages
        return result

    def timeForSimilarity(self, similarity, freqs=None):
        """Returns time exponent so that exp(q*time) diverges to right distance.

        Takes symbol freqs into account if specified; otherwise assumes equal.

        freqs: vector of frequencies, applied to each row successively.

        NOTE: harder to understand, but a factor of 5 faster than the naive
        version. The nested matrixmultiply calls have the same effect as
        exponentiating the matrix.
        """
        #if there's no change, the time is 0
        if similarity == 1:
            return 0.0
        #try fast version first, but if it fails we'll use the naive version.
        try:
            u, v, w = self._diagonalized
            if freqs is None:
                def similarity_f(t):
                    return abs(average(diagonal(dot(u, \
                    dot(diag(exp(v*t)), w)))) - similarity)
            else:
                def similarity_f(t):
                    return abs(sum(diagonal(dot(u, \
                    dot(diag(exp(v*t)), w)))*freqs) - similarity)
        except (TypeError, ValueError):
            #get here if diagonalization fails
            q = self._data
            if freqs is None:
                def similarity_f(t):
                    return abs(average(diagonal(expm(q)(t)))-similarity)
            else:
                def similarity_f(t):
                    return abs(sum(diagonal(expm(q)(t)*freqs))-similarity)
        return brent(similarity_f)

    def toSimilarProbs(self, similarity, freqs=None):
        """Returns Probs at specified divergence.

        Convenience wrapper for toProbs and timeForSimilarity.
        """
        return self.toProbs(self.timeForSimilarity(similarity, freqs))

    def random(cls, Alphabet, diags=None):
        """Makes random Q-matrix with specified diag elements and size.

        diags can be a single float, or vector of values with same number
        of chars as individual alphabet (e.g. list of 4 elements will act
        as elements for the 4 bases).
        """
        shape = Alphabet.Shape
        single_size = shape[0]
        if diags is None:
            diags = -randarray(single_size)
        else:
            diags = array(diags, 'd')
            #handle scalar case
            if not diags.shape:
                diags = reshape(diags, (1,))
            if len(diags) == 1:
                diags = repeat(diags, single_size)
        temp = randarray((single_size, single_size-1))
        temp *= ((-diags)/sum(temp, 1))[:,NewAxis]
        result = diag(diags)
        for r, row in enumerate(temp):
            result[r][:r] = row[:r]
            result[r][r+1:] = row[r:]
        return cls(result, Alphabet)

    random = classmethod(random)

    def hasNegOffDiags(self):
        """Returns True if any off-diagonal elements negative."""
        return has_neg_off_diags(self._data)

    def sumNegOffDiags(self):
        """Returns sum of negative off-diagonal elements."""
        return sum_neg_off_diags(self._data)

    def fixNegsDiag(self):
        """Returns copy of self w/o negative off-diags, using 'diag' heuristic.

        If a negative off-diagonal element is encountered, sets it to 0.

        Subtracts all the negative off-diagonals from the diagonal to preserve
        row sum = 0.
        """
        m = self._data.copy()
        #clip to 0
        m = choose(less(m, 0.), (m, 0.))
        for i, row in enumerate(m):
            row[i] = -sum(row)
        return self.__class__(m, self.Alphabet)

    def fixNegsEven(self):
        """Returns copy of self w/o negative off-diags, using 'even' heuristic.
        
        If a negative off-diagonal is encountered, sets it to 0.

        Distributes the negative score evenly among the other elements.
        """
        m = without_diag(self._data)
        for i, row in enumerate(m):
            is_neg = row < 0
            if any(is_neg):
                num_negs = sum(is_neg)
                sum_negs = sum(is_neg*row)
                is_not_neg = logical_not(is_neg)
                num_not_neg = sum(is_not_neg)
                new_row = (row + (sum_negs/(num_not_neg+1)))*is_not_neg
                m[i] = new_row
        return self.__class__(with_diag(m, -sum(m,1)), self.Alphabet)

    def _make_error_f(self, to_minimize):
        """Make error function whose minimization estimates q = ln(p)."""
        p = expm(self._data)(t=1)
        BIG = 1e10
        def result(q):
            new_q = reshape(q, (4,4))
            neg_sum = sum_neg_off_diags(new_q)
            p_new = expm(new_q)(t=1)
            return to_minimize(ravel(p), ravel(p_new)) - (BIG * neg_sum) \
                + (BIG * sum(abs(sum(new_q,1))))
        return result

    def fixNegsFmin(self, method=fmin, to_minimize=norm_diff, debug=False):
        """Uses an fmin method to find a good approximate q matrix.

        Possible values for method:
            
            fmin:           simplex method (the default)
            fmin_bfgs:      bfgs optimizer  #always produces negative elements!
            fmin_cg:        cg optimizer    #doesn't work!
            fmin_powell:    powell method   #doesn't work!
        """
        q = self._data
        #bail out if q is already ok to start with
        if not sum_neg_off_diags(q):
            return self
        err_f = self._make_error_f(to_minimize)
        initial_guess = q.copy()
        xmin = method(err_f, initial_guess.flat, disp=0)
        #disp=0 turns off messages
        new_q = reshape(xmin, self.Alphabet.Shape)[:]
        if debug:
            if sum_neg_off_diags(new_q):
                raise Exception, 'Made invalid Q matrix: %s' % q
        return self.__class__(new_q, self.Alphabet)

    def fixNegsConstrainedOpt(self, to_minimize=norm_diff, badness=1e6):
        """Uses constrained minimization to find approx q matrix.

        to_minimize: metric for comparing orig result and new result.

        badness: scale factor for penalizing negative off-diagonal values.
        """
        if not sum_neg_off_diags(self._data):
            return self
        q = ravel(without_diag(self._data))
        p = expm(self._data)(t=1)
        def err_f(q):
            new_q = reshape(array(q), (4,3))
            new_q = with_diag(new_q, -sum(new_q, 1))
            p_new = expm(new_q)(t=1)
            result = to_minimize(ravel(p), ravel(p_new))
            if q.min() < 0:
                result += -q.min() * badness
            return result
        a = array(q)
        xmin = fmin(func=err_f, x0=a, disp=0)
        r = reshape(xmin, (4,3))
        new_q = with_diag(r, -sum(r, 1))
        return self.__class__(new_q, self.Alphabet)

    def fixNegsReflect(self):
        """Fixes negative off-diagonals by subtracting m[i][j] from m[j][i].
        
        Specifically, if m[i][j] is negative, subtracts this value from
        m[i][j] and m[i][i] to keep the row total at 0, and then subtracts
        it from m[j][i] and m[j][j] to convert a negative flux in the forward
        direction into a positive flux in the reverse direction. If both
        m[i][j] and m[j][i] are negative, this algorithm converts them both
        into positive values, effectively exchanging the magnitudes of the
        changes and making the signs positive.

        NOTE: It's important to iterate over the original and make changes to
        the copy to avoid incorrect results in cases where both m[i][j] and
        m[j][i] are negative.
        """
        orig = self._data
        result = orig.copy()
        for i, row in enumerate(orig):
            for j, val in enumerate(row):
                #skip diagonal
                if i == j:
                    continue
                #only make changes if element < 0
                if val < 0:
                    result[i][j] -= val
                    result[i][i] += val
                    result[j][i] -= val
                    result[j][j] += val
        return self.__class__(result, self.Alphabet)

def goldman_q_rna_triple(seq1, seq2, outgroup):
    """Returns the Goldman rate matrix for seq1"""
    if len(seq1) != len(seq2) != len(outgroup):
        raise ValueError, "seq1,seq2 and outgroup are not the same length!"

    seq1 = ModelRnaSequence(seq1)
    seq2 = ModelRnaSequence(seq2)
    outgroup = ModelRnaSequence(outgroup)

    m = Counts.fromTriple(seq1, seq2, outgroup, RnaPairs)._data

    q = m / m.sum(axis=1)[:,NewAxis]
    new_diag = -(q.sum(axis=1) - diag(q))

    for i,v in enumerate(new_diag):
        q[i,i] = v

    return q

def goldman_q_dna_triple(seq1, seq2, outgroup):
    """Returns the Goldman rate matrix for seq1"""
    if len(seq1) != len(seq2) != len(outgroup):
        raise ValueError, "seq1,seq2 and outgroup are not the same length!"

    seq1 = ModelDnaSequence(seq1)
    seq2 = ModelDnaSequence(seq2)
    outgroup = ModelDnaSequence(outgroup)

    m = Counts.fromTriple(seq1, seq2, outgroup, DnaPairs)._data

    q = m / m.sum(axis=1)[:,NewAxis]
    new_diag = -(q.sum(axis=1) - diag(q))

    for i,v in enumerate(new_diag):
        q[i,i] = v

    return q

def goldman_q_dna_pair(seq1, seq2):
    """Returns the Goldman rate matrix"""
    if len(seq1) != len(seq2):
        raise ValueError, "seq1 and seq2 are not the same length!"

    seq1, seq2 = ModelDnaSequence(seq1), ModelDnaSequence(seq2)

    m = Counts.fromPair(seq1, seq2, DnaPairs,average=True)._data

    q = m / m.sum(axis=1)[:,NewAxis]
    new_diag = -(q.sum(axis=1) - diag(q))

    for i,v in enumerate(new_diag):
        q[i,i] = v

    return q

def goldman_q_rna_pair(seq1, seq2):
    """Returns the Goldman rate matrix"""
    if len(seq1) != len(seq2):
        raise ValueError, "seq1 and seq2 are not the same length!"

    seq1, seq2 = ModelRnaSequence(seq1), ModelRnaSequence(seq2)

    m = Counts.fromPair(seq1, seq2, RnaPairs,average=True)._data

    q = m / m.sum(axis=1)[:,NewAxis]
    new_diag = -(q.sum(axis=1) - diag(q))

    for i,v in enumerate(new_diag):
        q[i,i] = v

    return q

def make_random_from_file(lines):
    """Simulates array random() using values from an iterator."""
    def result(shape):
        size = product(shape)
        items = map(float, [lines.next() for s in range(size)])
        a = reshape(array(items), shape)
        return a
    return result


#randarray = make_random_from_file(open('/Users/rob/random.txt'))

def test_heuristics(p_range=None, num_to_do=71, heuristics=None):
    if p_range is None:
        p_range = [0.6]
    if heuristics is None:
        heuristics = ['fixNegsDiag', 'fixNegsEven', 'fixNegsReflect', 'fixNegsConstrainedOpt']
    num_heuristics = len(heuristics)
    print '\t'.join(['p'] + heuristics)
    for p in p_range:
        result = zeros((num_to_do, num_heuristics), Float64)
        has_nonzero = 0
        i = 0
        while i < num_to_do:
            curr_row = result[i]
            random_p = Probs.random(DnaPairs, p)
            q = random_p.toRates()
            if not q.hasNegOffDiags():
                continue
            has_nonzero += 1
            #print "P:"
            #print random_p._data
            #print "Q:"
            #print q._data
            i += 1
            for j, h in enumerate(heuristics):
                #print "HEURISTIC: ", h
                q_corr = getattr(q, h)()
                #print "CORRECTED Q: "
                #print q_corr._data
                p_corr = expm(q_corr._data)(t=1)
                #print "CORRECTED P:"
                #print p_corr
                dist = norm_diff(p_corr, random_p._data)
                #print "DISTANCE: ", dist
                curr_row[j] = dist
        averages = average(result)
        print p, '\t', '\t'.join(map(str, averages))

if __name__ == '__main__':
    test_heuristics()
                
        

#!/usr/bin/env python
#file cogent/maths/stats/util.py

"""Provides classes and utility methods for statistics.

Classes: Numbers, Freqs, NumberFreqs, SummaryStatistics.

Owner: Rob Knight rob@spot.colorado.edu

Status: Stable

Notes

SHOULD THIS FILE BE DELETED? Numbers is clearly superfluous since we're now
using numpy. Should some of the Freqs classes be kept? Should UnsafeFreqs
become Freqs, and the full MappedDict version of Freqs be deleted?
-- RK 8/2/05

The distinctions among the classes can be somewhat subtle, so here's a quick
guide to usage.

The core classes are Numbers, Freqs, NumberFreqs, and SummaryStatistics. For
each of the first three, there is also an Unsafe version (e.g. UnsafeFreqs),
which has the same interface but avoids overriding built-in methods for 
performance reasons. The performance differences can be very large (the unsafe
version can be an order of magnitude faster). However, the unsafe versions do
not validate their input, so methods will fail if the data are invalid.

SummaryStatistics holds a Count, Sum, Variance, Mean, StandardDeviation,
and SumSquares (all of these are optional). If initialized with only some of
these statistics, it will calculate the others when needed when it can (e.g.
it can calculate the Mean given the Sum and the Count). SummaryStatistics
is useful for holding information about a large data set when you need to
throw away the original data to free memory. The statistics functions all
work on properties, and most functions that work on Numbers or Freqs will
also work fine on the equivalent SummaryStatistics object.

Numbers holds a list of values that are all numbers (i.e. can be converted to
float -- it does not try to do anything with complex values). Numbers supports
the full list interface, and also supports methods like items (returns key, 
value pairs), toFixedWidth (for formatting), normalize, accumulate, 
firstIndexLessThan (and its relatives), randomSequence, round, and so forth.

UnsafeNumbers behaves like Numbers, except that it does not complain if you
initialize it with non-numeric values or insert such values later (with insert,
__setitem__, extend, and so forth). Also, UnsafeNumbers does not automatically
map strings to numbers during __init__.

Freqs holds a dict of key, value pairs where the keys are assumed to be 
separate categories. Statistics functions act on the categories: in other
words, Count returns the number of categories, Sum returns the number of
observations in all categories, Mean returns the mean number of observations
per category, and so on. Freqs (and its subclasses) notably supports addition
of counts from a sequence, a dict, a list of dicts, a sequence of key, value
pairs, and a sequence of sequences: in all cases, the counts are added (unlike
dict.update(), which overwrites the value with the last one read from that
key). Use +, +=, -, and -= to accumulate counts. Freqs supports all the stats
functions that Numbers does. Notable features of Freqs include randomSequence
(generates a random sequence of keys according to their frequencies), rekey
(maps the freqs onto new keys, given a conversion dictionary), normalize,
scale, round, and getSortedList (sorts items by key or value, ascending or
descending). Freqs refused to hold anything except non-negative floats (performs
conversion on addition), and checks this constraint on all methods that mutate
the object.

UnsafeFreqs is like Freqs except that it is non-validating. One important 
consequence is that __init__ behaves differently: UnsafeFreqs inherits the raw
dict __init__. This means that Freqs([('a',2),('b',3),('a',1)]) gives you an
object that compares equal to Freqs({'a':3,'b':3}) because it sums the repeated
keys, but UnsafeFreqs([('a',2),('b',3),('a',1)]) gives you an object that 
compares equal to Freqs({'a':1,'b':3}) because the last key encounted overwrites
the previous value for that key. Additionally, the UnsafeFreqs constructor 
can't accept all the things the Freqs constructor can -- the easiest workaround
is to create an empty UnsafeFreqs and use += to fill it with data, although
if you know the form in the data in advance it's much faster to use the 
appropriate method (e.g. fromTuples, fromDict) than to let += guess what you
gave it. UnsafeFreqs does not check that operations that mutate the dictionary
preserve validity.

NumberFreqs hold a dict of key, value pairs where the keys are assumed to be
numbers along an axis, and the values are assumed to be the counts at each
point. Thus, Count returns the number of _observations_ (not categories), Sum
returns the sum of key * value, Mean returns the mean value of the observations,
and so forth. An example of how this works is as follows:
    
    Freqs([1,2,2,1,1,1,3,3]) gives you the dict {1:4,2:2,3:2}. These values are
    interpreted as coming from 3 categories, so the Count is 3. There are 8
    observations, so the Sum is 8. The Mean across categories is 8/3.

    NumberFreqs([1,2,2,1,1,3,3]) gives you the dict {1:4,2:2,3:2}. These values
    are interpreted as 4 counts of 1, 2 counts of 2, and 2 counts of 3. All the
    values are treated as coming from the same category, so the Count is 8. The
    sum is calculated by weighting the value by the key, so is 1*4 + 2*2 + 3*2;
    thus, the Sum is 14. Consequently, the Mean of the values is 14/8.

Consequently, NumberFreqs is appropriate when you want to treat the distribution
of key, value pairs like a histogram in the keys, and find the average along
the key axis. Freqs is appropriate when you want to treat the distribution of
key, value pairs like a bar graph where the keys are separate categories, and
find the average along the value axis. This distinction between Freqs and
NumberFreqs holds for all the stats functions, which rely on the Sum, Count,
and other properties whose behavior differs between Freqs and NumberFreqs.

UnsafeNumberFreqs behaves like NumberFreqs except that it doesn't validate on
input or mutation of the dict. It's much faster, though.
"""
from __future__ import division
from cogent.util.misc import FunctionWrapper, MappedList, MappedDict, \
    ConstraintError
from cogent.util.table import Table
from numpy import array, sqrt, log2, e, floor, ceil
from random import choice, random
from operator import gt, ge, lt, le, add, sub

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit", "Gavin Huttley", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class SummaryStatisticsError(ValueError):
    """Raised when not possible to calculate a requested summary statistic."""
    pass

class SummaryStatistics(object):
    """Minimal statistics interface. Object is read-only once created."""
    def __init__(self, Count=None, Sum=None, Mean=None, StandardDeviation=None,
        Variance=None, SumSquares=None, Median=None):
        """Returns a new SummaryStatistics object."""
        self._count = Count
        self._sum = Sum
        self._mean = Mean
        self._standard_deviation = StandardDeviation
        self._variance = Variance
        self._sum_squares = SumSquares
        self._median = Median

    def __str__(self):
        """Returns string representation of SummaryStatistics object."""
        result = []
        for field in ["Count", "Sum", "Median", "Mean", "StandardDeviation", \
                       "Variance", "SumSquares"]:
            try:
                val = getattr(self, field)
                if not val:
                    continue
                result.append([field, val])
            except:
                pass
        if not result:
            return ''
        return str(Table("Statistic Value".split(), result,
            column_templates={'Value': "%.4g"}))
    
    def _get_count(self):
        """Returns Count if possible (tries to calculate as sum/mean)."""
        if self._count is None:
            try:
                self._count = self._sum/self._mean
            except (TypeError, ZeroDivisionError, FloatingPointError):
                raise SummaryStatisticsError, \
                    "Insufficient data to calculate count."
        return self._count
    Count = property(_get_count)

    def _get_sum(self):
        """Returns Sum if possible (tries to calculate as count*mean)."""
        if self._sum is None:
            try:
                self._sum = self._count * self._mean
            except TypeError:
                raise SummaryStatisticsError, \
                    "Insufficient data to calculate sum."
        return self._sum
    Sum = property(_get_sum)

    def _get_mean(self):
        """Returns Mean if possible (tries to calculate as sum/count)."""
        if self._mean is None:
            try:
                self._mean = self._sum / self._count
            except (TypeError, ZeroDivisionError, FloatingPointError):
                raise SummaryStatisticsError, \
                    "Insufficient data to calculate mean."
        return self._mean
    Mean = property(_get_mean)
    
    def _get_median(self):
        """Returns Median."""
        return self._median
    Median = property(_get_median)
    
    def _get_standard_deviation(self):
        """Returns StandardDeviation if possible (calc as sqrt(var)."""
        if self._standard_deviation is None:
            try:
                self._standard_deviation = sqrt(abs(self._variance))
            except TypeError:
                raise SummaryStatisticsError, \
                    "Insufficient data to calculate standard deviation."
        return self._standard_deviation
    StandardDeviation = property(_get_standard_deviation)

    def _get_variance(self):
        """Returns Variance if possible (calculates as stdev ** 2)."""
        if self._variance is None:
            try:
                self._variance = self._standard_deviation * \
                    self._standard_deviation
            except TypeError:
                raise SummaryStatisticsError, \
                    "Insufficient data to calculate variance."
        return self._variance
    Variance = property(_get_variance)
    
    def _get_sum_squares(self):
        """Returns SumSquares if possible."""
        if self._sum_squares is None:
            raise SummaryStatisticsError, \
                "Insufficient data to calculate sum of squares."
        return self._sum_squares
    SumSquares = property(_get_sum_squares)
    
    def __cmp__(self, other):
        """SummaryStatistics compares by count, then sum, then variance.
        
        Absent values compare as 0.
        """
        result = 0
        for attr in ['Count', 'Sum', 'Variance', 'SumSquares']:
            try:
                my_attr = getattr(self, attr)
            except SummaryStatisticsError:
                my_attr = 0
            try:
                other_attr = getattr(other, attr)
            except SummaryStatisticsError:
                other_attr = 0
            result = result or cmp(my_attr, other_attr)
        return result

class NumbersI(object):
    """Interface for Numbers, a list that performs numeric operations."""
    
    _is_sorted = False
    
    def isValid(self):
        """Checks that all items in self are numbers."""
        for i in self:
            if not (isinstance(i, int) or isinstance(i, float)):
                return False
        return True
    
    def items(self):
        """Returns list of (item, 1) tuples for each item in self.
        
        This is necessary because we want to delegate attribute accesses
        (specifically, method calls for stats functions) to a 
        Freqs object. Freqs tests whether an
        object is dictionary-like by calling items() on it, expecting an error
        if it doesn't have the method. However, NumericList delegates items()
        back to a new Freqs, calling the constructor in a
        cycle...

        The workaround, since we don't want to explicitly specify which items
        get passed on to Freqs, is just to give Numbers its
        own items() method.
        """
        return zip(self, [1] * len(self))
  
    def toFixedWidth(self, fieldwidth=10):
        """Returns string with elements mapped to fixed field width.
        
        Always converts to scientific notation. Minimum fieldwidth is 7,
        since it's necessary to account for '-xe-yyy'.

        Result has (fieldwidth - 7) significant figures of precision, or
        1 if fieldwidth is 7.
        """
        if fieldwidth < 7:
            raise ValueError, "toFixedWidth requres fieldwith of at least 8."
        if fieldwidth == 7:
            decimals = 0
        else:
            decimals = fieldwidth-8
        format = ''.join(["%+", str(fieldwidth), '.',str(decimals),'e'])
        return ''.join( [format % i for i in self])

    def normalize(self, x=None):
        """Normalizes items in Numbers by dividing by x (sum by default)."""
        if not self:
            return      #do nothing of empty
        if x is None:
            x = self.Sum
        if not x:       #do nothing if items are empty
            return
        x = float(x)
        for index, item in enumerate(self):
            self[index] = item/x

    def accumulate(self):
        """Converts self to cumulative sum, in place"""
        if self:
            curr = self[0]
            for i in xrange(1, len(self)):
                curr += self[i]
                self[i] = curr

    def firstIndexGreaterThan(self, value, inclusive=False, stop_at_ends=False):
        """Returns first index of self that is greater than value.

        inclusive: whether to use i > value or i >= value for the test.

        stop_at_ends: whether to return None or the last index in self if none 
        of the items in self are greater than the value.
        """
        if inclusive:
            operator = ge
        else:
             operator = gt
       
        for i, curr in enumerate(self):
            if operator(curr, value):
                return i
        #only get here if we didn't find anything greater
        if stop_at_ends:
            return i
        #default is to return None
             
    def firstIndexLessThan(self, value, inclusive=False, stop_at_ends=False):
        """Returns first index of self that is less than value.

        inclusive: whether to use i < value or i <= value for the test.

        stop_at_ends: whether to return None or the last index in self if none
        of the items in self are less than the value.
        """
        if inclusive:
            operator = le
        else:
            operator = lt
       
        for i, curr in enumerate(self):
            if operator(curr, value):
                return i
        #only get here if we didn't find anything greater
        if stop_at_ends:
            return i
        #default is to return None

    def lastIndexGreaterThan(self, value, inclusive=False, stop_at_ends=False):
        """Returns last index of self that is greater than value.

        inclusive: whether to use i > value or i >= value for the test.

        stop_at_ends: whether to return None or 0 if none of the items in self
        is greater than the value.
        """
        if inclusive:
            operator = ge
        else:
            operator = gt
       
        latest = None
       
        for i, curr in enumerate(self):
            if operator(curr, value):
                latest = i
        if stop_at_ends and (latest is None):
            return 0
        else:
            return latest

    def lastIndexLessThan(self, value, inclusive=False, stop_at_ends=False):
        """Returns last index of self that is less than value.

        inclusive: whether to use i < value or i <= value for the test.

        stop_at_ends: whether to return None or 0 if none of the items in self
        is less than the value.
        """
        if inclusive:
            operator = le
        else:
            operator = lt
       
        latest = None
       
        for i, curr in enumerate(self):
            if operator(curr, value):
                latest = i
        if stop_at_ends and (latest is None):
            return 0
        else:
            return latest
  
    def _get_sum(self):
        """Returns sum of items in self."""
        return sum(self)

    Sum = property(_get_sum)

    Count = property(list.__len__)
    
    def _get_sum_squares(self):
        """Returns sum of squares of items in self."""
        return sum([i*i for i in self])
    SumSquares = property(_get_sum_squares)

    def _get_variance(self):
        """Returns sample variance of items in self.
        
        Fault-tolerant: returns zero if one or no items.
        """
        if not self:
            return None
        total = self.Sum
        count = self.Count
        if count <= 1:  #no variance for a single item
            variance = 0.0
        else:
            variance = (self.SumSquares-(total*total)/count)/(count-1)
        return variance
    Variance = property(_get_variance)

    def _get_standard_deviation(self):
        """Returns sample standard deviation of items in self."""
        if not self:
            return None
        return sqrt(abs(self.Variance))
    StandardDeviation = property(_get_standard_deviation)

    def _get_mean(self):
        """Returns mean of items in self."""
        if not self:
            return None
        try:
            mean = self.Sum/self.Count
        except (ZeroDivisionError, FloatingPointError):
            mean = 0.0
        return mean
    Mean = property(_get_mean)
    
    def _get_mode(self):
        """Returns the most frequent item. If a tie, picks one at random.
        
        Usage: most_frequent = self.mode()

        """
        best = None
        best_count = 0
        for item, count in self.items():
            if count > best_count:
                best_count = count
                best = item
        return best
    Mode = property(_get_mode)
    
    def quantile(self, quantile):
        """Returns the specified quantile.
        
        Uses method type 7 from R. Only sorts on first call, so subsequent
        modifications may result in incorrect estimates unless directly sorted
        prior to using."""
        if not self._is_sorted:
            self.sort()
            self._is_sorted = True
        index = quantile * (len(self)-1)
        lo = int(floor(index))
        hi = int(ceil(index))
        diff = index - lo
        stat = (1-diff) * self[lo] + diff * self[hi]
        return stat
    
    def _get_median(self):
        """Returns the median"""
        return self.quantile(0.5)
    
    Median = property(_get_median)
    
    def summarize(self):
        """Returns summary statistics for self."""
        return SummaryStatistics(Count=self.Count, Sum=self.Sum, \
            Variance=self.Variance, Median = self.Median)
        
    def choice(self):
        """Returns random element from self."""
        return choice(self)

    def randomSequence(self, n):
        """Returns list of n random choices from self, with replacement."""
        return [choice(self) for i in range(n)]

    def subset(self, items, keep=True):
        """Retains (or deletes) everything in self contained in items, in place.

        For efficiency, items should be a dict.
        """
        if keep:
            self[:] = filter(items.__contains__, self)
        else:
            self[:] = filter(lambda x: x not in items, self)

    def copy(self):
        """Returns new copy of self, typecast to same class."""
        return self.__class__(self[:])

    def round(self, ndigits=0):
        """Rounds each item in self to ndigits."""
        self[:] = [round(i, ndigits) for i in self]

    #following properties/methods are handled by conversion to Freqs.
    # Uncertainty, Mode

    def _get_uncertainty(self):
        """Returns Shannon entropy of items in self."""
        return UnsafeFreqs().fromSeq(self).Uncertainty

    Uncertainty = property(_get_uncertainty)

    def _get_mode(self):
        """Returns most common element in self."""
        return UnsafeFreqs().fromSeq(self).Mode

    Mode = property(_get_mode)

class UnsafeNumbers(NumbersI, list):
    """Subclass of list that should only hold floating-point numbers.

    Usage: nl = Numbers(data)

    WARNING: UnsafeNumbers does not check that items added into it are really
    numbers, and performs almost no validation.
    """
    pass

class Numbers(NumbersI, MappedList):
    """Safe version of Numbers that validates on all list operations.

    For each item in data (which must be iterable), tests whether the item
    is a number and, if so, adds it to the Numbers.

    Note: this means we have to override _all_ the list methods that might
    potentially add new data to the list. This makes it much slower than 
    UnsafeNumbers, but impossible for it to hold invalid data.
    """
    Mask = FunctionWrapper(float)
    
    def __init__(self, data=None, Constraint=None, Mask=None):
        """Initializes a new Numbers object.

        Usage: nl = Numbers(data)

        For each item in data, tries to convert to a float. If successful,
        produces new Numbers with data.

        Note: this means that a single string of digits will be treated as
        a list of digits, _not_ as a single number. This might not be what
        you expected.

        Also, data must be iterable (so a 1-element list containing a number
        is OK, but a single number by itself is not OK).
        """
        if data is not None:
            data = map(float, data) #fails if any items are not floatable
        else:
            data = []
        MappedList.__init__(self, data, Constraint, Mask)
 
class FreqsI(object):
    """Interface for frequency distribution, i.e. a set of value -> count pairs.
    """
    RequiredKeys = {}

    def fromTuples(self, other, op=add, uses_key=False):
        """Adds counts to self inplace from list of key, val tuples.

        op: operator to apply to old and new values (default is add,
        but sub and mul might also be useful. Use names from the operator
        module).

        uses_key: whether or not op expects the key as the first argument,
        i.e. it gets (key, old, new) rather than just (old, new).

        Returns modified version of self as result.
        """
        for key, val in other:
            curr = self.get(key, 0)
            if uses_key:
                self[key] = op(key, curr, val)
            else:
                self[key] = op(curr, val)
        return self

    def newFromTuples(cls, other, op=add, uses_key=False):
        """Classmethod: returns new FreqsI object from tuples."""
        result = cls()
        return result.fromTuples(other, op, uses_key)

    newFromTuples = classmethod(newFromTuples)

    def fromDict(self, other, op=add, uses_key=False):
        """Adds counts to self inplace from dict of {key:count}.

        op: operator to apply to old and new values (default is add,
        but sub and mul might also be useful. Use names from the operator
        module).
        
        Returns modified version of self as result.
        """
        for key, val in other.items():
            curr = self.get(key, 0)
            if uses_key:
                self[key] = op(key, curr, val)
            else:
                self[key] = op(curr, val)
        return self

    def newFromDict(cls, other, op=add, uses_key=False):
        """Classmethod: returns new FreqsI object from single dict."""
        result = cls()
        return result.fromDict(other, op, uses_key)

    newFromDict = classmethod(newFromDict)

    def fromDicts(self, others, op=add, uses_key=False):
        """Adds counts to self inplace from list of dicts of {key:count}.

        op: operator to apply to old and new values (default is add,
        but sub and mul might also be useful. Use names from the operator
        module).

        Returns modified version of self as result.
        """
        for i in others:
            self.fromDict(i, op, uses_key)
        return self

    def newFromDicts(cls, others, op=add, uses_key=False):
        """Classmethod: returns new FreqsI object from single dict."""
        result = cls()
        return result.fromDicts(others, op, uses_key)

    newFromDicts = classmethod(newFromDicts)

        
    def fromSeq(self, seq, op=add, weight=1, uses_key=False):
        """Adds counts to self inplace from seq. Each item adds 'weight' counts.

        op: operator to apply to old and new values (default is add,
        but sub and mul might also be useful. Use names from the operator
        module).

        weight: increment to apply each time an item is found.

        Applies self[i] += op(curr, weight) for each item in seq.

        Returns modified version of self as result.
        """
        if uses_key:
            for i in seq:
                curr = self.get(i, 0)
                self[i] = op(i, curr, weight)
        else:
            for i in seq:
                curr = self.get(i, 0)
                self[i] = op(curr, weight)
        return self

    def newFromSeq(cls, seq, op=add, weight=1, uses_key=False):
        """Classmethod: returns new FreqsI object from single dict."""
        result = cls()
        r = result.fromSeq(seq, op, weight, uses_key)
        return result.fromSeq(seq, op, uses_key)

    newFromSeq = classmethod(newFromSeq)

    def fromSeqs(self, seqs, op=add, weight=1, uses_key=False):
        """Adds counts to self inplace from seq of sequences. Uses weight.

        op: operator to apply to old and new values (default is add).

        weight: increment to apply each time an item is found.
        """
        for s in seqs:
            self.fromSeq(s, op, weight, uses_key)
        return self

    def newFromSeqs(cls, seqs, op=add, weight=1, uses_key=False):
        """Classmethod: returns new FreqsI object from single dict."""
        result = cls()
        return result.fromSeqs(seqs, op, weight, uses_key)

    newFromSeqs = classmethod(newFromSeqs)


    def _find_conversion_function(self, data):
        """Figures out which conversion function to use for data."""
        # if the data is empty, it's easy...
        if not data:
            return None
        # if it has items, treat as a dict (fromDict only uses items() from the
        # dict interface).
        if hasattr(data, 'items'):
            return self.fromDict
        # if it's a string, treat it as a sequence
        if isinstance(data, str):
            return self.fromSeq
        # Otherwise, we need to check what it is. Must be a sequence of some
        # kind: elements could be dicts, tuples, or second-level sequences.
        # We know it's not empty, so we can get the first element
        first = data[0]
        #if the first item is a dict, assume they all are
        if isinstance(first, dict):
            return self.fromDicts
        # otherwise, if all items have two elements and the second is a number,
        # assume they're key-value pairs
        try:
            for key, value in data:
                v = float(value)
            # if it did work, data can be treated as a sequence of key-value
            # pairs. Note that if you _really_ have e.g. a sequence of pairs
            # of numbers that you want to add individually, you need to use
            # fromSeqs explicitly.
            return self.fromTuples
        except (TypeError, ValueError):
            # if that didn't work, data is either a sequence or a sequence of
            # sequences.
            # if first is iterable and not a string, treat as seq of seqs; 
            # otherwise, treat as seq.
            # Note that this means that lists of strings will always be treated
            # as though each string is a key that's being counted (e.g. if you
            # pass in a list of words, you'll get word frequencies rather than
            # character frequencies). If you want the character frequencies,
            # call fromSeqs explicitly -- there's no way to detect what's 
            # desired automatically.
            if isinstance(first, str):
                return self.fromSeq
            else:
                # if first item is iterable, treat as seq of seqs. otherwise,
                # treat as single seq of items.
                try:
                    i = iter(first)
                    return self.fromSeqs
                except TypeError:
                    return self.fromSeq
        # should never get here because of return values
        raise NotImplementedError, "Fell off end of _find_conversion_function"
                
    def isValid(self):
        """Checks presence of required keys, and that all vals are numbers."""
        for k in self.RequiredKeys:
            if k not in self:
                return False
        for v in self.values():
            if not (isinstance(v, float) or isinstance(v, int)):
                return False
            if v < 0:
                return False
        return True

    def __iadd__(self, other):
        """Adds items from other to self."""
        f = self._find_conversion_function(other)
        if f:   #do nothing if we got None, since it means other was empty
            f(other, op=add)
        return self

    def __add__(self, other):
        """Returns new Freqs object with counts from self and other."""
        result = self.copy()
        result += other
        return result

    def __isub__(self, other):
        """Subtracts items in other from self."""
        f = self._find_conversion_function(other)
        if f:   #do nothing if we got None, since it means other was empty
            f(other, op=sub)
        return self

    def __sub__(self, other):
        """Returns new Freqs containing difference between self and other."""
        result = self.copy()
        result -= other
        return result
          
    def __str__(self):
        """Prints the items of self out as tab-delimited text.
        
        Value, Frequency pairs are printed one pair to a line. Headers are 
        'Value' and 'Count'.
        """
        if self:
            lines = ["Value\tCount"]
            items = self.items()  #make and sort list of (key, value) pairs
            items.sort()
            for key, val in items:
                lines.append("\t".join([str(key), str(val)]))  #add pair
            return "\n".join(lines)
        else:
            return "Empty frequency distribution"

    def __delitem__(self, key):
        """May not delete key if it is required.
        
        WARNING: Will not work if your class doesn't subclass dict as well
        as FreqsI.
        """
        r = self.RequiredKeys
        if r and (key in r):
            raise KeyError, "May not delete required key %s" % key
        else:
            dict.__delitem__(self, key)

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
            result[new_key] = curr + val
        return result

    def purge(self):
        """If self.RequiredKeys is nonzero, deletes everything not in it."""
        req = self.RequiredKeys
        if req:
            for key in self.keys():
                if key not in req:
                    del self[key]

    def normalize(self, total=None, purge=True):
        """Converts all the counts into probabilities, i.e. normalized to 1.

        Ensures that all the required keys are present, if self.RequiredKeys
        is set.
        
        Does the transformation in place.

        If purge is True (the default), purges any keys that are not in
        self.RequiredKeys before normalizing.

        Can also pass in a number to divide by instead of using the total,
        which is useful for getting the average across n data sets.
        
        Usage: self.normalize()

        """
        if purge:
            self.purge()

        req = self.RequiredKeys
        if req:
            for r in req:
                if r not in self:
                    self[r] = 0
            
        if total is None:
            total = self.Sum
        if total != 0:          #avoid divide by zero
            for item, freq in self.items():
                f = float(freq)
                if f < 0:
                    raise ValueError, \
                    "Freqs.normalize(): found negative count!"
                self[item] = f/total

    def choice(self, prob):
        """If self is normalized, returns item corresponding to Pr(prob)."""
        sum = 0
        items = self.items()
        for item, freq in items:
            sum += freq
            if prob <= sum:
                return item
        return items[-1][0]    #return the last value if we run off the end

    def randomSequence(self, n):
        """Returns list of n random choices, with replacement.
        
        Will raise IndexError if there are no items in self.
        """
        num_items = self.Sum
        return [self.choice(random()*num_items) for i in xrange(n)]

    def subset(self, items, keep=True):
        """Deletes keys for all but items from self, in place."""
        if keep:
            to_delete = []
            for i in self:
                try: #fails if i is not a string but items is
                    delete = i not in items
                except TypeError: #i was wrong type, so can't be in items... 
                    to_delete.append(i)
                else:
                    if delete:
                        to_delete.append(i)
        else:
            to_delete = items
            
        for i in to_delete:
            try:
                del self[i]
            except KeyError:
                pass #don't care if it wasn't in the dictionary
    
    def scale(self, factor=1, offset=0):
        """Linear transform of values in freqs where val = facto*val + offset.
        
        Usage: f.scale(factor, offset)

        Does the transformation in place.
        """
        for k,v in self.items():
            self[k] = v * factor + offset

    def round(self, ndigits=0 ):
        """Rounds frequencies in Freqs to ndigits (default:0, i.e. integers).

        Usage: f.round()
        
        Does the transformation in place
        """
        for k,v in self.items():
            self[k] = round(v, ndigits)

    def expand(self, order=None, convert_to=None, scale=None):
        """Expands the Freqs into a sequence of symbols.

        Usage: f.expand(self, order=None, convert_to=list)

        order should be a sequence of symbols. Symbols that are not in f will
        be silently ignored (so it's safe to use on e.g. codon usage tables
        where some of the codons might not appear).

        Each item should appear in the list as many times as its frequency.

        convert_to should be a callable that takes an arbitrary sequence and
        converts it into the desired output format. Default is list, but
        ''.join is also popular.
        
        scale should be the number you want your frequencies multiplied with.
        Scaling only makes sense if your original freqs are fractions (and 
        otherwise it won't work anyway). The scaling and rounding are done
        on a copy of the original, so the original data is not changed.
        
        Calls round() on each frequency, so if your values are normalized you'll
        need to renormalize them (e.g. self.normalize(); self.normalize(1.0/100)
        to get percentages) or the counts will be zero for anything that's less
        frequent than 0.5. You _can_ use this to check whether any one symbol
        constitutes a majority, but it's probably more efficient to use 
        mode()...
        """
        if scale:
            if sum([round(scale*v) for v in self.values()]) != scale:
               raise ValueError,\
                    "Can't round to the desired number (%d)"%(scale)
            else:
                used_freq = self.copy()
                used_freq.scale(factor=scale)
                used_freq.round()
        else:
            used_freq = self
        if order is None:
            order = used_freq.keys()
        result = []
        for key in order:
            result.extend([key] * int(round(used_freq.get(key, 0))))
        if convert_to:
            return convert_to(result)
        else:
            return result

    def _get_count(self):
        """Calculates number of categories in the frequency distribution.
        
        Useful for other stats functions. Assumes that keys are categories.
        
        Usage: count = self.count()
        
        Note that for NumberFreqs, Count will instead return the total number
        of observations.
        """
        return len(self)
    Count = property(_get_count)
 
    def _get_sum(self):
        """Returns sum of items in self."""
        return sum(self.values())
    Sum = property(_get_sum)

    def _get_sum_squares(self):
        """Returns sum of squares of items in self."""
        return sum([i*i for i in self.values()])
    SumSquares = property(_get_sum_squares)

    def _get_variance(self):
        """Returns sample variance of counts in categories in self.
        
        Fault-tolerant: returns 0 if 0 or 1 items.
        """
        if not self:
            return None
        total = self.Sum
        count = self.Count
        if count <= 1: #no variance for a single item
            variance = 0.0
        else:
            variance = (self.SumSquares-(total*total)/count)/(count-1)
        return variance
    Variance = property(_get_variance)

    def _get_standard_deviation(self):
        """Returns sample standard deviation of items in self."""
        if not self:
            return None
        return sqrt(abs(self.Variance))
    StandardDeviation = property(_get_standard_deviation)

    def _get_mean(self):
        """Returns mean of items in self."""
        if not self:
            return None
        try:
            mean = self.Sum/self.Count
        except (ZeroDivisionError, FloatingPointError):
            mean = 0.0
        return mean
    Mean = property(_get_mean)

    def _get_uncertainty(self):
        """Returns the uncertainty of the Freqs.
        
        Calculates the Shannon uncertainty, defined as the sum of weighted
        log probabilities (multiplied by -1).
        
        Usage: H = self.uncertainty()

        """
        normalized = self.copy()
        normalized.normalize()
        total = 0
        for prob in normalized.values():
            if prob:
                total -= prob * log2(prob)
        return total
    Uncertainty = property(_get_uncertainty)

    def _get_mode(self):
        """Returns the most frequent item. If a tie, picks one at random.
        
        Usage: most_frequent = self.mode()

        """
        best = None
        best_count = 0
        for item, count in self.items():
            if count > best_count:
                best_count = count
                best = item
        return best
    Mode = property(_get_mode)


    def summarize(self):
        """Returns summary statistics for self."""
        return SummaryStatistics(Count=self.Count, Sum=self.Sum, \
            SumSquares=self.SumSquares, Variance=self.Variance)

    def getSortedList(self, descending=True, by_val=True):
        """Returns sorted list of tuples.

        descending: whether to sort highest to lowest (default True).
        by_val: whether to sort by val instead of key (default True).
        """
        if by_val:
            items = [(v, (k,v)) for k, v in self.items()]
        else:
            items = self.items()
        items.sort()
        if descending:
            items.reverse()
        if by_val:
            return [i[1] for i in items]
        else:
            return items

class UnsafeFreqs(FreqsI, dict):
    """Holds a frequency distribution, i.e. a set of categpory -> count pairs.

    Note: does not perform any validation. Use Freqs if data consistency
    is more important than speed.
    """
    pass

    def copy(self):
        """Returns copy of data in self, preserving class."""
        result = self.__class__()
        result.update(self)
        return result
 
def freqwatcher(x):
    """Checks frequencies are correct type and >= 0."""
    try:
        x = float(x)
    except:
        raise ConstraintError, "Could not convert frequency %s to float." % x
    if x >= 0:
        return x
    else:
        raise ConstraintError, "Got frequency %s < 0." % x

class Freqs(FreqsI, MappedDict):
    """Holds a frequency distribution, i.e. a set of category -> count pairs.
    
    Class data:
        ValueMask: function that transforms values before they are entered.
        RequiredKeys: keys that are automatically added with frequency 0 before
        frequencies are added.

    Performs (expensive) validation on many operations that change the 
    dictionary. Use UnsafeFreqs if speed is more important than validation.
    """
    ValueMask = FunctionWrapper(freqwatcher)
 
    def __init__(self, data=None, Constraint=None, Mask=None, ValueMask=None):
        """Passes on to superclass, but adds required keys if absent.
        
        Parameters (for polymorphism with MappedDict superclass):

        data:           data to load into self
        Constraint:     only items that Constraint __contains__ are allowed
        Mask:           function applied to keys before lookup
        ValueMask:      function applied to values before addition
        """
        super(Freqs, self).__init__(Constraint=Constraint, Mask=Mask, \
                ValueMask=ValueMask)
        self += data
        for key in self.RequiredKeys:
            if key not in self:
                self[key] = 0.0

class NumberFreqsI(FreqsI):
    """Interface for frequency distribution where values and counts are numbers.

    NOTE: In NumberFreqs (as opposed to Freqs), keys and values are assumed
    to be one axis. In other words, the data {1:5, 2:10} in Freqs is assumed
    to mean 5 items in category 1 and 10 items in category 2, for a mean
    of 7.5 items per category. In NumberFreqs, the same data would mean 5
    counts of 1 and 10 counts of 2, for a mean of (5*1 + 10*2)/15 = 1.66.
    """

    def isValid(self):
        """Returns True if all keys and values are numbers."""
        for item in self.items:
            for i in item:
                if not (isinstance(i, int) or isinstance(i, float)):
                    return False
        return True

    def _get_count(self):
        """Calculates sum of frequencies in the frequency distribution (i.e.
        number of occurrences). Useful for other stats functions.
        
        Usage: count = self.count()
        """
        return sum(self.values())
    Count = property(_get_count)

    def _get_sum(self):
        """Returns sum of items in self."""
        if not self:
            return None
        return sum([item*frequency for item,frequency in self.items()])
    Sum = property(_get_sum)

    def _get_sum_squares(self):
        """Returns sum of squares of items in self."""
        if not self:
            return None
        return sum([i*i*count for i, count in self.items()])
    SumSquares = property(_get_sum_squares)

    def _get_uncertainty(self):
        """Returns the uncertainty of the NumberFreqs.
        
        Calculates the Shannon uncertainty, defined as the sum of weighted
        log probabilities (multiplied by -1).

        Handled by conversion to Freqs, since numbers treated as categories.
        
        Usage: H = self.uncertainty()
        """
        f = UnsafeFreqs()
        f += self
        return f.Uncertainty
        
    Uncertainty = property(_get_uncertainty)
    
    def quantile(self, quantile):
        """Returns the specified quantile.
        
        Uses method type 7 from R."""
        def value_at_expanded_index(values, counts, index):
            cumsum = counts.cumsum()
            for i in range(cumsum.shape[0]):
                if cumsum[i] > index:
                    return values[i]
            
        vals = sorted(self.keys())
        counts = array([self[val] for val in vals])
        index = quantile * (counts.sum()-1)
        lo = int(floor(index))
        hi = int(ceil(index))
        diff = index - lo
        lo_val = value_at_expanded_index(vals, counts, lo)
        if diff != 0:
            hi_val = value_at_expanded_index(vals, counts, hi)
        else:
            hi_val = 0
        stat = (1-diff) * lo_val + diff * hi_val
        return stat
    
    def _get_median(self):
        """returns the median"""
        return self.quantile(0.5)
    
    Median = property(_get_median)


class UnsafeNumberFreqs(NumberFreqsI, dict):
    """Class holding freqs where keys and values are assumed to be numbers.

    Changes calculation of mean, standard deviation, etc. by assuming that
    the keys have weight proportional to their values (i.e. if the key is
    5 and the value is 3, it contributes 15 'units' rather than 3 to things
    like mean() and normalize()).

    Does not perform validation to check whether the keys and values are
    valid (will raise various exceptions depending on circumstances).
    """
    RequiredKeys = None
    
    def copy(self):
        """Returns copy of data in self, preserving class."""
        result = self.__class__()
        result.update(self)
        return result
 
class NumberFreqs(NumberFreqsI, MappedDict):
    """Class holding freqs where both keys and values are numbers.

    Mean, variance etc. assume that the data are frequencies of other
    numbers rather than treating each key as a separate category.
    
    Changes calculation of mean, standard deviation, etc. by assuming that
    the keys have weight proportional to their values (i.e. if the key is
    5 and the value is 3, it contributes 15 'units' rather than 3 to things
    like mean() and normalize()).

    Performs (expensive) validation to ensure that keys are floats and
    values are non-negative floats.
    
    All keys and values are automatically converted to float.
    """
    RequiredKeys = None
    Mask = FunctionWrapper(float)
    ValueMask = FunctionWrapper(freqwatcher)
 
    def __init__(self, data=None, Constraint=None, Mask=None, ValueMask=None):
        """Passes on to superclass, but adds required keys if absent.
        
        Parameters (for polymorphism with MappedDict superclass):

        data:           data to load into self
        Constraint:     only items that Constraint __contains__ are allowed
        Mask:           function applied to keys before lookup
        ValueMask:      function applied to values before addition
        """
        super(NumberFreqs, self).__init__(Constraint=Constraint, Mask=Mask, \
                ValueMask=ValueMask)
        self += data
        r = self.RequiredKeys
        if r:
            for key in r:
                if key not in self:
                    self[key] = 0.0


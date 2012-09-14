#!/usr/bin/env python
"""Generally useful utility classes and methods.
"""

import types
import sys
from time import clock
from datetime import datetime
from string import maketrans, strip
from random import randrange, choice, randint
from sys import maxint
from os import popen, remove, makedirs, getenv
from os.path import join, abspath, exists, isdir
from numpy import logical_not, sum
from cPickle import dumps, loads
from gzip import GzipFile
import hashlib
# import parse_command_line_parameters for backward compatibility
from cogent.util.option_parsing import parse_command_line_parameters

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Amanda Birmingham",
                    "Sandra Smit", "Zongzhi Liu", "Daniel McDonald",
                    "Kyle Bittinger", "Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def safe_md5(open_file, block_size=2**20):
    """Computes an md5 sum without loading the file into memory

    This method is based on the answers given in:
    http://stackoverflow.com/questions/1131220/get-md5-hash-of-a-files-without-open-it-in-python
    """
    md5 = hashlib.md5()
    data = True
    while data:
        data = open_file.read(block_size)
        if data:
            md5.update(data)
    return md5

def identity(x):
    """Identity function: useful for avoiding special handling for None."""
    return x

def if_(test, true_result, false_result):
    """Convenience function for one-line if/then/else with known return values.

    Note that both true_result and false_result are evaluated, which is not the
    case for the normal if/then/else, so don't use if one branch might fail.

    Additionally, true_result and false_result must be expressions, not 
    statements (so print and raise will not work, for example).
    """
    if test:
        return true_result
    else:
        return false_result

def iterable(item):
    """If item is iterable, returns item. Otherwise, returns [item].

    Useful for guaranteeing a result that can be iterated over.
    """
    try:
        iter(item)
        return item
    except TypeError:
        return [item]

def max_index(items):
    """Returns index of the largest item.

    items can be any sequence. If there is a tie, returns latest item.
    """
    return max([(item, index) for index, item in enumerate(items)])[1]

def min_index(items):
    """Returns index of the smallest item.

    items can be any sequence. If there is a tie, returns earliest item"""
    return min([(item, index) for index, item in enumerate(items)])[1]

def flatten(items):
    """Removes one level of nesting from items.

    items can be any sequence, but flatten always returns a list.
    """
    result = []
    for i in items:
        try:
            result.extend(i)
        except:
            result.append(i)
    return result

class DepthExceededError(Exception):
    pass

def deep_list(x):
    """Convert tuple recursively to list."""
    if isinstance(x, tuple):
        return map(deep_list, x)
    return x

def deep_tuple(x):
    """Convert list recursively to tuple."""
    if isinstance(x, list):
        return tuple(map(deep_tuple, x))
    return x

def between((min_, max_), number):
    """Same as: min_ <= number <= max_."""
    return min_ <= number <= max_

def combinate(items, n):
    """Returns combinations of items."""
    if n == 0: yield []
    else:
        for i in xrange(len(items) - n + 1):
            for cc in combinate(items[i + 1:], n - 1):
                yield [items[i]] + cc

def gzip_dump(object, filename, bin=2):
    """Saves a compressed object to file."""
    file = GzipFile(filename, 'wb')
    file.write(dumps(object, bin))
    try: # do not leave unlinked structures
        object.link()
    except AttributeError:
        pass
    file.close()

def gzip_load(filename):
    """Loads a compressed object from file."""
    file = GzipFile(filename, 'rb')
    buffer = []
    while True:
        data = file.read()
        if data == "":
            break
        buffer.append(data)
    buffer = "".join(buffer)
    object = loads(buffer)
    file.close()
    return object

def recursive_flatten_old(items, max_depth=None, curr_depth=0):
    """Removes all nesting from items, recursively.
    
    Note: Default max_depth is None, which removes all nesting (including
    unpacking strings). Setting max_depth unpacks a maximum of max_depth levels
    of nesting, but will not raise exception if the structure is not really
    that deep (instead, will just remove the nesting that exists). If max_depth
    is 0, will not remove any nesting (note difference from setting max_depth
    to None).
    """
    #bail out if greater than max_depth
    if max_depth is not None:
        if curr_depth > max_depth:
            raise DepthExceededError
    result = []
    for i in items:
        try:
            result.extend(recursive_flatten(i, max_depth, curr_depth + 1))
        except:
            result.append(i)
    return result

def curry(f, *a, **kw):
    """curry(f,x)(y) = f(x,y) or =lambda y: f(x,y)
    
    modified from python cookbook"""
    def curried(*more_a, **more_kw):
        return f(*(a + more_a), **dict(kw, **more_kw))

    ## make docstring for curried funtion
    curry_params = []
    if a:
        curry_params.extend([e for e in a])
    if kw:
        curry_params.extend(['%s=%s' % (k, v) for k, v in kw.items()])
    #str it to prevent error in join()
    curry_params = map(str, curry_params)

    try:
        f_name = f.func_name
    except:  #e.g.  itertools.groupby failed .func_name 
        f_name = '?'

    curried.__doc__ = ' curry(%s,%s)\n'\
        '== curried from %s ==\n %s'\
        % (f_name, ', '.join(curry_params), f_name, f.__doc__)

    return curried
#end curry

def is_iterable(obj):
    """return True if obj is iterable"""
    try:
        iter(obj)
    except TypeError, e:
        return False
    else:
        return True

def is_char(obj):
    """return True if obj is a char (str with lenth<=1)"""
    return isinstance(obj, basestring) and len(obj) <= 1

is_char_or_noniterable = lambda x: is_char(x) or\
        not is_iterable(x)

is_str_or_noniterable = lambda x: isinstance(x, basestring) or\
        not is_iterable(x)


def recursive_flatten(items, max_depth=None, curr_depth=1,
        is_leaf=is_char_or_noniterable):
    """Removes all nesting from items, recursively.

    Note: Default max_depth is None, which removes all nesting (including
    unpacking strings). Setting max_depth unpacks a maximum of max_depth levels
    of nesting, but will not raise exception if the structure is not really
    that deep (instead, will just remove the nesting that exists). If max_depth
    is 0, will not remove any nesting (note difference from setting max_depth
    to None).

    is_leaf: a predicate for 'leaf node'.  The default is_char_or_noniterable
    removes all nesting. is_str_or_noniterable removes all nesting sequences
    except strings. is_leaf=not_list_tuple removes only nesting list or tuple
    , which is considerably faster and recommended for general use.
    """
    result = []
    for i in items:
        if max_depth is not None and curr_depth > max_depth \
                or is_leaf(i):
            result.append(i)
        else:
            result.extend(recursive_flatten(i,
                    max_depth, curr_depth + 1, is_leaf))
    return result
#end recursive_flatten

def not_list_tuple(obj):
    """return False if obj is a list or a tuple"""
    return not isinstance(obj, (list, tuple))

list_flatten = curry(recursive_flatten, is_leaf=not_list_tuple)

def unflatten(data, row_width, keep_extras=False):
    """Converts items in data into a list of row_width-length lists.

    row_width must be an integer. Will raise error if zero.
    
    data can be any sequence type, but results will always be lists at the
    first level (including the common case of a list containing one sequence).

    Any items left over after the last complete row will be discarded. This
    means that the number of items in data need not be divisible by row_width.

    This function does _not_ reverse the effect of zip, since the lists it
    produces are not interleaved. If the list of lists were treated as a 2-d
    array, its transpose would be the reverse of the effect of zip (i.e. the
    original lists would be columns, not rows).
    """
    if row_width < 1:
        raise ValueError, "unflatten: row_width must be at least 1."
    result = []
    num_items = len(data)
    slices = num_items / row_width
    for s in xrange(slices):
        result.append(data[s * row_width:(s + 1) * row_width])
    if keep_extras:
        last_slice = data[slices * row_width:]
        if last_slice:
            result.append(last_slice)
    return result

def unzip(items):
    """Performs the reverse of zip, i.e. produces separate lists from tuples.

    items should be list of k-element tuples. Will raise exception if any tuples
    contain more items than the first one.

    Conceptually equivalent to transposing the matrix of tuples.

    Returns list of lists in which the ith list contains the ith element of each
    tuple.

    Note: zip expects *items rather than items, such that unzip(zip(*items))
    returns something that compares equal to items.

    Always returns lists: does not check original data type, but will accept
    any sequence.
    """
    if items:
        return map(list, zip(*items))
    else:
        return []

def select(order, items):
    """Returns the elements from items specified in order, a list of indices.

    Builds up a list containing the ith element in items for each item in order,
    which must be a list of valid keys into items.

    Example: vowels = select([0, 4, 8], 'abcdefghijklm')

    Can also be used to emulate Perl's hash slices.

    Example: reverse_vowel_freqs = select('ea', {'a':1,'b':5,'c':2,'d':4,'e':6})

    Return type is a list of whatever type the elements in items are.
    """
    return map(items.__getitem__, order)

def sort_order(items, cmpfunc=None):
    """Returns an array containing the sorted order of elements in items.

    The returned array contains indexes. Looking up items[i] for i in
    indexes returns a sorted list of the items in ascending order.

    Useful for returning just the n best items, etc.
    """
    indexed = [(item, index) for index, item in enumerate(items)]
    if cmpfunc is None:
        indexed.sort()
    else:
        indexed.sort(cmpfunc)
    return  [i[1] for i in indexed]

def find_all(text, pat):
    """Returns list of all overlapping occurrences of a pattern in a text.
    
    Each item in the (sorted) list is the index of one of the matches.
    """
    result = []
    last = 0
    try:
        while 1:
            curr = text.index(pat, last)
            result.append(curr)
            last = curr + 1
    except ValueError:  #raised when no more matches
        return result

def find_many(text, pats):
    """Returns sorted list of all occurrences of all patterns in text.

    Matches to all patterns are sorted together. Each item in the list is
    the index of one of the matches.

    WARNING: if pat is a single string, will search for the chars in the string
    individually. If this is not what you want (i.e. you want to search for the
    whole string), you should be using find_all instead; if you want to use
    find_many anyway, put the string in a 1-item list.
    """

    result = []
    for pat in pats:
        result.extend(find_all(text, pat))
    result.sort()
    return result



def unreserve(item):
    """Removes trailing underscore from item if it has one.

    Useful for fixing mutations of Python reserved words, e.g. class.
    """
    if hasattr(item, 'endswith') and item.endswith('_'):
        return item[:-1]
    else:
        return item

def add_lowercase(d):
    """Adds lowercase version of keys in d to itself. Converts vals as well.
    
    Should work on sequences of strings as well as strings.

    Now also works on strings and sets.
    """
    if hasattr(d, 'lower'):     #behaves like a string
        return d + d.lower()
    elif not hasattr(d, 'items'):   #not a dict
        items = list(d)
        return d.__class__(items + [i.lower() for i in items])

    #otherwise, assume dict-like behavior
    for key, val in d.items():
        try:
            new_key = key.lower()
        except:   #try to make tuple out of arbitrary sequence
            try:
                new_key = []
                for k in key:
                    try:
                        new_key.append(k.lower())
                    except:
                        new_key.append(k)
                new_key = tuple(new_key)
            except:
                new_key = key
        try:
            new_val = val.lower()
        except:
            new_val = val   #don't care if we couldn't convert it
        if new_key not in d:    #don't overwrite existing lcase keys
            d[new_key] = new_val
    return d

def extract_delimited(line, left, right, start_index=0):
    """Returns the part of line from first left to first right delimiter.
    
    Optionally begins searching at start_index.

    Note: finds the next complete field (i.e. if we start in an incomplete
    field, skip it and move to the next).
    """
    if left == right:
        raise TypeError, \
        "extract_delimited is for fields w/ different left and right delimiters"
    try:
        field_start = line.index(left, start_index)
    except ValueError:  #no such field
        return None
    else:
        try:
            field_end = line.index(right, field_start)
        except ValueError:  #left but no right delimiter: raise error
            raise ValueError, \
            "Found '%s' but not '%s' in line %s, starting at %s." \
            % (left, right, line, start_index)
    #if we got here, we found the start and end of the field
    return line[field_start + 1:field_end]

def caps_from_underscores(string):
    """Converts words_with_underscores to CapWords."""
    words = string.split('_')
    return ''.join([w.title() for w in words])

def InverseDict(d):
    """Returns inverse of d, setting keys to values and vice versa.

    Note: if more than one key has the same value, returns an arbitrary key
    for that value (overwrites with the last one encountered in the iteration).

    Can be invoked with anything that can be an argument for dict(), including
    an existing dict or a list of tuples. However, keys are always inserted in
    arbitrary order rather than input order.

    WARNING: will fail if any values are unhashable, e.g. if they are dicts or
    lists.
    """
    if isinstance(d, dict):
        temp = d
    else:
        temp = dict(d)
    return dict([(val, key) for key, val in temp.iteritems()])

def InverseDictMulti(d):
    """Returns inverse of d, setting keys to values and values to list of keys.

    Note that each value will _always_ be a list, even if one item.

    Can be invoked with anything that can be an argument for dict(), including
    an existing dict or a list of tuples. However, keys are always appended in
    arbitrary order, not the input order.

    WARNING: will fail if any values are unhashable, e.g. if they are dicts or
    lists.
    """
    if isinstance(d, dict):
        temp = d
    else:
        temp = dict(d)
    result = {}
    for key, val in temp.iteritems():
        if val not in result:
            result[val] = []
        result[val].append(key)
    return result

def DictFromPos(seq):
    """Returns dict mapping items to list of positions of those items in seq.

    Always assigns values as lists, even if the item appeared only once.

    WARNING: will fail if items in seq are unhashable, e.g. if seq is a list of
    lists.
    """
    result = {}
    for i, s in enumerate(seq):
        if s not in result:
            result[s] = []
        result[s].append(i)
    return result

def DictFromFirst(seq):
    """Returns dict mapping each item to index of its first occurrence in seq.

    WARNING: will fail if items in seq are unhashable, e.g. if seq is a list of
    lists.
    """
    result = {}
    for i, s in enumerate(seq):
        if s not in result:
            result[s] = i
    return result

def DictFromLast(seq):
    """Returns dict mapping each item to index of its last occurrence in seq.

    WARNING: will fail if items in seq are unhashable, e.g. if seq is a list of
    lists.
    """
    return dict([(item, index) for index, item in enumerate(seq)])

def DistanceFromMatrix(matrix):
    """Returns function(i,j) that looks up matrix[i][j].

    Useful for maintaining flexibility about whether a function is computed
    or looked up.

    Matrix can be a 2D dict (arbitrary keys) or list (integer keys).
    """
    def result(i, j):
        return matrix[i][j]
    return result

def PairsFromGroups(groups):
    """Returns dict such that d[(i,j)] exists iff i and j share a group.

    groups must be a sequence of sequences, e.g a list of strings.
    """
    result = {}
    for group in groups:
        for i in group:
            for j in group:
                result[(i, j)] = None
    return result

class ClassChecker(object):
    """Container for classes: 'if t in x == True' if t is the right class."""
    def __init__(self, *Classes):
        """Returns a new ClassChecker that accepts specified classes."""
        type_type = type(str)
        for c in Classes:
            if type(c) != type_type:
                raise TypeError, \
                "ClassChecker found non-type object '%s' in parameter list." % c
        self.Classes = list(Classes)

    def __contains__(self, item):
        """Returns True if item is a subclass of one of the classes in self."""
        for c in self.Classes:
            if isinstance(item, c):
                return True
        return False

    def __str__(self):
        """Informal string representation: returns list"""
        return str(self.Classes)

class Delegator(object):
    """Mixin class that forwards unknown attributes to a specified object.
    
    Handles properties correctly (this was somewhat subtle).

    WARNING: If you are delegating to an object that pretends to have every
    attribute (e.g. a MappedRecord), you _must_ bypass normal attribute access
    in __init__ of your subclasses to ensure that the properties are set in
    the object itself, not in the object to which it delegates. Alternatively,
    you can initialize with None so that unhandled attributes are set in self, 
    and then replace self._handler with your object right at the end of 
    __init__. The first option is probably safer and more general.

    Warning: will not work on classes that use __slots__ instead of __dict__.
    """
    def __init__(self, obj):
        """Returns a new Delegator that uses methods of obj.
        
        NOTE: It's important that this bypasses the normal attribute setting
        mechanism, or there's an infinite loop between __init__ and 
        __setattr__. However, subclasses should be able to use the normal
        mechanism with impunity.
        """
        self.__dict__['_handler'] = obj

    def __getattr__(self, attr):
        """Forwards unhandled attributes to self._handler.
        
        Sets _handler to None on first use if not already set.
        """
        handler = self.__dict__.setdefault('_handler', None)
        return getattr(handler, attr)

    def __setattr__(self, attr, value):
        """Forwards requests to change unhandled attributes to self._handler.
        
        This logic is rather complicated because of GenericRecord objects, which
        masquerade as having every attribute, which can be used as handlers for
        Delegators, which forward all requests to their handlers.

        Consequently, we need to check the following:

            1. Is attr in the object's __dict__? If so, set it in self.
            2. Does the handler have attr? If so, try to set it in handler.
            3. Does self lack the attr? If so, try to set it in handler.
            4. Did setting attr in the handler fail? If so, set it in self.
        """
        #if we're setting _handler, set it in dict directly (but complain if
        #it's self).
        if attr == '_handler':
            if value is self:
                raise ValueError, "Can't set object to be its own handler."
            self.__dict__['_handler'] = value
            return
        #check if the attribute is in this object's dict
        elif attr in self.__dict__:
            return object.__setattr__(self, attr, value)
        #then check if the class knows about it
        elif hasattr(self.__class__, attr):
            return object.__setattr__(self, attr, value)
        #then try to set it in the handler
        if hasattr(self._handler, attr) or not hasattr(self, attr):
            try:
                return setattr(self._handler, attr, value)
            except AttributeError:
                pass #will try to create the attribute on self
        return object.__setattr__(self, attr, value)

class FunctionWrapper(object):
    """Wraps a function to hide it from a class so that it isn't a method."""
    def __init__(self, Function):
        self.Function = Function
    def __call__(self, *args, **kwargs):
        return self.Function(*args, **kwargs)

class ConstraintError(Exception):
    """Raised when constraint on a container is violated."""
    pass

class ConstrainedContainer(object):
    """Mixin class providing constraint checking to a container.

    Container should have a Constraint property that __contains__ the items
    that will be allowed in the container. Can also have a Mask property that
    contains a function that will be applied to each item (a) on checking the
    item for validity, and (b) on inserting the item in the container. 

    WARNING: Because the Mask is evaluated both when the item is checked and 
    when it is inserted, any side-effects it has are applied _twice_. This 
    means that any Mask that mutates the object or changes global variables
    is unlikely to do what you want!
    """
    _constraint = None
    Mask = FunctionWrapper(identity)

    def _mask_for_new(self):
        """Returns self.Mask only if different from class data."""
        if self.Mask is not self.__class__.Mask:
            return self.Mask
        else:
            return None

    def __init__(self, Constraint=None, Mask=None):
        """Returns new ConstrainedContainer, incorporating constraint.
        
        WARNING: Does not perform validation. It is the subclass's 
        responsibility to perform validation during __init__ or __new__!
        """
        if Constraint is not None:
            self._constraint = Constraint
        if Mask is not None:
            self.Mask = Mask

    def matchesConstraint(self, constraint):
        """Returns True if all items in self are allowed."""
        #First checks if constraints are compatible. If not, or if the current
        #sequence has no constraint, does item by item search.

        #bail out if self or constraint is empty
        if not constraint or not self:
            return True
        #try checking constraints for compatibility
        if self.Constraint:
            try:
                constraint_ok = True
                for c in self.Constraint:
                    if c not in constraint:
                            constraint_ok = False
                            break
                if constraint_ok:
                    return True
            except TypeError:
                pass #e.g. tried to check wrong type item in string alphabet

        #get here if either self.Constraint is empty, or if we found an item
        #in self.Constraint that wasn't in the other constraint. In either case,
        #we need to check self item by item.
        if self:
            try:
                for i in self:
                    if i not in constraint:
                        return False
            except TypeError:   #e.g. tried to check int in string alphabet
                return False
        return True

    def otherIsValid(self, other):
        """Returns True if other has only items allowed in self.Constraint."""
        #First, checks other.Constrant for compatibility.
        #If other.Constraint is incompatible, checks items in other.
        mask = self.Mask
        constraint = self.Constraint
        if not constraint or not other:
            return True     #bail out if empty
        try:
            #if other has a constraint, check whether it's compatible
            other_constraint = other.Constraint
            if other_constraint:
                for c in map(mask, other_constraint):
                    if c not in constraint:
                        raise ConstraintError
                return True
        except (ConstraintError, AttributeError, TypeError):
            pass
        #get here if other doesn't have a constraint or if other's constraint
        #isn't valid on self's constraint.
        try:
            for item in map(mask, other):
                if item not in constraint:
                    return False
        except TypeError:
            return False    #e.g. tried to check int in str alphabet
        return True

    def itemIsValid(self, item):
        """Returns True if single item is in self.Constraint."""
        try:
            if (not self.Constraint) or self.Mask(item) in self.Constraint:
                return True
            else:
                return False
        except (TypeError, ConstraintError):  #wrong type or not allowed
            return False

    def sequenceIsValid(self, sequence):
        """Returns True if all items in sequence are in self.Constraint."""
        is_valid = self.itemIsValid
        for i in map(self.Mask, sequence):
            if not is_valid(i):
                return False
        return True

    def _get_constraint(self):
        """Accessor for constraint."""
        return self._constraint

    def _set_constraint(self, constraint):
        """Mutator for constraint."""
        if self.matchesConstraint(constraint):
            self._constraint = constraint
        else:
            raise ConstraintError, \
            "Sequence '%s' incompatible with constraint '%s'" % (self, constraint)

    Constraint = property(_get_constraint, _set_constraint)


class ConstrainedString(str, ConstrainedContainer):
    """String that is always valid on a specified constraint."""
    def __new__(cls, data, Constraint=None, Mask=None):
        """Constructor class method for validated ConstrainedString."""
        mask = Mask or cls.Mask
        if data == '':
            pass    #map can't handle an empty sequence, sadly...
        elif isinstance(data, str):
            data = ''.join(map(mask, data))
        else:
            try:
                data = str(map(mask, data))
            except (TypeError, ValueError):
                data = str(mask(data))
        new_string = str.__new__(cls, data)
        curr_constraint = Constraint or new_string.Constraint
        if curr_constraint and new_string:
            for c in new_string:
                try:
                    is_valid = c in curr_constraint
                except TypeError:
                    is_valid = False
                if not is_valid:
                    raise ConstraintError, \
                    "Character '%s' not in constraint '%s'" % (c, curr_constraint)
        return new_string

    def __init__(self, data, Constraint=None, Mask=None):
        """Constructor for validated ConstrainedString."""
        ConstrainedContainer.__init__(self, Constraint, Mask)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        if not self.otherIsValid(other):
            raise ConstraintError, \
            "Sequence '%s' doesn't meet constraint" % other
        result = self.__class__(str(self) + ''.join(map(self.Mask, other)), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__mul__(self, multiplier),
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__rmul__(self, multiplier),
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        result = self.__class__(str.__getslice__(self, *args, **kwargs), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __getitem__(self, index):
        """Make sure extended slice remembers the constraint."""
        items = str.__getitem__(self, index)
        if isinstance(index, slice):
            mask = self._mask_for_new()
            result = self.__class__(items, Constraint=self.Constraint)
            if mask:
                result.Mask = mask
            return result
        else:
            return items

class MappedString(ConstrainedString):
    """As for ConstrainedString, but maps __contained__ and __getitem__."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedString, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False


class ConstrainedList(ConstrainedContainer, list):
    """List that is always valid on a specified constraint."""

    def __init__(self, data=None, Constraint=None, Mask=None):
        """Constructor for validated ConstrainedList."""
        ConstrainedContainer.__init__(self, Constraint, Mask)
        if data:
            self.extend(data)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        result = self.__class__(list(self) + map(self.Mask, other) , \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __iadd__(self, other):
        """Adds other to self if constraint correct."""
        other = map(self.Mask, other)
        if self.otherIsValid(other):
            return list.__iadd__(self, other)
        else:
            raise ConstraintError, \
            "Sequence '%s' has items not in constraint '%s'" \
                % (other, self.Constraint)

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __setitem__(self, index, item):
        """Sets self[index] to item if item in constraint. Handles slices"""
        if isinstance(index, slice):
            if not self.otherIsValid(item):
                raise ConstraintError, \
                "Sequence '%s' contains items not in constraint '%s'." % \
                (item, self.Constraint)
            item = map(self.Mask, item)
        else:
            if not self.itemIsValid(item):
                raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                    (item, self.Constraint)
            item = self.Mask(item)
        list.__setitem__(self, index, item)

    def append(self, item):
        """Appends item to self."""
        if not self.itemIsValid(item):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (item, self.Constraint)
        list.append(self, self.Mask(item))

    def extend(self, sequence):
        """Appends sequence to self."""
        if self.otherIsValid(sequence):
            list.extend(self, map(self.Mask, sequence))
        else:
            raise ConstraintError, "Some items in '%s' not in constraint '%s'"\
                % (sequence, self.Constraint)

    def insert(self, position, item):
        """Inserts item at position in self."""
        if not self.itemIsValid(item):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (item, self.Constraint)
        list.insert(self, position, self.Mask(item))

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        result = self.__class__(list.__getslice__(self, *args, **kwargs), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __setslice__(self, start, end, sequence):
        """Make sure invalid data can't get into slice."""
        if self.otherIsValid(sequence):
            list.__setslice__(self, start, end, map(self.Mask, sequence))
        else:
            raise ConstraintError, \
                "Sequence '%s' has items not in constraint '%s'"\
                % (sequence, self.Constraint)

class MappedList(ConstrainedList):
    """As for ConstrainedList, but maps items on contains and getitem."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedList, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False

class ConstrainedDict(ConstrainedContainer, dict):
    """Dict containing only keys that are valid on a specified constraint.
    
    Default behavior when fed a sequence is to store counts of the items in
    that sequence, which is not the standard dict interface (should raise a
    ValueError instead) but which is surprisingly useful in practice.
    """
    ValueMask = FunctionWrapper(identity)

    def _get_mask_and_valmask(self):
        """Helper method to check whether Mask and ValueMask were set."""
        if self.Mask is self.__class__.Mask:
            mask = None
        else:
            mask = self.Mask

        if self.ValueMask is self.__class__.ValueMask:
            valmask = None
        else:
            valmask = self.ValueMask
        return mask, valmask

    def __init__(self, data=None, Constraint=None, Mask=None, ValueMask=None):
        """Constructor for validated ConstrainedDict."""
        ConstrainedContainer.__init__(self, Constraint, Mask)
        if ValueMask is not None:
            self.ValueMask = ValueMask
        if data:
            try:
                self.update(data)
            except (ValueError, TypeError):
                for d in map(self.Mask, iterable(data)):
                    curr = self.get(d, 0)
                    self[d] = curr + 1

    def __setitem__(self, key, value):
        """Sets self[key] to value if value in constraint."""
        if not self.itemIsValid(key):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (key, self.Constraint)
        key, value = self.Mask(key), self.ValueMask(value)
        dict.__setitem__(self, key, value)

    def copy(self):
        """Should return copy of self, including constraint."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(self, Constraint=self.Constraint, Mask=mask,
                ValueMask=valmask)

    def fromkeys(self, keys, value=None):
        """Returns new dictionary with same constraint as self."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(dict.fromkeys(keys, value),
            Constraint=self.Constraint, Mask=mask, ValueMask=valmask)

    def setdefault(self, key, default=None):
        """Returns self[key], setting self[key]=default if absent."""
        key, default = self.Mask(key), self.ValueMask(default)
        if key not in self:
            self[key] = default
        return self[key]

    def update(self, other):
        """Updates self with items in other.
        
        Implementation note: currently uses __setitem__, so no need to apply
        masks in this method.
        """
        if not hasattr(other, 'keys'):
            other = dict(other)
        for key in other:
            self[key] = other[key]

class MappedDict(ConstrainedDict):
    """As for ConstrainedDict, but maps keys on contains and getitem."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedDict, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False

    def __getitem__(self, item):
        """Ensure that getitem applies the mask."""
        return super(MappedDict, self).__getitem__(self.Mask(item))

    def get(self, item, default=None):
        """Ensure that get applies the mask."""
        return super(MappedDict, self).get(self.Mask(item), default)

    def has_key(self, item):
        """Ensure that has_key applies the mask."""
        return super(MappedDict, self).has_key(self.Mask(item))

def getNewId(rand_f=randrange):
    """Creates a random 12-digit integer to be used as an id."""

    NUM_DIGITS = 12
    return ''.join(map(str, [rand_f(10) for i in range(NUM_DIGITS)]))
#end function getNewId

def generateCombinations(alphabet, combination_length):
    """Returns an array of strings: all combinations of a given length.
    
    alphabet: a sequence (string or list) type object containing the
        characters that can be used to make the combinations.
    combination_length: a long-castable value (integer) specifying the
        length of desired combinations.
    
    comb is used as an abbreviation of combinations throughout.
    """

    found_combs = []
    num_combs = 0
    try:
        alphabet_len = len(alphabet)
        combination_length = long(combination_length)
    except TypeError, ValueError: #conversion failed
        raise RuntimeError, "Bad parameter: alphabet must be of sequence " + \
                            "type and combination_length must be castable " + \
                            "to long."
    #end parameter conversion try/catch

    #the number of combs is alphabet length raised to the combination length
    if combination_length != 0:
        num_combs = pow(alphabet_len, combination_length)
    #end if

    for curr_comb_num in xrange(num_combs):
        curr_digit = 0
        curr_comb = [0] * combination_length

        while curr_comb_num:
            curr_comb[curr_digit] = curr_comb_num % alphabet_len
            curr_comb_num = curr_comb_num / alphabet_len
            curr_digit += 1
        #end while

        #now translate the list of digits into a list of characters
        real_comb = []
        for position in curr_comb: real_comb.append(alphabet[position])
        found_combs.append("".join(real_comb))
    #next combination number

    return found_combs
#end generateCombinations

def toString(obj):
    """Public function to write a string of object's properties & their vals.
    
    This function looks only at the local properties/methods/etc of the
    object it is sent, and only examines public and first-level private
    (starts with _ but not __) entries.  It ignores anything that is a
    method, function, or class.  Any attribute whose value looks like a 
    printout of a memory address (starts with < and ends with >) has its 
    value replaced with the word "object".
    """

    ignored_types = [types.BuiltinFunctionType, types.BuiltinMethodType, \
                    types.ClassType, types.FunctionType, types.MethodType, \
                    types.UnboundMethodType]
    result = []
    for slot in obj.__dict__:
        if not slot.startswith("__"):
            ignore_attr = False
            attr = getattr(obj, slot)
            for ignored_type in ignored_types:
                if isinstance(attr, ignored_type): ignore_attr = True
            #next ignored type

            if not ignore_attr:
                attr_value = str(attr)
                if attr_value.startswith("<") and attr_value.endswith(">"):
                    attr_value = "object"
                #end if
                result.append(slot + ": " + attr_value)
            #end if
        #end if
    #next property

    return "; ".join(result)
#end toString

#A class for exceptions caused when param cannot be cast to nonneg int
class NonnegIntError(ValueError):
    """for exceptions caused when param cannot be cast to nonneg int"""

    def __init__(self, args=None):
        self.args = args
    #end __init__
#end NonnegIntError

def makeNonnegInt(n):
    """Public function to cast input to nonneg int and return, or raise err"""

    try:
        n = abs(int(n))
    except:
        raise NonnegIntError, n + " must be castable to a nonnegative int"
    #end try/except

    return n
#end makeNonnegInt

def reverse_complement(seq, use_DNA=True):
    """Public function to reverse complement DNA or RNA sequence string
    
    seq: a string
    use_DNA: a boolean indicating (if true) that A should translate to T.  
        If false, RNA is assumed (A translates to U).  Default is True.
    
    Returns a reverse complemented string.
    """
    bad_chars = set(seq) - set("ACGTUacgtu")
    if len(bad_chars) > 0:
        raise ValueError,\
         ("Only ACGTU characters may be passed to reverse_complement. Other "
          "characters were identified: %s. Use cogent.DNA.rc if you need to "
          "reverse complement ambiguous bases." % ''.join(bad_chars))
    #decide which translation to use for complementation
    if use_DNA:
        trans_table = maketrans("ACGTacgt", "TGCAtgca")
    else:
        trans_table = maketrans("ACGUacgu", "UGCAugca")
    #end if

    #complement the input sequence, then reverse
    complemented = seq.translate(trans_table)
    comp_list = list(complemented)
    comp_list.reverse()

    #join the reverse-complemented list and return
    return "".join(comp_list)
# The reverse_complement function was previously called revComp, but that 
# naming doesn't adhere to the PyCogent coding guidelines. Renamed, but 
# keeping the old name around to not break existing code.
revComp = reverse_complement
#end revComp

def timeLimitReached(start_time, time_limit):
    """Return true if more that time_limit has elapsed since start_time"""

    result = False
    curr_time = clock()
    elapsed = curr_time - start_time
    if elapsed > time_limit: result = True
    return result
#end _time_limit_reached

def not_none(seq):
    """Returns True if no item in seq is None."""
    for i in seq:
        if i is None:
            return False
    return True
#end not_none

def get_items_except(seq, indices, seq_constructor=None):
    """Returns all items in seq that are not in indices

    Returns the same type as parsed in except when a seq_constructor is set.
    """
    sequence = list(seq)
    index_lookup = dict.fromkeys(indices)
    result = [sequence[i] for i in range(len(seq)) \
                if i not in index_lookup]
    if not seq_constructor:
        if isinstance(seq, str):
            return ''.join(result)
        else:
            seq_constructor = seq.__class__
    return seq_constructor(result)
#end get_items_except

def NestedSplitter(delimiters=[None], same_level=False,
        constructor=strip, filter_=False):
    """return a splitter which return a list (maybe nested) from a str using 
    delimiters nestedly

    same_level -- if true, all the leaf items will be split whether there is 
    delimiters in it or not

    constructor: modify each splited fields.
    filter_: filter the splits if not False(default)
    
    Note: the line input in parser is expected to be a str, but without check
    """
    def parser(line, index=0):
        #split line with curr delimiter
        curr = delimiters[index]
        if isinstance(curr, (list, tuple)):
            try:
                delim, maxsplits = curr
            except ValueError:
                raise ValueError("delimiter tuple/list should be \
                        [delimiter_str, maxsplits]")
            if maxsplits < 0:
                result = line.rsplit(delim, -maxsplits)
            else:
                result = line.split(delim, maxsplits)
        else:
            result = line.split(curr)

        #modify splits if required
        if constructor:
            result = map(constructor, result)
        if filter_ != False: #allow filter(None,..) to rip off the empty items
            result = filter(filter_, result)

        #repeat recursively for next delimiter
        if index != len(delimiters) - 1: #not last delimiter
            result = [parser(f, index + 1) for f in result]

        #undo split if curr not in line and same_level==False
        #ignore the first delimiter
        if not same_level and index > 0 \
            and len(result) == 1 and isinstance(result[0], basestring):
            result = result[0]

        return result
    #parser.__doc__ = make_innerdoc(NestedSplitter, parser, locals())
    return parser
#end NestedSplitter

def app_path(app,env_variable='PATH'):
    """Returns path to an app, or False if app does not exist in env_variable
    
     This functions in the same way as which in that it returns
     the first path that contains the app.
    
    """
    # strip off " characters, in case we got a FilePath object
    app = app.strip('"')
    paths = getenv(env_variable).split(':')
    for path in paths:
        p = join(path,app)
        if exists(p):
            return p
    return False

#some error codes for creating a dir
def get_create_dir_error_codes():
    return {'NO_ERROR':      0,
            'DIR_EXISTS':    1,
            'FILE_EXISTS':   2,
            'OTHER_OS_ERROR':3}

def create_dir(dir_name, fail_on_exist=False, handle_errors_externally=False):
    """Create a dir safely and fail meaningful.

    dir_name: name of directory to create

    fail_on_exist: if true raise an error if dir already exists
    
    handle_errors_externally: if True do not raise Errors, but return
                   failure codes. This allows to handle errors locally and
                   e.g. hint the user at a --force_overwrite options.
                   
    returns values (if no Error raised):
    
         0:  dir could be safely made
         1:  directory already existed
         2:  a file with the same name exists          
         3:  any other unspecified OSError


    See qiime/denoiser.py for an example of how to use this mechanism.

    Note: Depending  of how thorough we want to be we could add tests,
          e.g. for testing actual write permission in an existing dir.
    """
    error_code_lookup = get_create_dir_error_codes()
    #pre-instanciate function with
    ror = curry(handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            #dir is there
            if fail_on_exist:
                return ror(error_code_lookup['DIR_EXISTS'])
            else:
                return error_code_lookup['DIR_EXISTS']
        else:
            #must be file with same name
            return ror(error_code_lookup['FILE_EXISTS'])
    else:
        #no dir there, try making it
        try:
            makedirs(dir_name)
        except OSError:
            return ror(error_code_lookup['OTHER_OS_ERROR'])
    
    return error_code_lookup['NO_ERROR']

def handle_error_codes(dir_name, supress_errors=False,
                       error_code=None):
    """Wrapper function for error_handling.

    dir_name: name of directory that raised the error
    suppress_errors: if True raise Errors, otherwise return error_code
    error_code: the code for the error
    """
    error_code_lookup = get_create_dir_error_codes()
    
    if error_code == None:
        error_code = error_code_lookup['NO_ERROR']
    
    error_strings = \
        {error_code_lookup['DIR_EXISTS'] :
          "Directory already exists: %s" % dir_name,
         error_code_lookup['FILE_EXISTS'] : 
          "File with same name exists: %s" % dir_name,
         error_code_lookup['OTHER_OS_ERROR']: 
          "Could not create output directory: %s. " % dir_name +
          "Check the permissions."}

    if error_code == error_code_lookup['NO_ERROR']:
        return error_code_lookup['NO_ERROR']
    if supress_errors:
        return error_code
    else:
        raise OSError, error_strings[error_code]

def remove_files(list_of_filepaths, error_on_missing=True):
    """Remove list of filepaths, optionally raising an error if any are missing
    """
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError, "Some filepaths were not accessible: %s" % '\t'.join(missing)

def get_random_directory_name(suppress_mkdir=False,\
    timestamp_pattern='%Y%m%d%H%M%S',\
    rand_length=20,\
    output_dir=None,\
    prefix='',
    suffix='',
    return_absolute_path=True):
    """Build a random directory name and create the directory 
    
        suppress_mkdir: only build the directory name, don't
         create the directory (default: False)
        timestamp_pattern: string passed to strftime() to generate
         the timestamp (pass '' to suppress the timestamp)
        rand_length: length of random string of characters
        output_dir: the directory which should contain the 
         random directory
        prefix: prefix for directory name
        suffix: suffix for directory name
        return_absolute_path: If False, returns the local (relative) path to the new directory
    """
    output_dir = output_dir or './'
    # Define a set of characters to be used in the random directory name
    chars = "abcdefghigklmnopqrstuvwxyz"
    picks = chars + chars.upper() + "0123456790"
    
    # Get a time stamp
    timestamp = datetime.now().strftime(timestamp_pattern)
        
    # Construct the directory name
    dirname = '%s%s%s%s' % (prefix,timestamp,\
                        ''.join([choice(picks) for i in range(rand_length)]),\
                        suffix)
    dirpath = join(output_dir,dirname)
    abs_dirpath = abspath(dirpath)

    # Make the directory
    if not suppress_mkdir:
        try:
            makedirs(abs_dirpath)
        except OSError:
            raise OSError,\
             "Cannot make directory %s. Do you have write access?" % dirpath
             
    # Return the path to the directory
    if return_absolute_path:
        return abs_dirpath
    return dirpath

def get_independent_coords(spans, random_tie_breaker=False):
    """returns non-overlapping spans. spans must have structure
        [(start, end, ..), (..)]. spans can be decorated with arbitrary data
        after the end entry.
    
    Arguments:
        - random_tie_breaker: break overlaps by randomly choosing the first
          or second span. Defaults to the first span.
    """
    
    if len(spans) <= 1:
        return spans
    
    last = spans[0]
    result = [last]
    for i in range(1, len(spans)):
        curr = spans[i]
        if curr[0] < last[1]:
            if random_tie_breaker:
                result[-1] = [last, curr][randint(0,1)]
            else:
                result[-1] = last
            continue
        
        result.append(curr)
        last = curr
    
    return result

def get_merged_overlapping_coords(start_end):
    """merges overlapping spans, assumes sorted by start"""
    result = [start_end[0]]
    prev_end = result[0][-1]
    for i in range(1, len(start_end)):
        curr_start, curr_end = start_end[i]
        # if we're beyond previous, add and continue
        if curr_start > prev_end:
            prev_end = curr_end
            result.append([curr_start, curr_end])
        elif curr_end > prev_end:
            prev_end = curr_end
            result[-1][-1] = prev_end
        else:
            pass # we lie completely within previous span
    
    return result

def get_run_start_indices(values, digits=None, converter_func=None):
    """returns starting index, value for all distinct values"""
    assert not (digits and converter_func), \
        'Cannot set both digits and converter_func'
    
    if digits is not None:
        converter_func = lambda x: round(x, digits)
    elif converter_func is None:
        converter_func = lambda x: x
    
    last_val = None
    for index, val in enumerate(values):
        val = converter_func(val)
        if val != last_val:
            yield [index, val]
        
        last_val = val
    
    return

def get_merged_by_value_coords(spans_value, digits=None):
    """returns adjacent spans merged if they have the same value. Assumes
    [(start, end, val), ..] structure and that spans_value is sorted in
    ascending order.
    
    Arguments:
        - digits: if None, any data can be handled and exact values are
          compared. Otherwise values are rounded to that many digits.
    """
    assert len(spans_value[0]) == 3, 'spans_value must have 3 records per row'
    
    starts, ends, vals = zip(*spans_value)
    indices_distinct_vals = get_run_start_indices(vals, digits=digits)
    data = []
    i = 0
    for index, val in indices_distinct_vals:
        start = starts[index]
        end = ends[index]
        prev_index = max(index-1, 0)
        try:
            data[-1][1] = ends[prev_index]
        except IndexError:
            pass
        
        data.append([start, end, val])
    
    if index < len(ends):
        data[-1][1] = ends[-1]
    
    return data

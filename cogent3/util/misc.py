#!/usr/bin/env python
"""Generally useful utility classes and methods.
"""

import types
import sys
from time import clock
from datetime import datetime
from random import randrange, choice, randint
from sys import maxsize
from os import popen, remove, makedirs, getenv
from os.path import join, abspath, exists, isdir
from tempfile import gettempdir
from numpy import logical_not, sum, array, float64, finfo
from pickle import dumps, loads
from gzip import GzipFile, open as gzip_open
from bz2 import open as bzip_open
import hashlib

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Amanda Birmingham",
               "Sandra Smit", "Zongzhi Liu", "Daniel McDonald",
               "Kyle Bittinger", "Marcin Cieslik"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


def adjusted_gt_minprob(probs, minprob=1e-6):
    """returns numpy array of probs scaled such that minimum is > minval
    
    result sums to 1 within machine precision"""
    assert 0 <= minprob < 1, "invalid minval %s" % minprob
    probs = array(probs, dtype=float64)
    if (probs > minprob).all():
        return probs
    
    total = probs.sum()
    smallest = probs.min()
    dim = probs.shape[0]
    # we need an adjustment that (smalvall+adj)/(adj + total) > minval
    # the following solves for this, then adds machine precision
    adj = -(smallest +  minprob * total) / (minprob * dim - 1)
    adj += finfo(float64).eps
    
    probs += adj
    probs /= probs.sum()
    return probs

def bytes_to_string(data):
    """returns a string if data is bytes, otherwise returns original"""
    if isinstance(data, bytes):
        data = data.decode('utf_8')
    return data


def open_(filename, mode='rt', **kwargs):
    """open that handles different compression"""
    op = {'gz': gzip_open, 'bz2': bzip_open}.get(
        filename.split('.')[-1], open)
    return op(filename, mode, **kwargs)


def get_format_suffixes(filename):
    """returns compression and/or file suffixes"""
    if '.' not in filename:
        return None, None

    compression_suffixes = ("bz2", "gz")
    filename = filename.split(".")
    suffixes = filename[1:] if len(filename) == 2 else filename[-2:]
    if suffixes[-1] in compression_suffixes:
        cmp_suffix = suffixes.pop(-1)
    else:
        cmp_suffix = None

    if suffixes:
        suffix = suffixes[-1]
    else:
        suffix = None
    return suffix, cmp_suffix


class FilePath(str):
    """ Hold paths for proper handling

        Paths in this sense are filenames, directory paths, or filepaths.
        Some examples include:
         file.txt
         ./path/to/file.txt
         ./path/to/dir/
         /path/to/file.txt
         .
         /

        The purpose of this class is to allow all paths to be handled the
         same since they sometimes need to be treated differently than
         simple strings. For example, if a path has a space in it, and it
         is being passed to system, it needs to be wrapped in quotes. But,
         you wouldn't want it as a string wrapped in quotes b/c, e.g.,
         isabs('"/absolute/path"') == False, b/c the first char is a ", not
         a /.

        * This would make more sense to call Path, but that conflicts with
            the ResultPath.Path attribute. I'm not sure what to do about this
            and want to see what others think. Once finalized, a global
            replace should take care of making the switch.

    """
    def __new__(cls, path):
        try:
            return str.__new__(cls, path.strip('"'))
        except AttributeError:
            return str.__new__(cls, '')

    def __str__(self):
        """ wrap self in quotes, or return the empty string if self == '' """
        if self == '':
            return ''
        return ''.join(['"', self, '"'])

    def __add__(self, other):
        return FilePath(''.join([self, other]))


def get_tmp_filename(tmp_dir=gettempdir(), prefix="tmp", suffix=".txt",
                     result_constructor=FilePath):
    """ Generate a temporary filename and return as a FilePath object

        tmp_dir: the directory to house the tmp_filename (default: '/tmp')
        prefix: string to append to beginning of filename (default: 'tmp')
            Note: It is very useful to have prefix be descriptive of the
            process which is creating the temporary file. For example, if
            your temp file will be used to build a temporary blast database,
            you might pass prefix=TempBlastDB
        suffix: the suffix to be appended to the temp filename (default '.txt')
        result_constructor: the constructor used to build the result filename
            (default: cogent3.app.parameters.FilePath). Note that joining
            FilePath objects with one another or with strings, you must use
            the + operator. If this causes trouble, you can pass str as the
            the result_constructor.
    """
    # check not none
    if not tmp_dir:
        tmp_dir = ""
    # if not current directory, append "/" if not already on path
    elif not tmp_dir.endswith("/"):
        tmp_dir += "/"

    chars = "abcdefghigklmnopqrstuvwxyz"
    picks = chars + chars.upper() + "0123456790"
    return result_constructor(tmp_dir) + result_constructor(prefix) +\
        result_constructor("%s%s" %
                           (''.join([choice(picks) for i in range(20)]), suffix))


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


def curry(f, *a, **kw):
    """curry(f,x)(y) = f(x,y) or = lambda y: f(x,y)

    modified from python cookbook"""
    def curried(*more_a, **more_kw):
        return f(*(a + more_a), **dict(kw, **more_kw))

    # make docstring for curried funtion
    curry_params = []
    if a:
        curry_params.extend([e for e in a])
    if kw:
        curry_params.extend(['%s=%s' % (k, v) for k, v in list(kw.items())])
    # str it to prevent error in join()
    curry_params = list(map(str, curry_params))

    try:
        f_name = f.__name__
    except:  # e.g.  itertools.groupby failed .func_name
        f_name = '?'

    curried.__doc__ = ' curry(%s,%s)\n'\
        '== curried from %s ==\n %s'\
        % (f_name, ', '.join(curry_params), f_name, f.__doc__)

    return curried
# end curry


def is_iterable(obj):
    """return True if obj is iterable"""
    try:
        iter(obj)
    except TypeError as e:
        return False
    else:
        return True


def is_char(obj):
    """return True if obj is a char (str with lenth<=1)"""
    return isinstance(obj, str) and len(obj) <= 1


def is_char_or_noniterable(x): return is_char(x) or\
    not is_iterable(x)


def is_str_or_noniterable(x): return isinstance(x, str) or\
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
# end recursive_flatten


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
        raise ValueError("unflatten: row_width must be at least 1.")
    result = []
    num_items = len(data)
    slices = num_items // row_width
    for s in range(slices):
        result.append(data[s * row_width:(s + 1) * row_width])
    if keep_extras:
        last_slice = data[slices * row_width:]
        if last_slice:
            result.append(last_slice)
    return result


def select(order, items):
    """Returns the elements from items specified in order, a list of indices.

    Builds up a list containing the ith element in items for each item in order,
    which must be a list of valid keys into items.

    Example: vowels = select([0, 4, 8], 'abcdefghijklm')

    Can also be used to emulate Perl's hash slices.

    Example: reverse_vowel_freqs = select('ea', {'a':1,'b':5,'c':2,'d':4,'e':6})

    Return type is a list of whatever type the elements in items are.
    """
    return list(map(items.__getitem__, order))


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
    except ValueError:  # raised when no more matches
        return result


def add_lowercase(d):
    """Adds lowercase version of keys in d to itself. Converts vals as well.

    Should work on sequences of strings as well as strings.

    Now also works on strings and sets.
    """
    if hasattr(d, 'lower'):  # behaves like a string
        return d + d.lower()
    elif not hasattr(d, 'items'):  # not a dict
        items = list(d)
        return d.__class__(items + [i.lower() for i in items])

    # otherwise, assume dict-like behavior
    for key, val in list(d.items()):
        try:
            new_key = key.lower()
        except:  # try to make tuple out of arbitrary sequence
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
            new_val = val  # don't care if we couldn't convert it
        if new_key not in d:  # don't overwrite existing lcase keys
            d[new_key] = new_val
    return d


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
    return dict([(val, key) for key, val in temp.items()])


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
    for key, val in temp.items():
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
                raise TypeError(
                    "ClassChecker found non-type object '%s' in parameter list." % c)
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
        # if we're setting _handler, set it in dict directly (but complain if
        # it's self).
        if attr == '_handler':
            if value is self:
                raise ValueError("Can't set object to be its own handler.")
            self.__dict__['_handler'] = value
            return
        # check if the attribute is in this object's dict
        elif attr in self.__dict__:
            return object.__setattr__(self, attr, value)
        # then check if the class knows about it
        elif hasattr(self.__class__, attr):
            return object.__setattr__(self, attr, value)
        # then try to set it in the handler
        if hasattr(self._handler, attr) or not hasattr(self, attr):
            try:
                return setattr(self._handler, attr, value)
            except AttributeError:
                pass  # will try to create the attribute on self
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

    Container should have a constraint property that __contains__ the items
    that will be allowed in the container. Can also have a mask property that
    contains a function that will be applied to each item (a) on checking the
    item for validity, and (b) on inserting the item in the container.

    WARNING: Because the mask is evaluated both when the item is checked and
    when it is inserted, any side-effects it has are applied _twice_. This
    means that any mask that mutates the object or changes global variables
    is unlikely to do what you want!
    """
    _constraint = None
    mask = FunctionWrapper(identity)

    def _mask_for_new(self):
        """Returns self.mask only if different from class data."""
        if self.mask is not self.__class__.mask:
            return self.mask
        else:
            return None

    def __init__(self, constraint=None, mask=None):
        """Returns new ConstrainedContainer, incorporating constraint.

        WARNING: Does not perform validation. It is the subclass's
        responsibility to perform validation during __init__ or __new__!
        """
        if constraint is not None:
            self._constraint = constraint
        if mask is not None:
            self.mask = mask

    def matches_constraint(self, constraint):
        """Returns True if all items in self are allowed."""
        # First checks if constraints are compatible. If not, or if the current
        # sequence has no constraint, does item by item search.

        # bail out if self or constraint is empty
        if not constraint or not self:
            return True
        # try checking constraints for compatibility
        if self.constraint:
            try:
                constraint_ok = True
                for c in self.constraint:
                    if c not in constraint:
                        constraint_ok = False
                        break
                if constraint_ok:
                    return True
            except TypeError:
                pass  # e.g. tried to check wrong type item in string alphabet

        # get here if either self.constraint is empty, or if we found an item
        # in self.constraint that wasn't in the other constraint. In either case,
        # we need to check self item by item.
        if self:
            try:
                for i in self:
                    if i not in constraint:
                        return False
            except TypeError:  # e.g. tried to check int in string alphabet
                return False
        return True

    def other_is_valid(self, other):
        """Returns True if other has only items allowed in self.constraint."""
        # First, checks other.Constrant for compatibility.
        # If other.constraint is incompatible, checks items in other.
        mask = self.mask
        constraint = self.constraint
        if not constraint or not other:
            return True  # bail out if empty
        try:
            # if other has a constraint, check whether it's compatible
            other_constraint = other.constraint
            if other_constraint:
                for c in map(mask, other_constraint):
                    if c not in constraint:
                        raise ConstraintError
                return True
        except (ConstraintError, AttributeError, TypeError):
            pass
        # get here if other doesn't have a constraint or if other's constraint
        # isn't valid on self's constraint.
        try:
            for item in map(mask, other):
                if item not in constraint:
                    return False
        except TypeError:
            return False  # e.g. tried to check int in str alphabet
        return True

    def item_is_valid(self, item):
        """Returns True if single item is in self.constraint."""
        try:
            if (not self.constraint) or self.mask(item) in self.constraint:
                return True
            else:
                return False
        except (TypeError, ConstraintError):  # wrong type or not allowed
            return False

    def sequence_is_valid(self, sequence):
        """Returns True if all items in sequence are in self.constraint."""
        is_valid = self.item_is_valid
        for i in map(self.mask, sequence):
            if not is_valid(i):
                return False
        return True

    def _get_constraint(self):
        """Accessor for constraint."""
        return self._constraint

    def _set_constraint(self, constraint):
        """Mutator for constraint."""
        if self.matches_constraint(constraint):
            self._constraint = constraint
        else:
            raise ConstraintError(
                "Sequence '%s' incompatible with constraint '%s'" % (self, constraint))

    constraint = property(_get_constraint, _set_constraint)


class ConstrainedString(str, ConstrainedContainer):
    """String that is always valid on a specified constraint."""
    def __new__(cls, data, constraint=None, mask=None):
        """Constructor class method for validated ConstrainedString."""
        mask = mask or cls.mask
        if data == '':
            pass  # map can't handle an empty sequence, sadly...
        elif isinstance(data, str):
            data = ''.join(map(mask, data))
        else:
            try:
                data = str(list(map(mask, data)))
            except (TypeError, ValueError):
                data = str(mask(data))
        new_string = str.__new__(cls, data)
        curr_constraint = constraint or new_string.constraint
        if curr_constraint and new_string:
            for c in new_string:
                try:
                    is_valid = c in curr_constraint
                except TypeError:
                    is_valid = False
                if not is_valid:
                    raise ConstraintError(
                        "Character '%s' not in constraint '%s'" % (c, curr_constraint))
        return new_string

    def __init__(self, data, constraint=None, mask=None):
        """Constructor for validated ConstrainedString."""
        ConstrainedContainer.__init__(self, constraint, mask)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        if not self.other_is_valid(other):
            raise ConstraintError(
                "Sequence '%s' doesn't meet constraint" % other)
        result = self.__class__(str(self) + ''.join(map(self.mask, other)),
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__mul__(self, multiplier),
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__rmul__(self, multiplier),
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        result = self.__class__(str.__getslice__(self, *args, **kwargs),
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __getitem__(self, index):
        """Make sure extended slice remembers the constraint."""
        items = str.__getitem__(self, index)
        if isinstance(index, slice):
            mask = self._mask_for_new()
            result = self.__class__(items, constraint=self.constraint)
            if mask:
                result.mask = mask
            return result
        else:
            return items


class MappedString(ConstrainedString):
    """As for ConstrainedString, but maps __contained__ and __getitem__."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedString, self).__contains__(self.mask(item))
        except (TypeError, ValueError):
            return False


class ConstrainedList(ConstrainedContainer, list):
    """List that is always valid on a specified constraint."""

    def __init__(self, data=None, constraint=None, mask=None):
        """Constructor for validated ConstrainedList."""
        ConstrainedContainer.__init__(self, constraint, mask)
        if data:
            self.extend(data)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        result = self.__class__(list(self) + list(map(self.mask, other)),
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __iadd__(self, other):
        """Adds other to self if constraint correct."""
        other = list(map(self.mask, other))
        if self.other_is_valid(other):
            return list.__iadd__(self, other)
        else:
            raise ConstraintError("Sequence '%s' has items not in constraint '%s'"
                                  % (other, self.constraint))

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
                                constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __setitem__(self, index, item):
        """Sets self[index] to item if item in constraint. Handles slices"""
        if isinstance(index, slice):
            if not self.other_is_valid(item):
                raise ConstraintError("Sequence '%s' contains items not in constraint '%s'." %
                                      (item, self.constraint))
            item = list(map(self.mask, item))
        else:
            if not self.item_is_valid(item):
                raise ConstraintError("Item '%s' not in constraint '%s'" %
                                      (item, self.constraint))
            item = self.mask(item)
        list.__setitem__(self, index, item)

    def __setslice__(self, start, end, sequence):
        """Make sure invalid data can't get into slice."""
        if self.other_is_valid(sequence):
            list.__setslice__(self, start, end, list(map(self.mask, sequence)))
        else:
            raise ConstraintError("Sequence '%s' has items not in constraint '%s'"
                                  % (sequence, self.constraint))

    def append(self, item):
        """Appends item to self."""
        if not self.item_is_valid(item):
            raise ConstraintError("Item '%s' not in constraint '%s'" %
                                  (item, self.constraint))
        list.append(self, self.mask(item))

    def extend(self, sequence):
        """Appends sequence to self."""
        if self.other_is_valid(sequence):
            list.extend(self, list(map(self.mask, sequence)))
        else:
            raise ConstraintError("Some items in '%s' not in constraint '%s'"
                                  % (sequence, self.constraint))

    def insert(self, position, item):
        """Inserts item at position in self."""
        if not self.item_is_valid(item):
            raise ConstraintError("Item '%s' not in constraint '%s'" %
                                  (item, self.constraint))
        list.insert(self, position, self.mask(item))

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        # to be deleted in py3
        val = list.__getslice__(self, *args, **kwargs)
        result = self.__class__(val, constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __getitem__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        if len(args) == 1 and type(args[0]) == int and not kwargs:
            val = list.__getitem__(self, args[0])
            return val

        val = list.__getitem__(self, *args, **kwargs)
        result = self.__class__(val, constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result


class MappedList(ConstrainedList):
    """As for ConstrainedList, but maps items on contains and getitem."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedList, self).__contains__(self.mask(item))
        except (TypeError, ValueError):
            return False


class ConstrainedDict(ConstrainedContainer, dict):
    """Dict containing only keys that are valid on a specified constraint.

    Default behavior when fed a sequence is to store counts of the items in
    that sequence, which is not the standard dict interface (should raise a
    ValueError instead) but which is surprisingly useful in practice.
    """
    value_mask = FunctionWrapper(identity)

    def _get_mask_and_valmask(self):
        """Helper method to check whether mask and value_mask were set."""
        if self.mask is self.__class__.mask:
            mask = None
        else:
            mask = self.mask

        if self.value_mask is self.__class__.value_mask:
            valmask = None
        else:
            valmask = self.value_mask
        return mask, valmask

    def __init__(self, data=None, constraint=None, mask=None, value_mask=None):
        """Constructor for validated ConstrainedDict."""
        ConstrainedContainer.__init__(self, constraint, mask)
        if value_mask is not None:
            self.value_mask = value_mask
        if data:
            try:
                self.update(data)
            except (ValueError, TypeError):
                for d in map(self.mask, iterable(data)):
                    curr = self.get(d, 0)
                    self[d] = curr + 1

    def __setitem__(self, key, value):
        """Sets self[key] to value if value in constraint."""
        if not self.item_is_valid(key):
            raise ConstraintError("Item '%s' not in constraint '%s'" %
                                  (key, self.constraint))
        key, value = self.mask(key), self.value_mask(value)
        dict.__setitem__(self, key, value)

    def copy(self):
        """Should return copy of self, including constraint."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(self, constraint=self.constraint, mask=mask,
                              value_mask=valmask)

    def fromkeys(self, keys, value=None):
        """Returns new dictionary with same constraint as self."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(dict.fromkeys(keys, value),
                              constraint=self.constraint, mask=mask, value_mask=valmask)

    def setdefault(self, key, default=None):
        """Returns self[key], setting self[key]=default if absent."""
        key, default = self.mask(key), self.value_mask(default)
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
            return super(MappedDict, self).__contains__(self.mask(item))
        except (TypeError, ValueError):
            return False

    def __getitem__(self, item):
        """Ensure that getitem applies the mask."""
        return super(MappedDict, self).__getitem__(self.mask(item))

    def get(self, item, default=None):
        """Ensure that get applies the mask."""
        return super(MappedDict, self).get(self.mask(item), default)

    def has_key(self, item):
        """Ensure that has_key applies the mask."""
        return self.mask(item) in super(MappedDict, self)


def to_string(obj):
    """Public function to write a string of object's properties & their vals.

    This function looks only at the local properties/methods/etc of the
    object it is sent, and only examines public and first-level private
    (starts with _ but not __) entries.  It ignores anything that is a
    method, function, or class.  Any attribute whose value looks like a
    printout of a memory address (starts with < and ends with >) has its
    value replaced with the word "object".
    """

    ignored_types = [types.BuiltinFunctionType, types.BuiltinMethodType,
                     type, types.FunctionType, types.MethodType]
    result = []
    for slot in obj.__dict__:
        if not slot.startswith("__"):
            ignore_attr = False
            attr = getattr(obj, slot)
            for ignored_type in ignored_types:
                if isinstance(attr, ignored_type):
                    ignore_attr = True
            # next ignored type

            if not ignore_attr:
                attr_value = str(attr)
                if attr_value.startswith("<") and attr_value.endswith(">"):
                    attr_value = "object"
                # end if
                result.append(slot + ": " + attr_value)
            # end if
        # end if
    # next property

    return "; ".join(result)
# end to_string

# A class for exceptions caused when param cannot be cast to nonneg int


class NonnegIntError(ValueError):
    """for exceptions caused when param cannot be cast to nonneg int"""

    def __init__(self, args=None):
        self.args = args
    # end __init__
# end NonnegIntError


def reverse_complement(seq, use_DNA=True):
    """Public function to reverse complement DNA or RNA sequence string

    seq: a string
    use_DNA: a boolean indicating (if true) that A should translate to T.
        If false, RNA is assumed (A translates to U).  Default is True.

    Returns a reverse complemented string.
    """
    bad_chars = set(seq) - set("ACGTUacgtu")
    if len(bad_chars) > 0:
        raise ValueError("Only ACGTU characters may be passed to reverse_complement. Other "
                         "characters were identified: %s. Use cogent3.DNA.rc if you need to "
                         "reverse complement ambiguous bases." % ''.join(bad_chars))
    # decide which translation to use for complementation
    if use_DNA:
        trans_table = str.maketrans("ACGTacgt", "TGCAtgca")
    else:
        trans_table = str.maketrans("ACGUacgu", "UGCAugca")
    # end if

    # complement the input sequence, then reverse
    complemented = seq.translate(trans_table)
    comp_list = list(complemented)
    comp_list.reverse()

    # join the reverse-complemented list and return
    return "".join(comp_list)
# end revComp


def timeLimitReached(start_time, time_limit):
    """Return true if more that time_limit has elapsed since start_time"""

    result = False
    curr_time = clock()
    elapsed = curr_time - start_time
    if elapsed > time_limit:
        result = True
    return result
# end _time_limit_reached


def not_none(seq):
    """Returns True if no item in seq is None."""
    for i in seq:
        if i is None:
            return False
    return True
# end not_none


def NestedSplitter(delimiters=[None], same_level=False,
                   constructor=str.strip, filter_=False):
    """return a splitter which return a list (maybe nested) from a str using
    delimiters nestedly

    same_level -- if true, all the leaf items will be split whether there is
    delimiters in it or not

    constructor: modify each splited fields.
    filter_: filter the splits if not False(default)

    Note: the line input in parser is expected to be a str, but without check
    """
    def parser(line, index=0):
        # split line with curr delimiter
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

        # modify splits if required
        if constructor:
            result = list(map(constructor, result))
        # allow filter(None,..) to rip off the empty items
        if filter_ is not False:
            result = list(filter(filter_, result))

        # repeat recursively for next delimiter
        if index != len(delimiters) - 1:  # not last delimiter
            result = [parser(f, index + 1) for f in result]

        # undo split if curr not in line and same_level==False
        # ignore the first delimiter
        if not same_level and index > 0 \
                and len(result) == 1 and isinstance(result[0], str):
            result = result[0]

        return result
    # parser.__doc__ = make_innerdoc(NestedSplitter, parser, locals())
    return parser
# end NestedSplitter


def get_create_dir_error_codes():
    return {'NO_ERROR': 0,
            'DIR_EXISTS': 1,
            'FILE_EXISTS': 2,
            'OTHER_OS_ERROR': 3}


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
    # pre-instanciate function with
    ror = curry(handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            # dir is there
            if fail_on_exist:
                return ror(error_code_lookup['DIR_EXISTS'])
            else:
                return error_code_lookup['DIR_EXISTS']
        else:
            # must be file with same name
            return ror(error_code_lookup['FILE_EXISTS'])
    else:
        # no dir there, try making it
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

    if error_code is None:
        error_code = error_code_lookup['NO_ERROR']

    error_strings = \
        {error_code_lookup['DIR_EXISTS']:
         "Directory already exists: %s" % dir_name,
         error_code_lookup['FILE_EXISTS']:
         "File with same name exists: %s" % dir_name,
         error_code_lookup['OTHER_OS_ERROR']:
         "Could not create output directory: %s. " % dir_name +
         "Check the permissions."}

    if error_code == error_code_lookup['NO_ERROR']:
        return error_code_lookup['NO_ERROR']
    if supress_errors:
        return error_code
    else:
        raise OSError(error_strings[error_code])


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
        raise OSError("Some filepaths were not accessible: %s" %
                      '\t'.join(missing))


def get_random_directory_name(suppress_mkdir=False,
                              timestamp_pattern='%Y%m%d%H%M%S',
                              rand_length=20,
                              output_dir=None,
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
    dirname = '%s%s%s%s' % (prefix, timestamp,
                            ''.join([choice(picks)
                                     for i in range(rand_length)]),
                            suffix)
    dirpath = join(output_dir, dirname)
    abs_dirpath = abspath(dirpath)

    # Make the directory
    if not suppress_mkdir:
        try:
            makedirs(abs_dirpath)
        except OSError:
            raise OSError(
                "Cannot make directory %s. Do you have write access?" % dirpath)

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
                result[-1] = [last, curr][randint(0, 1)]
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
            pass  # we lie completely within previous span

    return result


def get_run_start_indices(values, digits=None, converter_func=None):
    """returns starting index, value for all distinct values"""
    assert not (digits and converter_func), \
        'Cannot set both digits and converter_func'

    if digits is not None:
        def converter_func(x): return round(x, digits)
    elif converter_func is None:
        def converter_func(x): return x

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

    starts, ends, vals = list(zip(*spans_value))
    indices_distinct_vals = get_run_start_indices(vals, digits=digits)
    data = []
    i = 0
    for index, val in indices_distinct_vals:
        start = starts[index]
        end = ends[index]
        prev_index = max(index - 1, 0)
        try:
            data[-1][1] = ends[prev_index]
        except IndexError:
            pass

        data.append([start, end, val])

    if index < len(ends):
        data[-1][1] = ends[-1]

    return data

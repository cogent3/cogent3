"""Generally useful utility classes and methods."""

from __future__ import annotations

import contextlib
import inspect
import os
import re
import typing
import warnings
from random import choice
from typing import Any
from urllib.parse import urlparse
from warnings import warn

import numpy
from numpy import array, finfo, float64, floating, integer, ndarray, zeros

if typing.TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Callable


def _adjusted_gt_minprob_vector(probs, minprob):
    # operates on a 1D numpy vector
    total = probs.sum()
    smallest = probs.min()
    if smallest > minprob:
        # nothing to do
        return probs

    dim = probs.shape[0]
    # we need an adjustment that (small_val + adj) / (n * adj + total) > minprob
    # the following solves for this, then adds machine precision
    adj = -(smallest + minprob * total) / (minprob * dim - 1)
    adj += finfo(float64).eps

    probs += adj
    probs /= probs.sum()
    return probs


def adjusted_gt_minprob(probs, minprob=1e-6):
    """returns numpy array of probs scaled such that minimum is > minval

    result sums to 1 within machine precision

    if 2D array, assumes row-order"""
    assert 0 <= minprob < 1, f"invalid minval {minprob}"
    probs = array(probs, dtype=float64)
    if (probs > minprob).all():
        return probs

    if probs.ndim == 1:
        probs = _adjusted_gt_minprob_vector(probs, minprob)
    else:
        for i in range(probs.shape[0]):
            probs[i] = _adjusted_gt_minprob_vector(probs[i], minprob)

    return probs


def adjusted_within_bounds(value, lower, upper, eps=1e-7, action="warn"):
    """returns value such that lower <= value <= upper

    Parameters
    ----------
    value
        number, converted to float64
    lower
        lower bound
    upper
        upper bound
    eps : float
        if value lies within eps of either lower/upper, it's returned inside
        this interval by machine precision
    action : str
        'warn', 'raise' (ValueError), 'ignore'. What happens if value lies further than eps
        from either bound
    """
    if lower <= value <= upper:
        return value

    assert action in ("warn", "raise", "ignore"), f"Unknown action {action!r}"

    value = float64(value)
    eps = float64(eps) + finfo(float64).eps
    err_msg = f"value[{value}] not within lower[{lower}]/upper[{upper}] bounds"
    wrn_msg = f"value[{value}] forced within lower[{lower}]/upper[{upper}] bounds"

    if value < lower and (lower - value) <= eps:
        value = lower
    elif value > upper and (value - upper) <= eps:
        value = upper
    elif (lower > value or value > upper) and action == "raise":
        raise ValueError(err_msg)
    else:
        warn(wrn_msg, category=UserWarning, stacklevel=2)
        value = upper if value > upper else lower

    return value


def bytes_to_string(data):
    """returns a string if data is bytes, otherwise returns original"""
    if isinstance(data, bytes):
        data = data.decode("utf_8")
    elif isinstance(data, list | tuple):
        data = "".join(str(d) for d in data)
    else:
        data = str(data)
    return data


_wout_period = re.compile(r"^\.")


def iterable(item):
    """If item is iterable, returns item. Otherwise, returns [item].

    Useful for guaranteeing a result that can be iterated over.
    """
    try:
        iter(item)
        return item
    except TypeError:
        return [item]


def curry(f, *a, **kw):
    """curry(f,x)(y) = f(x,y) or = lambda y: f(x,y)

    modified from python cookbook"""

    def curried(*more_a, **more_kw):
        return f(*(a + more_a), **dict(kw, **more_kw))

    # make docstring for curried funtion
    curry_params = []
    if a:
        curry_params.extend(list(a))
    if kw:
        curry_params.extend([f"{k}={v}" for k, v in list(kw.items())])
    # str it to prevent error in join()
    curry_params = list(map(str, curry_params))

    try:
        f_name = f.__name__
    except:  # e.g.  itertools.groupby failed .func_name
        f_name = "?"

    curried.__doc__ = " curry({},{})\n== curried from {} ==\n {}".format(
        f_name,
        ", ".join(curry_params),
        f_name,
        f.__doc__,
    )

    return curried


# end curry


def is_iterable(obj) -> bool:
    """return True if obj is iterable"""
    try:
        iter(obj)
    except TypeError:
        return False
    else:
        return True


def is_char(obj):
    """return True if obj is a char (str with lenth<=1)"""
    return isinstance(obj, str) and len(obj) <= 1


def is_char_or_noniterable(x):
    return is_char(x) or not is_iterable(x)


def is_url(text: str) -> bool:
    _urls = re.compile("^(http[s]*|file)")
    r = urlparse(text)
    return _urls.search(r.scheme) is not None


def recursive_flatten(
    items,
    max_depth=None,
    curr_depth=1,
    is_leaf=is_char_or_noniterable,
):
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
        if (max_depth is not None and curr_depth > max_depth) or is_leaf(i):
            result.append(i)
        else:
            result.extend(recursive_flatten(i, max_depth, curr_depth + 1, is_leaf))
    return result


def not_list_tuple(obj) -> bool:
    """return False if obj is a list or a tuple"""
    return not isinstance(obj, list | tuple)


list_flatten = curry(recursive_flatten, is_leaf=not_list_tuple)


def add_lowercase(d):
    """Adds lowercase version of keys in d to itself. Converts vals as well.

    Should work on sequences of strings as well as strings.

    Now also works on strings and sets.
    """
    if hasattr(d, "lower"):  # behaves like a string
        return d + d.lower()
    if not hasattr(d, "items"):  # not a dict
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


def DistanceFromMatrix(matrix):
    """Returns function(i,j) that looks up matrix[i][j].

    Useful for maintaining flexibility about whether a function is computed
    or looked up.

    Matrix can be a 2D dict (arbitrary keys) or list (integer keys).
    """

    def result(i, j):
        return matrix[i][j]

    return result


class ClassChecker:
    """Container for classes: 'if t in x == True' if t is the right class."""

    def __init__(self, *Classes) -> None:
        """Returns a new ClassChecker that accepts specified classes."""
        type_type = type(str)
        for c in Classes:
            if type(c) != type_type:
                msg = f"ClassChecker found non-type object '{c}' in parameter list."
                raise TypeError(
                    msg,
                )
        self.Classes = list(Classes)

    def __contains__(self, item) -> bool:
        """Returns True if item is a subclass of one of the classes in self."""
        return any(isinstance(item, c) for c in self.Classes)

    def __str__(self) -> str:
        """Informal string representation: returns list"""
        return str(self.Classes)


class Delegator:
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

    def __init__(self, obj) -> None:
        """Returns a new Delegator that uses methods of obj.

        NOTE: It's important that this bypasses the normal attribute setting
        mechanism, or there's an infinite loop between __init__ and
        __setattr__. However, subclasses should be able to use the normal
        mechanism with impunity.
        """
        self.__dict__["_handler"] = obj

    def __getattr__(self, attr):
        """Forwards unhandled attributes to self._handler.

        Sets _handler to None on first use if not already set.
        """
        handler = self.__dict__.setdefault("_handler", None)
        return getattr(handler, attr)

    def __setattr__(self, attr, value) -> None:
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
        if attr == "_handler":
            if value is self:
                msg = "Can't set object to be its own handler."
                raise ValueError(msg)
            self.__dict__["_handler"] = value
            return None
        # check if the attribute is in this object's dict
        if attr in self.__dict__ or hasattr(self.__class__, attr):
            return object.__setattr__(self, attr, value)
        # then try to set it in the handler
        if hasattr(self._handler, attr) or not hasattr(self, attr):
            try:
                return setattr(self._handler, attr, value)
            except AttributeError:
                pass  # will try to create the attribute on self
        return object.__setattr__(self, attr, value)


class FunctionWrapper:
    """Wraps a function to hide it from a class so that it isn't a method."""

    def __init__(self, Function=None) -> None:
        self.Function = Function if Function is not None else self._default

    def _default(self, x):
        return x

    def __call__(self, *args, **kwargs):
        return self.Function(*args, **kwargs)


class ConstraintError(Exception):
    """Raised when constraint on a container is violated."""


class ConstrainedContainer:
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
    mask = FunctionWrapper()

    def _mask_for_new(self):
        """Returns self.mask only if different from class data."""
        if self.mask is not self.__class__.mask:
            return self.mask
        return None

    def __init__(self, constraint=None, mask=None) -> None:
        """Returns new ConstrainedContainer, incorporating constraint.

        WARNING: Does not perform validation. It is the subclass's
        responsibility to perform validation during __init__ or __new__!
        """
        if constraint is not None:
            self._constraint = constraint
        if mask is not None:
            self.mask = mask

    def matches_constraint(self, constraint) -> bool:
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

    def other_is_valid(self, other) -> bool:
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

    def item_is_valid(self, item) -> bool | None:
        """Returns True if single item is in self.constraint."""
        try:
            return bool(not self.constraint or self.mask(item) in self.constraint)
        except (TypeError, ConstraintError):  # wrong type or not allowed
            return False

    def sequence_is_valid(self, sequence) -> bool:
        """Returns True if all items in sequence are in self.constraint."""
        is_valid = self.item_is_valid
        return all(is_valid(i) for i in map(self.mask, sequence))

    def _get_constraint(self):
        """Accessor for constraint."""
        return self._constraint

    def _set_constraint(self, constraint) -> None:
        """Mutator for constraint."""
        if self.matches_constraint(constraint):
            self._constraint = constraint
        else:
            msg = f"Sequence '{self}' incompatible with constraint '{constraint}'"
            raise ConstraintError(
                msg,
            )

    constraint = property(_get_constraint, _set_constraint)


class ConstrainedList(ConstrainedContainer, list):
    """List that is always valid on a specified constraint."""

    def __init__(self, data=None, constraint=None, mask=None) -> None:
        """Constructor for validated ConstrainedList."""
        ConstrainedContainer.__init__(self, constraint, mask)
        if data:
            self.extend(data)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        result = self.__class__(
            list(self) + list(map(self.mask, other)),
            constraint=self.constraint,
        )
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __iadd__(self, other):
        """Adds other to self if constraint correct."""
        other = list(map(self.mask, other))
        if self.other_is_valid(other):
            return list.__iadd__(self, other)
        msg = f"Sequence '{other}' has items not in constraint '{self.constraint}'"
        raise ConstraintError(
            msg,
        )

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier, constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier, constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result

    def __setitem__(self, index, item) -> None:
        """Sets self[index] to item if item in constraint. Handles slices"""
        if isinstance(index, slice):
            if not self.other_is_valid(item):
                msg = f"Sequence '{item}' contains items not in constraint '{self.constraint}'."
                raise ConstraintError(
                    msg,
                )
            item = list(map(self.mask, item))
        else:
            if not self.item_is_valid(item):
                msg = f"Item '{item}' not in constraint '{self.constraint}'"
                raise ConstraintError(
                    msg,
                )
            item = self.mask(item)
        list.__setitem__(self, index, item)

    def __setslice__(self, start, end, sequence):
        """Make sure invalid data can't get into slice."""
        if self.other_is_valid(sequence):
            list.__setslice__(self, start, end, list(map(self.mask, sequence)))
        else:
            msg = (
                f"Sequence '{sequence}' has items not in constraint '{self.constraint}'"
            )
            raise ConstraintError(
                msg,
            )

    def append(self, item) -> None:
        """Appends item to self."""
        if not self.item_is_valid(item):
            msg = f"Item '{item}' not in constraint '{self.constraint}'"
            raise ConstraintError(
                msg,
            )
        list.append(self, self.mask(item))

    def extend(self, sequence) -> None:
        """Appends sequence to self."""
        if self.other_is_valid(sequence):
            list.extend(self, list(map(self.mask, sequence)))
        else:
            msg = f"Some items in '{sequence}' not in constraint '{self.constraint}'"
            raise ConstraintError(
                msg,
            )

    def insert(self, position, item) -> None:
        """Inserts item at position in self."""
        if not self.item_is_valid(item):
            msg = f"Item '{item}' not in constraint '{self.constraint}'"
            raise ConstraintError(
                msg,
            )
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
            return list.__getitem__(self, args[0])

        val = list.__getitem__(self, *args, **kwargs)
        result = self.__class__(val, constraint=self.constraint)
        mask = self._mask_for_new()
        if mask:
            result.mask = mask
        return result


class MappedList(ConstrainedList):
    """As for ConstrainedList, but maps items on contains and getitem."""

    def __contains__(self, item) -> bool:
        """Ensure that contains applies the mask."""
        try:
            return super().__contains__(self.mask(item))
        except (TypeError, ValueError):
            return False


class ConstrainedDict(ConstrainedContainer, dict):
    """Dict containing only keys that are valid on a specified constraint.

    Default behavior when fed a sequence is to store counts of the items in
    that sequence, which is not the standard dict interface (should raise a
    ValueError instead) but which is surprisingly useful in practice.
    """

    value_mask = FunctionWrapper()

    def _get_mask_and_valmask(self):
        """Helper method to check whether mask and value_mask were set."""
        mask = None if self.mask is self.__class__.mask else self.mask

        if self.value_mask is self.__class__.value_mask:
            valmask = None
        else:
            valmask = self.value_mask
        return mask, valmask

    def __init__(self, data=None, constraint=None, mask=None, value_mask=None) -> None:
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

    def __setitem__(self, key, value) -> None:
        """Sets self[key] to value if value in constraint."""
        if not self.item_is_valid(key):
            msg = f"Item '{key}' not in constraint '{self.constraint}'"
            raise ConstraintError(msg)
        key, value = self.mask(key), self.value_mask(value)
        dict.__setitem__(self, key, value)

    def copy(self):
        """Should return copy of self, including constraint."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(
            self,
            constraint=self.constraint,
            mask=mask,
            value_mask=valmask,
        )

    def fromkeys(self, keys, value=None):
        """Returns new dictionary with same constraint as self."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(
            dict.fromkeys(keys, value),
            constraint=self.constraint,
            mask=mask,
            value_mask=valmask,
        )

    def setdefault(self, key, default=None):
        """Returns self[key], setting self[key]=default if absent."""
        key, default = self.mask(key), self.value_mask(default)
        if key not in self:
            self[key] = default
        return self[key]

    def update(self, other) -> None:
        """Updates self with items in other.

        Implementation note: currently uses __setitem__, so no need to apply
        masks in this method.
        """
        if not hasattr(other, "keys"):
            other = dict(other)
        for key in other:
            self[key] = other[key]


class MappedDict(ConstrainedDict):
    """As for ConstrainedDict, but maps keys on contains and getitem."""

    def __contains__(self, item) -> bool:
        """Ensure that contains applies the mask."""
        try:
            return super().__contains__(self.mask(item))
        except (TypeError, ValueError):
            return False

    def __getitem__(self, item):
        """Ensure that getitem applies the mask."""
        return super().__getitem__(self.mask(item))

    def get(self, item, default=None):
        """Ensure that get applies the mask."""
        return super().get(self.mask(item), default)

    def has_key(self, item):
        """Ensure that has_key applies the mask."""
        return self.mask(item) in super()


def NestedSplitter(
    delimiters=None,
    same_level=False,
    constructor=str.strip,
    filter_=False,
):
    """return a splitter which return a list (maybe nested) from a str using
    delimiters nestedly

    same_level -- if true, all the leaf items will be split whether there is
    delimiters in it or not

    constructor: modify each splited fields.
    filter_: filter the splits if not False(default)

    Note: the line input in parser is expected to be a str, but without check
    """

    delimiters = delimiters or [None]

    def parser(line, index=0):
        # split line with curr delimiter
        curr = delimiters[index]
        if isinstance(curr, list | tuple):
            try:
                delim, maxsplits = curr
            except ValueError:
                msg = "delimiter tuple/list should be \
                        [delimiter_str, maxsplits]"
                raise ValueError(
                    msg,
                )
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
        if (
            not same_level
            and index > 0
            and len(result) == 1
            and isinstance(result[0], str)
        ):
            result = result[0]

        return result

    # parser.__doc__ = make_innerdoc(NestedSplitter, parser, locals())
    return parser


def get_independent_coords(spans, random_tie_breaker=False):
    """returns non-overlapping spans. spans must have structure
        [(start, end, ..), (..)]. spans can be decorated with arbitrary data
        after the end entry.

    Parameters
    ----------
    random_tie_breaker
        break overlaps by randomly choosing the first
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
                result[-1] = choice((last, curr))
            else:
                result[-1] = last
            continue

        result.append(curr)
        last = curr

    return result


def get_merged_overlapping_coords(start_end):
    """merges overlapping spans, assumes sorted by start"""
    result = [list(start_end[0])]
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
    assert not (digits and converter_func), "Cannot set both digits and converter_func"

    if digits is not None:

        def converter_func(x):
            return round(x, digits)

    elif converter_func is None:

        def converter_func(x):
            return x

    last_val = None
    for index, val in enumerate(values):
        val = converter_func(val)
        if val != last_val:
            yield [index, val]

        last_val = val


def get_merged_by_value_coords(spans_value, digits=None):
    """returns adjacent spans merged if they have the same value. Assumes
    [(start, end, val), ..] structure and that spans_value is sorted in
    ascending order.

    Parameters
    ----------
    digits
        if None, any data can be handled and exact values are
        compared. Otherwise values are rounded to that many digits.

    """
    assert len(spans_value[0]) == 3, "spans_value must have 3 records per row"

    starts, ends, vals = list(zip(*spans_value, strict=False))
    indices_distinct_vals = get_run_start_indices(vals, digits=digits)
    data = []
    for index, val in indices_distinct_vals:
        start = starts[index]
        end = ends[index]
        prev_index = max(index - 1, 0)
        with contextlib.suppress(IndexError):
            data[-1][1] = ends[prev_index]

        data.append([start, end, val])

    if index < len(ends):
        data[-1][1] = ends[-1]

    return data


def get_object_provenance(obj):
    """returns string of complete object provenance"""
    # algorithm inspired by Greg Baacon's answer to
    # https://stackoverflow.com/questions/2020014/get-fully-qualified-class
    # -name-of-an-object-in-python
    if isinstance(obj, type) or inspect.isfunction(obj):
        mod = obj.__module__
        name = obj.__name__
    else:
        mod = obj.__class__.__module__
        name = obj.__class__.__name__

    return name if mod is None or mod == "builtins" else f"{mod}.{name}"


def extend_docstring_from(source, pre=False):
    def docstring_inheriting_decorator(dest):
        parts = [source.__doc__ or "", dest.__doc__ or ""]
        # trim leading/trailing blank lines from parts
        for i, part in enumerate(parts):
            part = part.split("\n")
            if not part[0].strip():
                part.pop(0)
            if part and not part[-1].strip():
                part.pop(-1)

            parts[i] = "\n".join(part)

        if pre:
            parts.reverse()
        dest.__doc__ = "\n".join(parts)
        return dest

    return docstring_inheriting_decorator


_doc_block = re.compile(
    r"^\s*(Parameters|Notes|Raises)",
    flags=re.IGNORECASE | re.MULTILINE,
)


def docstring_to_summary_rest(text: str) -> tuple[str, str]:
    """separates the summary at the start of a docstring from the rest

    Notes
    -----
    Assumes numpydoc style.
    """
    if not text:
        return "", ""

    pos = _doc_block.search(text)
    if pos is None:
        return text, ""

    summary = text[: pos.start()].rstrip()
    text = text[pos.start() :]
    return summary, text.lstrip("\n").rstrip(" ")


def ascontiguousarray(
    source_array: numpy.ndarray,
    dtype: numpy.dtype | None = None,
) -> numpy.ndarray:
    if source_array is not None:
        return numpy.ascontiguousarray(source_array, dtype=dtype)
    return source_array


def get_setting_from_environ(environ_var, params_types):
    """extract settings from environment variable

    Parameters
    ----------
    environ_var : str
        name of an environment variable
    params_types : dict
        {param name: type}, values will be cast to type

    Returns
    -------
    dict

    Notes
    -----
    settings must of form 'param_name1=param_val,param_name2=param_val2'
    """
    var = os.environ.get(environ_var, None)
    if var is None:
        return {}

    var = var.split(",")
    result = {}
    for item in var:
        item = item.split("=")
        if len(item) != 2 or item[0] not in params_types:
            continue

        name, val = item
        try:
            val = params_types[name](val)
            result[name] = val
        except Exception:
            warnings.warn(
                f"could not cast {name}={val} to type {params_types[name]}, skipping",
                stacklevel=2,
            )

    return result


def in_jupyter() -> bool:
    """whether code is being executed within a jupyter notebook"""
    val = True
    try:
        # primitive approach, just check whether the following function
        # is in the namespace
        get_ipython
    except NameError:
        val = False

    return val


def get_true_spans(arr: ndarray, absolute_pos: bool = True) -> ndarray:
    """
    returns array of [[start position, length],...] of True runs

    Parameters
    ----------
    arr
        bool array
    absolute_pos
        the absolute position in self, otherwise the position
        represents the index in the underlying series of False values
        between which the True run occurs.

    Notes
    -----
    Designed for use with individual records from Alignment.get_gap_array(). This
    method returns a bool array with the gap state as True.
    """
    result = zeros((arr.shape[0], 2), dtype=int)
    true_run = False
    num_runs = 0
    cum_sum = 0
    for i, v in enumerate(arr):
        if v and true_run:
            result[num_runs][1] += 1
            cum_sum += 1
        elif not v and true_run:
            # we just ended a run of False values
            num_runs += 1
            true_run = False
        elif v:
            result[num_runs][:] = [i if absolute_pos else i - cum_sum, 1]
            true_run = True
            cum_sum += 1

    if true_run:
        # finished in a run
        num_runs += 1

    return result[:num_runs]


def get_first_value(data: typing.Iterable[typing.Any]):
    """return first value from a container

    Parameters
    ----------
    data
        any iterable

    Returns
    -------
    Returns the first element of data. If data is a dict, it will return
    the first element of data.values().
    """
    with contextlib.suppress(StopIteration):
        if hasattr(data, "values"):
            return next(iter(data.values()))

        with contextlib.suppress(TypeError):
            return next(iter(data))


class negate_condition:
    """negates the result of calling func."""

    def __init__(self, func: Callable[[typing.Any], bool]) -> None:
        self.func = func

    def __call__(self, *args: Any, **kwargs: Any) -> bool:
        return not bool(self.func(*args, **kwargs))


def is_number(v) -> bool:
    """built in or numpy number type"""
    return isinstance(v, (int, float, integer, floating))


def is_int(v) -> bool:
    """built in or numpy int type"""
    return isinstance(v, (int, integer))


def is_float(v) -> bool:
    """built in or numpy float type"""
    return isinstance(v, (float, floating))

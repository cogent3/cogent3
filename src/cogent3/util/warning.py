#!/usr/bin/env python
import functools

from typing import Any, Callable, List, Tuple
from warnings import catch_warnings, simplefilter
from warnings import warn as _warn


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jai Ram Rideout", "Richard Morris"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def deprecated(_type, old, new, version, reason=None, stack_level=3):
    """a convenience function for deprecating classes, functions, arguments.

    Parameters
    ----------
    _type
        should be one of class, method, function, argument
    old, new
        the old and new names
    version
        the version by which support for the old name will be
        discontinued
    reason
        why, and what choices users have
    stack_level
        as per warnings.warn

    """
    msg = (
        f"use {_type} {new} instead of {old}, support discontinued in version {version}"
    )
    if reason is not None:
        msg = f"{msg}\n{reason}"

    with catch_warnings():
        simplefilter("always")
        _warn(msg, DeprecationWarning, stacklevel=stack_level)


def discontinued(_type, name, version, reason=None, stack_level=3):
    """convenience func to warn about discontinued attributes

    Parameters
    ----------
    _type
        should be one of class, method, function, argument
    name
        the attributes name
    version
        the version by which support for the old name will be
        discontinued
    reason
        why, and what choices users have
    stack_level
        as per warnings.warn
    """
    msg = (
        f"{_type} {name} is discontinued, support will be stopped in version {version}"
    )
    if reason is not None:
        msg = f"{msg}\n{reason}"

    with catch_warnings():
        simplefilter("always")
        _warn(msg, DeprecationWarning, stacklevel=stack_level)


def deprecated_args(
    mapping: List[Tuple[str, str]],
    version: str,
    reason: str,
) -> Callable[..., Any]:
    """
    A decorator that marks specific arguments of a function as deprecated.

    The decorator accepts a list of 2-tuples specifying the mapping of old argument names to new argument names.
    When the decorated function is called with any of the old argument names, they will be replaced with their
    corresponding new names in the kwargs dictionary.

    Parameters
    ----------
    mapping : List[Tuple[str, str]]
        A list of 2-tuples specifying the mapping of old argument names to new argument names.
    version : str
        A string indicating the version when the deprecated arguments will be removed.
    reason : str
        A string providing a reason for deprecation.

    Returns
    -------
    Callable[..., Any]
        The decorated function.

    Warnings
    --------
    DeprecationWarning
        A warning will be raised when the decorated function is called with any of the deprecated arguments.

    Examples
    --------
    Here's an example of how to use the `deprecated_args` decorator to mark the argument `old_name` as deprecated
    and replace it with the new name `new_name`.

    >>> @deprecated_args(mapping=[('old_name', 'new_name')], version='2.0', reason='Use new_name instead')
    >>> def my_function(new_name):
    >>>     # do something here

    When `my_function` is called with the argument `old_name`, a warning will be raised indicating that the argument
    is deprecated and should be replaced with `new_name`.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
            message = []
            for old, new in mapping:
                if old in kwargs:
                    kwargs[new] = kwargs.pop(old)
                    message.append(f"{old}->{new}")
            if message:
                msg = f'The following parameters of `{func.__name__}` are deprecated. [{", ".join(message)}] will be removed in {version or "a future release"}.{f"  {reason}" if reason else ""}'
                _warn(msg, DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)

        return wrapper

    return decorator

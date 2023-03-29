#!/usr/bin/env python
from typing import Any, Callable, Union
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


def deprecate(
    alternate: Union[Callable[..., Any], None] = None,
    version: str = None,
    reason: str = None,
) -> Callable[..., Any]:
    """
    Decorator to mark a function as deprecated and optionally to be replaced by an
    alternate function.

    Args:
        alternate (Callable[..., Any], optional): The alternate function to be called
            instead of the deprecated function. Defaults to None.
        version (str, optional): Specifies when the deprecated function will be removed.
            Defaults to None.
        reason (str, optional): The reason for deprecating the function. Defaults to None.

    The decorator will call the alternate function when the decorated function is called.
    If an alternate function is not provided, or is not a callable function
    then the original function will be called instead.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:

            deprecated(
                _type="function",
                old=func.__name__,
                new=alternate.__name__ if alternate else None,
                version=version,
                reason=reason,
            )

            return (
                alternate(*args, **kwargs)  # call the alternate function
                if alternate and callable(alternate)
                else func(*args, **kwargs)  # call the original function
            )

        return wrapper

    return decorator

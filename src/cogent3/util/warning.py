import functools
import inspect

from typing import Any, Callable, List, Optional, Tuple, Union
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
    msg = f"{_type} {old} which will be removed in version {version}, use {new} instead"
    if reason is not None:
        msg = f"{msg}\nreason={reason!r}"

    with catch_warnings():
        simplefilter("always")
        _warn(msg, DeprecationWarning, stacklevel=stack_level)


def discontinued(_type, old, version, reason=None, stack_level=3):
    """convenience func to warn about discontinued attributes

    Parameters
    ----------
    _type
        should be one of class, method, function, argument
    old
        the attributes name
    version
        the version by which support for the old name will be
        discontinued
    reason
        why, and what choices users have
    stack_level
        as per warnings.warn
    """
    msg = f"{_type} {old} is discontinued and will be removed in version {version}"
    if reason is not None:
        msg = f"{msg}\nreason={reason!r}"

    with catch_warnings():
        simplefilter("always")
        _warn(msg, DeprecationWarning, stacklevel=stack_level)


def deprecated_args(
    version: str,
    reason: str,
    arguments: List[Union[Tuple[str, str], str]] = None,
) -> Callable[..., Any]:
    """
    A decorator that marks specific arguments of a function as deprecated.

    The decorator accepts a list of 2-tuples specifying the mapping of old argument names to new argument names.
    When the decorated function is called with any of the old argument names, they will be replaced with their
    corresponding new names in the kwargs dictionary.

    Parameters
    ----------
    version : str
        The version when the old arguments will be removed in calver format, e.g. 'YYYY.MM'
    reason : str
        Reason for deprecation or guidance on what to do
    mapping : List[Union[Tuple[str, str],str]]
        A list of discontinued old argument names, or deprecated old and replacement new argument names.

    Returns
    -------
    Callable[..., Any]
        The decorated function.

    Warnings
    --------
    DeprecationWarning
        A warning will be raised when the decorated function is called for each deprecated argument
        used in the calling function.

    Examples
    --------
    Here's an example of how to use the `deprecated_args` decorator to mark the argument `old_name` as deprecated
    and replace it with the new name `new_name`.

    >>> @deprecated_args('2.0', 'Use new_name instead',[('old_arg', 'new_arg'),'discontinued_arg'],)
    >>> def my_function(new_name):
    >>>     # do something here

    When `my_function` is called with the argument `old_arg`, a warning will be raised indicating that the argument
    is deprecated and should be replaced with `new_arg`, and `discontinued_arg` is to be deprecated without replacement.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
            for argument in arguments:
                if isinstance(argument, tuple):
                    old, new = argument
                    if old in kwargs:
                        kwargs[new] = kwargs.pop(old)
                        deprecated("argument", old, new, version, reason)
                else:
                    if argument in kwargs:
                        discontinued("argument", argument, version, reason)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def deprecated_callable(
    version: str,
    reason: str,
    new: Optional[str] = None,
    is_discontinued: bool = False,
    stack_level=2,
) -> Callable:
    """
    A decorator that marks callables (function or method) as deprecated or discontinued..
    Parameters
    ----------
    version : str
        The version when it will be removed in calver format, e.g. 'YYYY.MM'
    reason : str
        Reason for deprecation or guidance on what to do
    new : str
        If the callable is being replaced, this is the replacement, e.g. 'ClassName.new_method()'
    is_discontinued : bool
        If True the callable is being discontinued.
    stack_level
        as per warnings.warn

    Returns
    -------
    Callable
        The decorated callable.

    Warnings
    --------
    DeprecationWarning
        A warning will be raised when the decorated function is called.

    Examples
    --------
    Here's an example of how to use the `deprecated_callable` decorator to mark the function `my_function` as deprecated
    in favour of a new function.

    >>> @deprecated_callable(version='2023.6', reason='function rename', new='a_function')
    >>> def my_function(arg): pass

    """

    def decorator(func: Callable) -> Callable:
        sig = set(inspect.signature(func).parameters)
        _type = "method" if sig & {"self", "cls", "klass"} else "function"
        old = func.__name__
        params = dict(
            _type=_type,
            old=old,
            version=version,
            reason=reason,
            stack_level=stack_level,
        )
        if is_discontinued:
            depr_func = discontinued
        else:
            params["new"] = new
            depr_func = deprecated

        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Callable:
            depr_func(**params)
            return func(*args, **kwargs)

        return wrapper

    return decorator

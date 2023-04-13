import functools
import inspect

from typing import Any, Callable, List, Optional, Tuple
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
        msg = f"{msg}\n{reason}"

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


def deprecated_callable(
    version: str,
    reason: str,
    new: Optional[str] = None,
    is_discontinued: bool = False,
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
        _type = "method" if inspect.ismethod(func) else "function"
        old = func.__name__
        params = dict(
            _type=_type, old=old, version=version, reason=reason, stack_level=2
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

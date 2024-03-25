import functools
import inspect

from typing import Any, Callable, List, Optional, Tuple
from warnings import catch_warnings, simplefilter
from warnings import warn as _warn


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


_discontinued = discontinued  # renamed to avoid name clash with discontinued argument in deprecated args decorator


def deprecated_args(
    version: str,
    reason: str,
    old_new: List[Tuple[str, str]] = None,
    discontinued: List[str] = None,
    stack_level=2,
) -> Callable[..., Any]:
    """
    A decorator that marks specific arguments of a function as deprecated.

    The decorator accepts a list of 2-tuples specifying the mapping of old
    argument names to new argument names. When the decorated function is
    called with any of the old argument names, they will be replaced with their
    corresponding new names in the kwargs dictionary.

    Parameters
    ----------
    version : str
        The version when the old arguments will be removed in calver
        format, e.g. 'YYYY.MM'
    reason : str
        Reason for deprecation or guidance on what to do
    old-new : List[Tuple[str, str]]
        A list of deprecated old and replacement new argument names.
    discontinued : List[str]
        Names of single or multiple arguments to be discontinued. This should
        only be applied to arguments that have no effect.

    Returns
    -------
    Callable[..., Any]
        The decorated function.

    Warnings
    --------
    DeprecationWarning
        A warning will be raised when the decorated function is called for
        each deprecated argument used in the calling function.

    Examples
    --------
    To use, change the signature of the function / method by removing the
    deprecated / discontinued arguments. Apply the decorator to the function,
    indicating the old and new the argument names.

    >>> @deprecated_args('2024.1',
    ...                  'Use new_name instead',
    ...                  old_new=[('old_arg', 'new_arg')],
    ...                  discontinued='discontinued_arg',
    ... )
    >>> def my_function(new_name):
    >>>     # do something here

    When `my_function` is called with the argument `old_arg`, a warning
    will be raised indicating that the argument is deprecated and should
    be replaced with `new_arg`, and `discontinued_arg` is to be
    discontinued.
    """

    discontinued = [discontinued] if isinstance(discontinued, str) else discontinued
    old_args = dict(old_new).keys() if old_new else set()

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
            if old_args & kwargs.keys():
                for old, new in old_new:
                    if old in kwargs:
                        kwargs[new] = kwargs.pop(old)
                        deprecated(
                            "argument",
                            old,
                            new,
                            version,
                            reason,
                            stack_level=stack_level,
                        )
            if discontinued:
                for dropped in discontinued:
                    if dropped in kwargs:
                        _discontinued(
                            "argument",
                            dropped,
                            version,
                            reason,
                            stack_level=stack_level,
                        )

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
        if is_discontinued and old == "__init__":
            # we're really deprecating a class, so get that name
            old = func.__qualname__.split(".")[-2]
            _type = "class"

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

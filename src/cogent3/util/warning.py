#!/usr/bin/env python
from warnings import catch_warnings, simplefilter
from warnings import warn as _warn


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jai Ram Rideout"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def deprecated(_type, old, new, version, reason=None, stack_level=2):
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
        _warn(msg, DeprecationWarning, stack_level)


def discontinued(_type, name, version, reason=None, stack_level=2):
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
        _warn(msg, DeprecationWarning, stack_level)

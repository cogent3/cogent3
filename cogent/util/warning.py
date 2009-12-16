#!/usr/bin/env python
from warnings import warn as _warn

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def deprecated(_type, old, new, version, stack_level=2):
    """a convenience function for deprecating classes, functions, arguments.
    
    Arguments:
        - _type should be one of class, method, function, argument
        - old, new: the old and new names
        - version: the version by which support for the old name will be
          discontinued
        - stack_level: as per warnings.warn"""
    msg = "use %s %s instead of %s, support discontinued in version %s" % \
        (_type, new, old, version)
    _warn(msg, DeprecationWarning, stack_level)

def discontinued(_type, name, version, stack_level=2):
    """convenience func to warn about discontinued attributes
    Arguments:
        - _type should be one of class, method, function, argument
        - name: the attributes name
        - version: the version by which support for the old name will be
          discontinued
        - stack_level: as per warnings.warn"""
    msg = "%s %s is discontinued, support will be stopped in version %s" %\
        (_type, name, version)
    _warn(msg, DeprecationWarning, stack_level)

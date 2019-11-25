#!/usr/bin/env python
"""Compiled modules may be out of date or missing"""

import os
import sys


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class ExpectedImportError(ImportError):
    pass


def fail(msg):
    print(msg, file=sys.stderr)
    raise ExpectedImportError


def importVersionedModule(name, exec_globals, min_version, alt_desc):
    if "COGENT3_PURE_PYTHON" in os.environ:
        fail('Not using compiled module "%s".  Will use %s.' % (name, alt_desc))
    try:
        m = __import__(name, exec_globals)
    except ImportError:
        fail('Compiled module "%s" not found.  Will use %s.' % (name, alt_desc))
    version = getattr(m, "version_info", (0, 0))
    desc = ".".join(str(n) for n in version)
    min_desc = ".".join(str(n) for n in min_version)
    max_desc = str(min_version[0]) + ".x"
    if version < min_version:
        fail(
            'Compiled module "%s" is too old as %s < %s. '
            "Will use %s." % (name, desc, min_desc, alt_desc)
        )
    if version[0] > min_version[0]:
        fail(
            'Compiled module "%s" is too new as %s > %s. '
            "Will use %s." % (name, desc, max_desc, alt_desc)
        )
    return m

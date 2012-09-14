#!/usr/bin/env python
"""Compiled modules may be out of date or missing"""

import os, sys

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class ExpectedImportError(ImportError):
    pass
        
def fail(msg):
    print >>sys.stderr, msg
    raise ExpectedImportError

def importVersionedModule(name, globals, min_version, alt_desc):
    if os.environ.has_key('COGENT_PURE_PYTHON'):
        fail('Not using compiled module "%s".  Will use %s.' % 
                (name, alt_desc))
    try:
        m = __import__(name, globals)
    except ImportError:
        fail('Compiled module "%s" not found.  Will use %s.' % 
                (name, alt_desc))
    version = getattr(m, 'version_info', (0, 0))
    desc = '.'.join(str(n) for n in version)
    min_desc = '.'.join(str(n) for n in min_version)
    max_desc = str(min_version[0])+'.x'
    if version < min_version:
        fail('Compiled module "%s" is too old as %s < %s. '
                'Will use %s.' % (name, desc, min_desc, alt_desc))
    if version[0] > min_version[0]:
        fail('Compiled module "%s" is too new as %s > %s. '
                'Will use %s.' % (name, desc, max_desc, alt_desc))
    return m

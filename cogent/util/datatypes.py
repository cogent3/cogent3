#!/usr/bin/env python

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class ImmutableDictionary(dict):
    def _immutable(self, *args, **kw):
        raise TypeError("%ss are immutable" % type(self).__name__)
    __setitem__ = __delitem__ = _immutable
    update = clear = pop = popitem = setdefault = _immutable

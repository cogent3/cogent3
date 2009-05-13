#!/usr/bin/env python
"""
This takes doctest files and turns them into standalone scripts.
"""
import doctest, sys, os

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

for filename in sys.argv[1:]:
    print filename, 
    (name, suffix) = os.path.splitext(filename)
    if suffix != '.rst':
        print 'not a .rst file'
        continue
    f = open(filename,'r')
    s = ''.join(f.readlines())
    f.close()
    
    s = doctest.script_from_examples(s)
    f = open(name+'.py','w')
    f.write(s)
    f.close()
    print '->', name+'.py'

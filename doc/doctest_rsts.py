#!/usr/bin/env python
import doctest, os, sys

"""
This will doctest all files ending with .rst in this directory.
"""

fl = sys.argv[1:]
if not fl:
    # find all files that end with rest
    fl = [fname for fname in os.listdir(os.getcwd()) if fname.endswith('.rst')]

    
for test in fl:
    doctest.testfile(test, 
    optionflags = doctest.ELLIPSIS, verbose=True)

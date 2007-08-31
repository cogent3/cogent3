#!/usr/bin/env python
import doctest, os, sys

"""
This will doctest all files ending with .rest in this directory.
"""

fl = sys.argv[1:]
if not fl:
    # find all files that end with rest
    fl = [fname for fname in os.listdir(os.getcwd()) if fname.endswith('.rest')]

    
for test in fl:
    doctest.testfile(test, 
    optionflags = doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS, verbose=True)

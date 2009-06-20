#!/usr/bin/env python

"""An example of how to distribute jobs across multiple cpu's. Note that this example works even if on a single cpu since the parallel assigns a `fake` communicator in that instance.
"""

from cogent.util import parallel

__author__ = "Gavin Huttleu"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

JOBS = range(12) # nonsense jobs, just a list of numbers to be printed.

# we divide up the CPUs into (at most) 12 groups of (at least) 1 CPU.
(comm, leftover) = parallel.getSplitCommunicators(len(JOBS))
# and set the cpu's available to lower levels
parallel.push(leftover)
try:
    for job in JOBS:
        if not job % comm.size == comm.rank:
            continue
        print "My ID=%d, my message=%s" % (comm.rank, JOBS[job])
finally:
    # always restore the original parallel context
    parallel.pop(leftover)
    

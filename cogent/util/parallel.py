#!/usr/bin/env python
from __future__ import with_statement
import os, sys
from contextlib import contextmanager
import warnings

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley",
                "Matthew Wakefield", "Edward Lang"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin Huttley"
__status__ = "Production"

# A flag to control if excess CPUs are worth a warning.
inefficiency_forgiven = False

class _FakeCommunicator(object):
    """Looks like a 1-cpu MPI communicator, but isn't"""
    def Get_rank(self):
        return 0
    def Get_size(self):
        return 1
    def Split(self, colour, key=0):
        return (self, self)
    def allreduce(self, value, op=None):
        return value
    def allgather(self, value):
        return [value]
    def bcast(self, obj, source):
        return obj
    def Bcast(self, array, source):
        pass
    def Barrier(self):
        pass

class _FakeMPI(object):
    # required MPI module constants
    SUM = MAX = DOUBLE = 'fake'   

if os.environ.get('DONT_USE_MPI', 0):
    print >>sys.stderr, 'Not using MPI'
    MPI = None
else:
    try:
        from mpi4py import MPI
    except ImportError:
        warnings.warn('Not using MPI as mpi4py not found', stacklevel=2)
        MPI = None
    else:
        size = MPI.COMM_WORLD.Get_size()
        if size == 1:
            MPI = None

if MPI is None:
    def get_processor_name():
        return os.environ.get('HOSTNAME', 'one')
    _ParallelisationStack = [_FakeCommunicator()]
    MPI = _FakeMPI()
else:
    get_processor_name = MPI.Get_processor_name
    _ParallelisationStack = [MPI.COMM_WORLD]
        
def sync_random(r):
    if _ParallelisationStack[-1].Get_size() > 1:
        state = _ParallelisationStack[-1].bcast(r.getstate(), 0)
        r.setstate(state)

def getCommunicator():
    return _ParallelisationStack[-1]

def getSplitCommunicators(jobs):
    comm = getCommunicator()
    assert jobs > 0
    (size, rank) = (comm.Get_size(), comm.Get_rank())
    group_count = min(jobs, size)
    while size % group_count:
        group_count += 1
    if group_count == 1:
        (next, sub) = (_FakeCommunicator(), comm)
    elif group_count == size:
        (next, sub) = (comm, _FakeCommunicator())
    else:
        next = comm.Split(rank // group_count, rank)
        sub = comm.Split(rank % group_count, rank)
    return (next, sub)

@contextmanager
def mpi_context(comm):
    _ParallelisationStack.append(comm)
    try:
        yield
    finally:
        popped = _ParallelisationStack.pop()    
        assert popped is comm

@contextmanager
def mpi_split(jobs):
    (next, sub) = getSplitCommunicators(jobs)
    with mpi_context(sub):
        yield next

def map(f,s):
    result = []
    with mpi_split(len(s)) as comm:
        (size, rank) = (comm.Get_size(), comm.Get_rank())
        for start in range(0, len(s), size):
            chunk = s[start:start+size]
            if rank < len(chunk):
                local_result = f(chunk[rank])
            else:
                local_result = None
            split_results = comm.allgather(local_result)[:len(chunk)]
            result.extend(split_results)
    return result

output_cpu = getCommunicator().Get_rank() == 0

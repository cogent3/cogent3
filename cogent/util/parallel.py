#!/usr/bin/env python
from __future__ import with_statement
import os, sys
from contextlib import contextmanager
import warnings
import threading
import multiprocessing
import multiprocessing.pool

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley",
                "Matthew Wakefield", "Edward Lang"]
__license__ = "GPL"
__version__ = "1.5.3"
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

FAKE_MPI_COMM = _FakeCommunicator()

class _FakeMPI(object):
    # required MPI module constants
    SUM = MAX = DOUBLE = 'fake'
    COMM_WORLD = FAKE_MPI_COMM

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
    USING_MPI = False
    MPI = _FakeMPI()
else:
    USING_MPI = True


class ParallelContext(object):
    """A parallel context encapsulates the number of CPUs available and the
    mechanism by which they communicate.  All contexts offer the lowest-common-
    denominator of parallel mechanisms - parallel imap."""
    pass


class NonParallelContext(ParallelContext):
    """This dummy parallel context is used when there is only one CPU available
    """
    size = 1
    
    def getCommunicator(self):
        return FAKE_MPI_COMM

    def split(self, jobs):
        return (self, self)
    
    def imap(self, f, s, chunksize=None):
        for element in s:
            yield f(element)

        
NONE = NonParallelContext()

class UnFlattened(list):
    pass
    
class MPIParallelContext(ParallelContext):
    """This parallel context divides the available CPUs into groups of equal 
    size.  Inner levels of potential parallelism can then further subdivide
    those groups.  It helps to have a CPU count which is divisible by the
    task count."""
    
    def __init__(self, comm=None):
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm
        self.size = comm.Get_size()
    
    def getCommunicator(self):
        return self.comm
    
    def split(self, jobs):
        assert jobs > 0
        size = self.size
        group_count = min(jobs, size)
        while size % group_count:
            group_count += 1
        if group_count == 1:
            (next, sub) = (NONE, self)
        elif group_count == size:
            (next, sub) = (self, NONE)
        else:
            rank = self.comm.Get_rank()
            klass = type(self)
            next = klass(self.comm.Split(rank // group_count, rank))
            sub = klass(self.comm.Split(rank % group_count, rank))
        return (next, sub)
        
    def imap(self, f, s, chunksize=1):
        comm = self.comm
        (size, rank) = (comm.Get_size(), comm.Get_rank())
        ordinals = range(0, len(s), size*chunksize)
        # ensure same number of allgather calls in every process
        for start in ordinals:
            start += rank
            local_results = UnFlattened([f(x) for x in s[start:start+chunksize]])
            for results in comm.allgather(local_results):
                # mpi4py allgather has a nasty inconsistancy about flattening
                # lists of simple values
                if isinstance(results, UnFlattened):
                    for result in results:
                        yield result
                else:
                    yield results


# Helping MultiprocessingParallelContext map unpicklable functions
_FUNCTIONS = {}
class PicklableAndCallable(object):
    def __init__(self, key):
        self.key = key
        self.func = None
    def __call__(self, *args, **kw):
        if self.func is None:
            try:
                self.func = _FUNCTIONS[self.key]
            except KeyError:
                raise RuntimeError
        return self.func(*args, **kw)
    
class MultiprocessingParallelContext(ParallelContext):
    """At the outermost opportunity, this parallel context delegates all
    work to a multiprocessing.Pool.  
    Subprocesses may also make pools if the outer pool is more than half idle.
    
    Ideally the pool would be reused for later tasks, but cogent code mostly 
    uses map() with functions defined in local scopes, which are unpicklable,
    so that is hacked around and pools are only ever used for one map() call"""
    
    def __init__(self, size=None):
        if size is None:
            size = multiprocessing.cpu_count()
        self.size = size
    
    def getCommunicator(self):
        return FAKE_MPI_COMM

    def _subContext(self, size):
        if size == 1:
            return NONE
        elif size == self.size:
            return self
        else:
            return type(self)(size)
        
    def split(self, jobs):
        assert jobs > 0
        group_count = min(self.size, jobs)
        remaining = self.size // group_count
        next = self._subContext(group_count)
        sub = self._subContext(remaining)
        return (next, sub)

    def _initWorkerProcess(self):
        from cogent.util import progress_display
        progress_display.CURRENT.context = progress_display.NULL_CONTEXT

    def imap(self, f, s, chunksize=1):
        key = id(f)
        _FUNCTIONS[key] = f
        f = PicklableAndCallable(id(f))
        pool = multiprocessing.Pool(self.size, self._initWorkerProcess)
        for result in pool.imap(f, s, chunksize=chunksize):
            yield result
        del _FUNCTIONS[key]
        pool.close()


class ContextStack(threading.local):
    """This singleton object holds the current and enclosing parallel contexts."""
    
    def __init__(self):
        # Because this is a thread.local, any secondary threads will see this 
        # default and so not attempt to use MPI/multiprocessing:
        self.stack = []
        self.top = NONE
    
    def setInitial(self, context):
        """The real initialiser.  Should be called once from the main thread."""
        assert self.stack == [] and self.top is NONE
        self.top = context
        
    @contextmanager
    def pushed(self, context):
        """Temporarily enter a pre-existing parallel context"""
        self.stack.append(self.top)
        try:
            self.top = context
            yield
        finally:
            self.top = self.stack.pop()

    @contextmanager
    def split(self, jobs=None):
        """Divide the available CPUs up into groups to handle 'jobs' independent
        tasks.  If jobs << CPUs so that there are multiple CPUS per job, leave
        that reduced number of CPUs available to any nested parallelism 
        opportunities within each job"""
        if jobs is None:
            jobs = self.top.size
        (next, sub) = self.top.split(jobs)
        with self.pushed(sub):
            yield next

    def imap(self, f, s, chunksize=1):
        """Like itertools.imap(f,s) only parallel."""
        chunks = (len(s)-1) // chunksize + 1
        with self.split(chunks) as next:
            for element in next.imap(f, s, chunksize=chunksize):
                yield element

    def map(self, f, s):
        return list(self.imap(f, s))
    
    def getCommunicator(self):
        """For code needing an MPI communicator interface.  If not
        using MPI this will be a dummy communicator of 1 CPU."""
        return self.top.getCommunicator()
        
    def getContext(self):
        return self.top

CONTEXT = ContextStack()


if os.environ.get('COGENT_CPUS', False):
    try:
        cpus = int(os.environ['COGENT_CPUS'])
    except ValueError:
        cpus = None
    else:
        assert cpus > 0
else:
    cpus = 1

if USING_MPI:
    CONTEXT.setInitial(MPIParallelContext())
elif cpus > 1:
    CONTEXT.setInitial(MultiprocessingParallelContext(cpus))


getContext = CONTEXT.getContext
getCommunicator = CONTEXT.getCommunicator
parallel_context = CONTEXT.pushed
split = CONTEXT.split
imap = CONTEXT.imap
map = CONTEXT.map

def use_multiprocessing(cpus=None):
    CONTEXT.setInitial(MultiprocessingParallelContext(cpus))

def sync_random(r):
    # Only matters with MPI
    comm = getCommunicator()
    state = comm.bcast(r.getstate(), 0)
    r.setstate(state)


    

#!/usr/bin/env python

import os
import sys
from contextlib import contextmanager
import warnings
import threading
import multiprocessing
import concurrent.futures as concurrentfutures

__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Sheng Han Moses Koh"]
__license__ = "GPL"
__version__ = "3.0a2"

class _FakeCommunicator(object):
    """Looks like a 1-cpu MPI communicator, but isn't"""

    def Get_rank(self):
        return 0

    def Get_size(self):
        if USING_MPI:
            return MPI.INFO_ENV.get('maxprocs', 1)
        else:
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
    INFO_ENV = {'maxprocs':multiprocessing.cpu_count()}

if os.environ.get('DONT_USE_MPI', 0):
    warnings.warn('Not using MPI', stacklevel=2)
    MPI = None
else:
    try:
        from mpi4py import MPI
        from mpi4py import futures as MPIfutures
    except ImportError:
        warnings.warn('Not using MPI as mpi4py not found', stacklevel=2)
        MPI = None
    else:
        size = MPI.INFO_ENV.get("maxprocs", 1)
        if size == 1:
            MPI = None

MPI = _FakeMPI()
if MPI is None:
    USING_MPI = False
else:
    USING_MPI = True

# Helping ProcessPoolExecutor map unpicklable functions
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

def imap(f, s, max_workers=None, use_mpi=False):
    if use_mpi:
        if not max_workers:
            raise ValueError
        assert max_workers < MPI.INFO_ENV.get("maxprocs")
        with MPIfutures.MPIPoolExecutor(max_workers) as executor:
            for result in executor.map(f, s):
                yield result
    else:
        # If max_workers is not defined, get number of all processes available
        # minus 1 to leave for master process
        if not max_workers:
            max_workers = multiprocessing.cpu_count() - 1
        assert max_workers < multiprocessing.cpu_count()
        #max_workers = max_workers or (multiprocessing.cpu_count() - 1)
        key = id(f)
        _FUNCTIONS[key] = f
        f = PicklableAndCallable(id(f))
        with concurrentfutures.ProcessPoolExecutor(max_workers) as executor:
            for result in executor.map(f, s):
                yield result

def map(f, s, max_workers=None, use_mpi=False):
    return list(imap(f, s, max_workers, use_mpi))

def get_communicator():
    return FAKE_MPI_COMM

#!/usr/bin/env python

import os
import sys
import time
import math
import threading
import multiprocessing
import warnings
import concurrent.futures as concurrentfutures

__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Sheng Han Moses Koh"]
__license__ = "GPL"
__version__ = "3.0a2"

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

COMM = MPI.COMM_WORLD


def generateRandomSeed(use_mpi):
    global ran_seed
    if use_mpi:
        rank = COMM.Get_rank()
    else:
        processName = multiprocessing.current_process().name
        rank = int(processName[-1])
    ran_seed = int(time.time()) + rank
if MPI is not None:
    generateRandomSeed(True)
    USING_MPI = True
else:
    USING_MPI = False

# Helping ProcessPoolExecutor map unpicklable functions
_FUNCTIONS = {}


class PicklableAndCallable():

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
    # If max_workers is not defined, get number of all processes available
    # minus 1 to leave for master process
    if use_mpi:
        if not USING_MPI:
            raise RuntimeError
        if not max_workers:
            max_workers = multiprocessing.cpu_count() - 1
        with MPIfutures.MPIPoolExecutor(max_workers) as executor:
            for result in executor.map(f, s):
                yield result
    else:
        if not max_workers:
            max_workers = multiprocessing.cpu_count() - 1
        assert max_workers < multiprocessing.cpu_count()
        key = id(f)
        _FUNCTIONS[key] = f
        f = PicklableAndCallable(id(f))
        with concurrentfutures.ProcessPoolExecutor(max_workers,
                                    initializer=generateRandomSeed,
                                    initargs=([False])) as executor:
            for result in executor.map(f, s):
                yield result


def map(f, s, max_workers=None, use_mpi=False):
    return list(imap(f, s, max_workers, use_mpi))

Distribution of tasks across multiple CPUs using ``multiprocess``
=================================================================

.. note:: Using multiprocess requires no extra commands on invoking the script!

In the ``multiprocess`` case,  ``foo()`` is called in a *temporary* subprocesses. Communication from a subprocess back into the top process is via the value that ``foo()`` returns.

I illustrate the use of ``parallel.map`` here with an example that collects both the process id and the integer.

.. doctest::
    
    >>> from cogent.util import parallel
    >>> import os, time
    >>> parallel.use_multiprocessing(2)
    >>> def foo(val):
    ...     # do some real intensive work in here, which I've simulated by
    ...     # using sleep
    ...     time.sleep(0.01)
    ...     return os.getpid(), val
    >>> data = range(20)
    >>> result = parallel.map(foo, data)

For the purpose of testing the code is executing correctly, I'll check that there are 2 pid's returned.

.. doctest::
    
    >>> pids, nums = zip(*result)
    >>> len(set(pids)) == 2
    True

The ``result`` list is a series of tuples of process id's and the integer, the latter not necessarily being in order.

.. we don't test the following output since the pid will vary between runs

.. doctest::
    
    >>> print result # doctest: +SKIP
    [(7332, 0), (7333, 1), (7332, 2), (7333, 3), (7332, 4), (7333, 5), (7332, 6), (7333, 7), (7332, 8), (7333, 9), (7332, 10), (7333, 11), (7332, 12), (7333, 13), (7332, 14), (7333, 15), (7332, 16), (7333, 17), (7332, 18), (7333, 19)]


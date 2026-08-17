.. _parallel:

Parallel computations
=====================

``cogent3`` supports parallel computation explicitly for the case where the same calculations need to be performed on many different data sets. As an example, consider the case of aligning all the one-to-one orthologs of protein coding genes sampled from 100 vertebrate species where the data for each gene is stored in a separate text file. These files are used as input for an alignment algorithm that will produce a corresponding output file. In other words, applying the alignment algorithm to ``"homologs1.fasta"`` produces ``"aligned-homologs1.fasta"``.

We could perform the alignments in serial, one after the other, on one CPU core of a single computer. But what if we have 18,000 such files? If we had 18,000 CPUs then we could assign one alignment task to each file and be done in the same time as aligning a single file! This case is an example of "data parallelism" or "data level parallelism".

There are multiple algorithmic approaches to solving parallel computation problems. The approach ``cogent3`` adopts is that of a master process and helper (or worker) processes. The master process splits the work up amongst the available CPU cores. Using our alignment example, the master process assigns sets of files to each worker CPU core. Each worker then performs the alignment step on its designated files and returns each alignment to the master process.

.. warning:: 
    
    It is not always faster to split tasks between processes. You should see a performance gain if the calculation time per task of the worker is significantly greater than the time it will take the master process to deal with the result -- in our example, the time it takes to write the alignment to file.

    While the alignment problem indicated above stipulated writing all results to separate files, this is not always a good idea. It can prove very inefficient if the individual alignment files are small. In such a case, storing the result in a single file (e.g. as a ``sqlitedb`` database) is better.

Parallel computation on a single computer
-----------------------------------------

This is the simplest case to implement, requires no additional software installs and will work with standalone scripts or within Jupyter notebooks. For this use case, ``cogent3.util.parallel`` uses the Python standard library ``concurrent.futures`` module.

Using ``app.apply_to()``
^^^^^^^^^^^^^^^^^^^^^^^^

If you are using a composed ``cogent3`` :ref:`app <apps>` **with** a writer, then the simplest approach is to use the ``apply_to()`` method. The conditions of parallel execution are controlled using the keyword arguments ``parallel`` and ``par_kw``. The former indicates parallel execution is to be undertaken. The latter is how additional arguments are provided to ``parallel.as_completed()``. For instance, using 4 workers would be specified as:

.. code-block:: python
    
    results = app.apply_to(data, parallel=True, par_kw=dict(max_workers=4))

.. note:: If you are using mpi, set ``par_kw=dict(max_workers=4, use_mpi=True)``.

Using ``app.as_completed()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are using a composed ``cogent3`` :ref:`app <apps>` **without** a writer, then use the ``as_completed()`` method. The arguments are the same as for ``apply_to()`` but as this method returns a generator, you use a builtin type to execute the call.

.. code-block:: python
    
    results = list(app.as_completed(data, parallel=True, par_kw=dict(max_workers=4)))

.. note:: If you are using mpi, set ``par_kw=dict(max_workers=4, use_mpi=True)``.

Directly using ``cogent3.util.parallel.as_completed()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This function enables distribution of calculations across CPUs.

.. warning:: This function delivers results in the order completed! This is different to ``cogent3.util.parallel.map()`` which delivers results in the order provided.

The demo script shown below calculates a small number of prime numbers by splitting chunks of numbers across the provided cores. The key line is

.. code-block:: python

    result = parallel.as_completed(is_prime, PRIMES, max_workers=4)

The first argument, ``is_prime``, is the function to be called with values from the data, ``PRIMES``. The ``max_workers`` argument indicates how many worker processes to use. The elements of ``PRIMES`` will be broken into ``max_workers`` number of equal sized chunks. Each such chunk is applied to ``is_prime`` on a separate CPU. In this case, the returned results will be a series of ``bool`` values.

.. note:: If you don't specify ``max_workers``, all available CPUs will be used.

.. literalinclude:: demo-multiprocess-parallel.py

Parallel computation on multiple computers
------------------------------------------

On systems consisting of multiple separate computers, we use the mpi4py_ bindings to the message passing interface (MPI) standard. Specifically, ``cogent3.util.parallel.map(..., use_mpi=True, ...)`` uses the `mpi4py futures`_ module of mpi4py_. This module is modelled after that of ``concurrent.futures`` but using it has some important differences.

First, you must install additional software. You will need to install a tool implementing the MPI system (e.g. `openmpi <https://www.open-mpi.org/>`_) and the MPI python bindings library ``mpi4py``. To install ``openmpi``, you can use conda, homebrew or your preferred package manager. You can just pip install ``mpi4py``.

Second, as described in the documentation on `mpi4py futures`_, you need to write your code slightly differently. We provide an example that runs on a supercomputer. To execute a program on this facility, we submit a "job" to a "queuing system" (e.g. `PBS <https://en.wikipedia.org/wiki/Portable_Batch_System>`_) which controls the scheduling of our job with the computing resources we requested (how many CPUs, how much RAM, etc..). There are many such job control systems and the specifics of how to select the resources your job needs can vary between them. In general, however, our experience is the user writes two scripts.

1. a script performing the computations you actually care about
2. a bash script for the queuing system setting out the job parameters and invoking (1)

The example code presented below is based on the ``mpi4py`` demo script for computing prime numbers. In addition to validating the prime numbers, it also prints out the "MPI rank" of the processor [1]_. The script relies on the environment variable, ``PBS_NCPUS`` [2]_, to establish the number of CPUs that are available. It prints to stdout, the rank of each processor [3]_. 

To execute this script as part of a PBS job script you need to use the following command::

$ mpiexec -n $PBS_NCPUS python3 -m mpi4py.futures demo-mpi-parallel.py

.. note::

    To execute it directly with 4 CPUs do::
    
    $ PBS_NCPUS=4 mpiexec -n 4 python3 -m mpi4py.futures demo-mpi-parallel.py

The ``-n`` argument tells ``mpiexec`` to use this number of CPUs.

In the ``demo-mpi-parallel.py`` script, the key line is

.. code-block:: python

    result = parallel.map(is_prime, PRIMES, use_mpi=True, max_workers=PBS_CPUS)

The ``use_mpi`` argument invokes the correct back end, otherwise the interface is the same as described above.

.. note:: You can use ``mpi`` for parallel execution on a single computer. This can be useful for checking your code prior to migrating to a larger system.

.. literalinclude:: demo-mpi-parallel.py

.. _mpi4py futures: https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html
.. _mpi4py: https://mpi4py.readthedocs.io/

.. [1] On MPI, the main process has rank 0, all others have rank > 0.
.. [2] This environment variable is created by the PBS system on executing the job script.
.. [3] You can check your execution of the script is correct by validating you get all the ranks up to one minus the number of CPUs you requested.

Easy distribution of tasks across multiple CPUs
===============================================

.. warning:: This example requires execution on 2 CPUs. It can be run using: $ mpirun -n 2 python path/to/parallel_tasks.rst

All I'm doing here is illustrating the use of ``parallel.map`` with the simplest example I could come up with. I create a list which will have an integer appended to it -- hardly useful, but hopefully demonstrates how a series of data can be mapped onto a function in parallel. In this case, the data is just a list of numbers.

.. doctest::
    
    >>> from cogent.util import parallel
    >>> passed_indices = []
    >>> series = range(20)
    >>> result = parallel.map(passed_indices.append, series)
    >>> assert passed_indices == range(0,20,2) or passed_indices == range(1,20,2)

The result is either the list of even numbers up to (but not including) 20 or the list of odd numbers.

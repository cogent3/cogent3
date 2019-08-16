import time

from unittest import TestCase, main

import numpy

from cogent3.util import parallel


__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Sheng Han Moses Koh"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def get_process_value(n):
    # Sleep to accommodate Windows process creation overhead
    time.sleep(1)
    return (parallel.get_rank(), n)


def get_ranint(n):
    numpy.random.seed(n)
    return numpy.random.randint(1, 10)


class ParallelTests(TestCase):
    def test_create_processes(self):
        """Procressor pool should create multiple distingue processes"""
        index = [2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = parallel.map(get_process_value, index, max_workers=None, use_mpi=False)
        resultProcesses = [v[0] for v in result]
        resultValues = [v[1] for v in result]
        self.assertEqual(sorted(list(resultValues)), index)
        self.assertNotEqual(len(set(resultProcesses)), 1)

    def test_random_seeding(self):
        """Random seed should be set every function call"""
        # On Windows process ids are not guaranteed to be sequential(1,2,3,4...)
        # thus they cannot be used for reproducibility
        index1 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
        index2 = [2, 2, 2, 2, 2, 2, 2, 2, 2]
        result1 = parallel.map(get_ranint, index1, max_workers=1, use_mpi=False)
        result2 = parallel.map(get_ranint, index2, max_workers=1, use_mpi=False)
        self.assertEqual(result1[0], result2[0])
        self.assertNotEqual(result1, result2)


if __name__ == "__main__":
    main()

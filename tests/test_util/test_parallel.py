import multiprocessing
import sys
import time

from unittest import TestCase, main, skipIf

import numpy

from cogent3.util import parallel


__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Sheng Han Moses Koh"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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


def check_is_master_process(n):
    return parallel.is_master_process()


class ParallelTests(TestCase):
    def test_create_processes(self):
        """Procressor pool should create multiple distingue processes"""
        max_worker_count = multiprocessing.cpu_count() - 1
        index = list(range(max_worker_count))
        result = parallel.map(get_process_value, index, max_workers=None, use_mpi=False)
        result_processes = [v[0] for v in result]
        result_values = [v[1] for v in result]
        self.assertEqual(sorted(list(result_values)), index)
        self.assertEqual(len(set(result_processes)), max_worker_count)

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

    @skipIf(sys.version_info[1] < 7, "method exclusive to Python 3.7 and above")
    def test_is_master_process(self):
        """
        is_master_process() should return False
        for all child processes
        """
        index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        master_processes = 0
        for result in parallel.imap(
            check_is_master_process, index, max_workers=None, use_mpi=False
        ):
            if result:
                master_processes += 1
        self.assertEqual(master_processes, 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python


import os
import sys

from unittest import TestCase, main

from cogent3.maths.optimisers import MaximumEvaluationsReached, maximise


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def quartic(x):
    # Has global maximum at -4 and local maximum at 2
    # http://www.wolframalpha.com/input/?i=x**2*%283*x**2%2B8*x-48%29
    # Scaled down 10-fold to avoid having to change init_temp
    return x ** 2 * (3 * x ** 2 + 8 * x - 48)


class NullFile(object):
    def write(self, x):
        pass

    def isatty(self):
        return False


def quiet(f, *args, **kw):
    # Checkpointer still has print statements
    orig = sys.stdout
    try:
        sys.stdout = NullFile()
        result = f(*args, **kw)
    finally:
        sys.stdout = orig
    return result


def MakeF():
    evals = [0]
    last = [0]

    def f(x):
        evals[0] += 1
        last[0] = x
        # Scaled down 10-fold to avoid having to change init_temp
        return -0.1 * quartic(x)

    return f, last, evals


class OptimiserTestCase(TestCase):
    def _test_optimisation(self, target=-4, xinit=1.0, bounds=None, **kw):
        bounds = bounds or ([-10, 10])
        f, last, evals = MakeF()

        x = quiet(maximise, f, [xinit], bounds, **kw)
        self.assertEqual(x, last[0])  # important for Calculator
        error = abs(x[0] - target)
        self.assertTrue(error < 0.0001, (kw, x, target, x))

    def test_global(self):
        # Should find global minimum
        self._test_optimisation(local=False, seed=1)

    def test_bounded(self):
        # Global minimum out of bounds, so find secondary one
        # numpy.seterr('raise')
        self._test_optimisation(bounds=([0.0], [10.0]), target=2, seed=1)

    def test_local(self):
        # Global minimum not the nearest one
        self._test_optimisation(local=True, target=2)

    def test_limited(self):
        self.assertRaises(
            MaximumEvaluationsReached, self._test_optimisation, max_evaluations=5
        )

    # def test_limited_warning(self):
    #     """optimiser warning if max_evaluations exceeded"""
    #     self._test_optimisation(max_evaluations=5, limit_action='warn')

    def test_get_max_eval_count(self):
        """return the evaluation count from optimisation"""
        f, last, evals = MakeF()
        x, e = quiet(
            maximise, f, xinit=[1.0], bounds=([-10, 10]), return_eval_count=True
        )
        # picking arbitrary numerical value
        self.assertGreaterEqual(e, 10)

    def test_checkpointing(self):
        filename = "checkpoint.tmp.pickle"
        if os.path.exists(filename):
            os.remove(filename)
        self._test_optimisation(filename=filename, seed=1, init_temp=10)
        self._test_optimisation(filename=filename, seed=1, init_temp=10)
        self.assertRaises(
            Exception,
            self._test_optimisation,
            filename=filename,
            seed=1,
            init_temp=3.21,
        )
        if os.path.exists(filename):
            os.remove(filename)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Simulated annealing optimiser. Derives from basic optimiser class.

The simulated annealing optimiser is a translation into Python of the fortran
program simman.f authored by Bill Goffe (bgoffe@whale.st.usm.edu). The original
citation is "Global Optimization of Statistical Functions with Simulated
Annealing," Goffe, Ferrier and Rogers, Journal of Econometrics, vol. 60, no. 1/2,
Jan./Feb. 1994, pp. 65-100.
"""

import random

from collections import deque

import numpy

from cogent3.util import checkpointing


__author__ = "Andrew Butterfield and Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Andrew Butterfield", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class AnnealingSchedule(object):
    """Responsible for the shape of the simulated annealing temperature profile"""

    def __init__(self, temp_reduction, initial_temp, temp_iterations, step_cycles):
        if initial_temp < 0.0:
            raise ValueError("Initial temperature not +ve")
        self.T = self.initial_temp = initial_temp
        self.temp_reduction = temp_reduction
        self.temp_iterations = temp_iterations
        self.step_cycles = step_cycles
        self.dwell = temp_iterations * step_cycles

    def checkSameConditions(self, other):
        for attr in [
            "temp_reduction",
            "initial_temp",
            "temp_iterations",
            "step_cycles",
        ]:
            if getattr(self, attr) != getattr(other, attr):
                raise ValueError(f"Checkpoint file ignored - {attr} different")

    def roundsToReach(self, T):
        from math import log

        return int(-log(self.initial_temp / T) / log(self.temp_reduction)) + 1

    def cool(self):
        self.T = self.temp_reduction * self.T

    def willAccept(self, newF, oldF, random_series):
        deltaF = newF - oldF
        return deltaF >= 0 or random_series.uniform(0.0, 1.0) < numpy.exp(
            deltaF / self.T
        )


class AnnealingHistory(object):
    """Keeps the last few results, for convergence testing"""

    def __init__(self, sample=4):
        self.sample_size = sample
        # self.values = deque([None]*sample, sample) Py2.6
        self.values = deque([None] * sample)

    def note(self, F):
        self.values.append(F)
        # Next 2 lines not required once above Py2.6 line is uncommented
        if len(self.values) > self.sample_size:
            self.values.popleft()

    def minRemainingRounds(self, tolerance):
        last = self.values[-1]
        return max(
            [0]
            + [
                i + 1
                for (i, v) in enumerate(self.values)
                if v is None or abs(v - last) > tolerance
            ]
        )


class AnnealingState(object):
    def __init__(self, X, function, random_series):
        self.random_series = random_series
        self.NFCNEV = 1
        self.VM = numpy.ones(len(X), float)
        self.setX(X, function(X))
        (self.XOPT, self.FOPT) = (X, self.F)
        self.NACP = [0] * len(X)
        self.NTRY = 0

    def setX(self, X, F):
        self.X = numpy.array(X, float)
        self.F = F

    def step(self, function, accept_test):
        # One attempted move in each dimension
        X = self.X
        self.NTRY += 1
        for H in range(len(X)):
            self.NFCNEV += 1

            current_value = X[H]
            X[H] += self.VM[H] * self.random_series.uniform(-1.0, 1.0)
            F = function(X)

            if accept_test(F, self.F, self.random_series):
                self.NACP[H] += 1
                self.F = F
                if F > self.FOPT:
                    (self.FOPT, self.XOPT) = (F, X.copy())
            else:
                X[H] = current_value

    def adjustStepSizes(self):
        # Adjust velocity in each dimension to keep acceptance ratios near 50%
        if self.NTRY == 0:
            return
        for I in range(len(self.X)):
            RATIO = (self.NACP[I] * 1.0) / self.NTRY
            if RATIO > 0.6:
                self.VM[I] *= 1.0 + (2.0 * ((RATIO - 0.6) / 0.4))
            elif RATIO < 0.4:
                self.VM[I] /= 1.0 + (2.0 * ((0.4 - RATIO) / 0.4))
            self.NACP[I] = 0
        self.NTRY = 0


class AnnealingRun(object):
    def __init__(self, function, X, schedule, random_series):
        self.history = AnnealingHistory()
        self.schedule = schedule
        self.state = AnnealingState(X, function, random_series)
        self.test_count = 0

    def checkFunction(self, function, xopt, checkpointing_filename):
        if len(xopt) != len(self.state.XOPT):
            raise ValueError(
                "Number of parameters in checkpoint file '%s' (%s) "
                "don't match current function (%s)"
                % (checkpointing_filename, len(self.state.XOPT), len(xopt))
            )
        # if f(x) != g(x) then f isn't g.
        then = self.state.FOPT
        now = function(self.state.XOPT)
        if not numpy.allclose(now, then, 1e-8):
            raise ValueError(
                "Function to optimise doesn't match checkpoint file "
                "'%s': F=%s now, %s in file." % (checkpointing_filename, now, then)
            )

    def run(self, function, tolerance, checkpointer, show_remaining):
        state = self.state
        history = self.history
        schedule = self.schedule

        est_anneal_remaining = schedule.roundsToReach(tolerance / 10) + 3
        while True:
            min_history_remaining = history.minRemainingRounds(tolerance)
            if min_history_remaining == 0:
                break
            self.save(checkpointer)
            remaining = max(min_history_remaining, est_anneal_remaining)
            est_anneal_remaining += -1

            for i in range(self.schedule.dwell):
                show_remaining(
                    remaining + 1 - i / self.schedule.dwell,
                    state.FOPT,
                    schedule.T,
                    state.NFCNEV,
                )
                state.step(function, self.schedule.willAccept)
                self.test_count += 1
                if self.test_count % schedule.step_cycles == 0:
                    state.adjustStepSizes()

            history.note(state.F)
            state.setX(state.XOPT, state.FOPT)
            schedule.cool()

        self.save(checkpointer, final=True)

        return state

    def save(self, checkpointer, final=False):
        msg = "Number of function evaluations = %d; current F = %s" % (
            self.state.NFCNEV,
            self.state.FOPT,
        )
        checkpointer.record(self, msg, final)


class SimulatedAnnealing(object):
    """Simulated annealing optimiser for bounded functions"""

    def __init__(self, filename=None, interval=None, restore=True):
        """
        Set the checkpointing filename and time interval.

        Parameters
        ----------
        filename
            name of the file to which data will be written. If None, no
            checkpointing will be done.
        interval
            time expressed in seconds
        restore
            flag to restore from this filename or not. will be set to 0 after
            restoration

        """
        self.checkpointer = checkpointing.Checkpointer(filename, interval)
        self.restore = restore

    def maximise(
        self,
        function,
        xopt,
        show_remaining,
        random_series=None,
        seed=None,
        tolerance=None,
        temp_reduction=0.5,
        init_temp=5.0,
        temp_iterations=5,
        step_cycles=20,
    ):
        """Optimise function(xopt).

        Parameters
        ----------
        show_progress
            whether the function values are printed as
            the optimisation proceeds. Default is True.
        tolerance
            the error condition for termination, default is 1E
        temp_reduction
            the factor by which the annealing
            "temperature" is reduced, default is 0.5
        temp_iterations
            the number of iterations before a
            temperature reduction, default is 5
        step_cycles
            the number of cycles after which the step size
            is modified, default is 20

        Returns optimised parameter vector xopt
        """
        if tolerance is None:
            tolerance = 1e-6

        if len(xopt) == 0:
            return xopt

        random_series = random_series or random.Random()
        if seed is not None:
            random_series.seed(seed)

        schedule = AnnealingSchedule(
            temp_reduction, init_temp, temp_iterations, step_cycles
        )

        if self.restore and self.checkpointer.available():
            run = self.checkpointer.load()
            run.checkFunction(function, xopt, self.checkpointer.filename)
            run.schedule.checkSameConditions(schedule)
        else:
            run = AnnealingRun(function, xopt, schedule, random_series)
        self.restore = False

        result = run.run(
            function,
            tolerance,
            checkpointer=self.checkpointer,
            show_remaining=show_remaining,
        )

        return result.XOPT

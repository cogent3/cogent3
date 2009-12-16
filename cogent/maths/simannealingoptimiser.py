#!/usr/bin/env python
"""
Simulated annealing optimiser. Derives from basic optimiser class.

The simulated annealing optimiser is a translation into Python of the fortran
program simman.f authored by Bill Goffe (bgoffe@whale.st.usm.edu). The original
citation is "Global Optimization of Statistical Functions with Simulated
Annealing," Goffe, Ferrier and Rogers, Journal of Econometrics, vol. 60, no. 1/2,
Jan./Feb. 1994, pp. 65-100.
"""

from optimiser import OptimiserBase

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
import time

__author__ = "Andrew Butterfield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Andrew Butterfield", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class AnnealingSchedule(object):
    """Responsible for the shape of the simulated annealing temperature profile"""
    
    def __init__(self, temp_reduction, initial_temp, temp_iterations, step_cycles):
        if initial_temp < 0.0 :
            raise RuntimeError, "Initial temperature not +ve"
        self.T = self.initial_temp = initial_temp
        self.temp_reduction = temp_reduction
        self.temp_iterations = temp_iterations
        self.step_cycles = step_cycles
        self.dwell = temp_iterations * step_cycles
    
    def checkSameConditions(self, other):
        for attr in ['temp_reduction', 'initial_temp', 'temp_iterations', 'step_cycles']:
            if getattr(self, attr) != getattr(other, attr):
                raise ValueError('Checkpoint file ignored - %s different' % attr)
    
    def cool(self):
        self.T = self.temp_reduction * self.T
    
    def willAccept(self, newF, oldF, random_series):
        deltaF = newF - oldF
        return deltaF >= 0 or random_series.uniform(0.0, 1.0) < numpy.exp(deltaF / self.T)
    

class AnnealingHistory(object):
    """Keeps the last few results, for convergence testing"""
    
    def __init__(self, sample=4):
        self.values = [None] * sample
        self.i = 0
    
    def note(self, F):
        self.values[self.i] = F
        self.i = (self.i + 1) % len(self.values)
    
    def hasConverged(self, tolerance):
        return None not in self.values and max(self.values) - min(self.values) < tolerance
    

class AnnealingState(object):
    def __init__(self, X, function, random_series):
        self.random_series = random_series
        self.NFCNEV = 1
        self.VM = numpy.ones(len(X), Float)
        self.setX(X, function(X))
        (self.XOPT, self.FOPT) = (X, self.F)
        self.NACP = [0] * len(X)
        self.NTRY = 0
        self.elapsed_time = 0
    
    def setX(self, X, F):
        self.X = numpy.array(X, Float)
        self.F = F
    
    def step(self, function, accept_test):
        # One attempted move in each dimension
        t0 = time.time()
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
        self.elapsed_time += time.time() - t0
    
    def adjustStepSizes(self):
        # Adjust velocity in each dimension to keep acceptance ratios near 50%
        if self.NTRY == 0:
            return
        for I in range(len(self.X)):
            RATIO = (self.NACP[I]*1.0) / self.NTRY
            if RATIO > 0.6:
                self.VM[I] *= (1.0 + (2.0 * ((RATIO-0.6)/0.4)))
            elif RATIO < 0.4:
                self.VM[I] /= (1.0 + (2.0 * ((0.4 - RATIO)/0.4)))
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
                "Number of parameters in checkpoint file '%s' (%s) " \
                "don't match current function (%s)" % (
                    checkpointing_filename, len(self.state.XOPT), len(xopt)))
        # if f(x) != g(x) then f isn't g.
        then = self.state.FOPT
        now = function(self.state.XOPT)
        if not numpy.allclose(now, then, 1e-8):
            raise ValueError(
                "Function to optimise doesn't match checkpoint file " \
                "'%s': F=%s now, %s in file." % (
                    checkpointing_filename, now, then))
    
    def run(self, function, tolerance, max_iterations, checkpointer,
                show_progress):
        state = self.state
        history = self.history
        schedule = self.schedule
        
        while not history.hasConverged(tolerance):
            if show_progress:
                print "Outer loop = %d" % self.test_count
            
            self.save(checkpointer)
            
            for i in range(self.schedule.dwell):
                state.step(function, self.schedule.willAccept)
                self.test_count += 1
                if max_iterations and self.test_count >= max_iterations:
                    raise MaximumEvaluationsReached(state)
                if self.test_count % schedule.step_cycles == 0:
                    state.adjustStepSizes()
            
            history.note(state.F)
            if show_progress:
                print "\tF = %f EVALS = %s" % (state.FOPT, state.NFCNEV)
            state.setX(state.XOPT, state.FOPT)
            schedule.cool()
        
        self.save(checkpointer, final=True)
        
        return state
    
    def save(self, checkpointer, final=False):
        msg = "Number of function evaluations = %d; current F = %s" % \
                (self.state.NFCNEV, self.state.FOPT)
        checkpointer.record(self, msg, final)
    

class MaximumEvaluationsReached(Exception):
    # Used to pass out the results when iteration has to stop early
    """FORCED EXIT from SimulatedAnnealing:
Too many function evaluations, results are likely to be poor.
You can increase max_evaluations or decrease tolerance."""


class SimulatedAnnealing(OptimiserBase):
    """Simulated annealing optimiser for bounded functions
    """
    # this is a maximiser
    algorithm_direction = +1
    
    def _setdefaults(self):
        """set all the conditions for the sim annealing algorithm to default values"""
        self.setConditions(tolerance = 1E-6, temp_reduction = 0.5, init_temp=5.0,
                temp_iterations = 5, step_cycles = 20, max_evaluations=1e100)
    
    def setConditions(self, tolerance = None, temp_reduction = None, init_temp=None,
                temp_iterations = None, step_cycles = None, max_evaluations=None):
        """Set the conditions that control the optimisation.
        
        Arguments:
            - tolerance: the error condition for termination, default is
              1E-6
            - temp_reduction: the factor by which the annealing
              "temperature" is reduced, default is 0.5
            - temp_iterations: the number of iterations before a
              temperature reduction, default is 5
            - step_cycles: the number of cycles after which the step size
              is modified, default is 20
            - max_evaluations: the maximum number of function
              evaluations, default is 1E100. Note that a full run across
              the vector will be always be performed, with the outcome that
              the program will excape number of evaluations is greater than
              or equal to max_evaluations.
        """
        
        for (attr, value) in locals().items():
            if value is not None:
                setattr(self, attr, value)
    
    def runInner(self, function, xopt, show_progress, random_series):
        """Optimise the vector within the bounds specified by the base class.
        
        Arguments:
            - show_progress: whether the function values are printed as
              the optimisation proceeds. Default is True.
        
        Returns function value, parameter vector, evaluation count
        """
        if len(xopt) == 0:
            return function(xopt), xopt, 0, 0.0
        
        schedule = AnnealingSchedule(
            self.temp_reduction, self.init_temp, self.temp_iterations, self.step_cycles)
        max_iterations = (self.max_evaluations-1) / len(xopt) + 1
        
        if self.restore and self.checkpointer.available():
            run = self.checkpointer.load()
            run.checkFunction(function, xopt, self.checkpointer.filename)
            run.schedule.checkSameConditions(schedule)
        else:
            run = AnnealingRun(function, xopt, schedule, random_series)
        self.restore = False
        
        try:
            result = run.run(
                function,
                self.tolerance,
                max_iterations,
                checkpointer = self.checkpointer,
                show_progress = show_progress)
        except MaximumEvaluationsReached, detail:
            print detail.__doc__
            result = detail.args[0]
        return result.FOPT, result.XOPT, result.NFCNEV, result.elapsed_time
    

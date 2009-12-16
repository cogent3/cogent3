#!/usr/bin/env python
"""
Provides services for parametric bootstrapping. These include
the ability to estimate probabilities or estimate confidence intervals.

The estimation of probabilities is done by the EstimateProbability class.
Functions that provide ParameterController objects for the 'null' and
'alternative' cases are provided to the constructor. Numerous aspects of the
bootstrapping can be controlled such as the choice of numerical optimiser, and
the number of sample's from which to estimate the probability. This class
can be run in serial or in parallel (at the level of each random sample).

An observed Likelihood Ratio (LR) statistic is estimated using the provided
'observed' data. Random data sets are simulated under the null model and the
LR estimated from these. The probability of the observed LR is taken as the
number of sample LR's that were >= to the observed.

Confidence interval estimation can be done using the EstimateConfidenceIntervals
class. Multiple statistics associated with an analysis can be evaluated
simultaneously. Similar setup, and parallelisation options as provided by
the EstimateProbability class.

"""

from cogent.util import parallel
from cogent.maths import optimisers

import random

__author__ = "Gavin Huttley and Andrew Butterfield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley","Andrew Butterfield", "Matthew Wakefield",
                    "Edward Lang", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class ParametricBootstrapCore(object):
    """Core parametric bootstrap services."""
    
    def __init__(self):
        """Constructor for core parametric bootstrap services class."""
        self._numreplicates = 10
        self.seed = None
        self.results = []
    
    def setNumReplicates(self, num):
        self._numreplicates = num
    
    def setSeed(self, seed):
        self.seed = seed
    
    def run(self, show_progress=False, **opt_args):
        # Sets self.observed and self.results (a list _numreplicates long) to
        # whatever is returned from self.simplify([LF result from each PC]).
        # self.simplify() is used as the entire LF result might not be picklable
        # for MPI. Subclass must provide self.alignment and
        # self.parameter_controllers
        alignment_randomness = random.Random()
        alignment_randomness.seed(self.seed)
        parallel.sync_random(alignment_randomness)
        
        opt_args['show_progress'] = show_progress
        if 'random_series' not in opt_args and not opt_args.get('local', None):
            opt_args['random_series'] = random.Random()
        
        self._observed = []
        starting_points = []
        for pc in self.parameter_controllers:
            pc.setAlignment(self.alignment)
            lf = pc.optimise(return_calculator=True, **opt_args)
            starting_points.append(lf)
            self._observed.append(lf)
        self.observed = self.simplify(*self.parameter_controllers)

        (parallel_context, parallel_subcontext) = parallel.getSplitCommunicators(self._numreplicates)
        parallel.push(parallel_subcontext)
        try:
            # pc needs to be told that it can't use all the CPUs it did initially
            # which seems less than ideal (it should be able to discover that itself), 
            # but here's the 2 line workaround: 
            for pc in self.parameter_controllers:
                pc.setParamRule('parallel_context', is_const=True, value=parallel_subcontext)
            cpu_count = parallel_context.Get_size()
            this_cpu = parallel_context.Get_rank()
            local_results = []
            for i in range(self._numreplicates):
                # making the simulated alignment on every node keeps alignment_randomness in sync
                # across all cpus.  This is only a good way of doing it so long as making the alignment
                # is much faster than optimising the likelihood function.
                self.parameter_controllers[0].setAlignment(self.alignment)
                simalign = self.parameter_controllers[0].simulateAlignment(random_series=alignment_randomness)
                #print "simulated", simalign[:5]
                if i % cpu_count == this_cpu:
                    for (pc, starting_point) in zip(self.parameter_controllers, starting_points):
                            pc.real_par_controller.updateFromCalculator(starting_point) # reset
                            pc.setAlignment(simalign)
                            pc.optimise(**opt_args)
                    local_result = self.simplify(*self.parameter_controllers)
                    local_results.append(local_result)
        finally:
            parallel.pop(parallel_subcontext)
        
        split_results = parallel_context.allgather(local_results)
        
        self.results = []
        for results in split_results:
            self.results.extend(results)
    

class EstimateProbability(ParametricBootstrapCore):
    # 2 parameter controllers, LR
    
    def __init__(self, null_parameter_controller, alt_parameter_controller, alignment):
        ParametricBootstrapCore.__init__(self)
        self.alignment = alignment
        self.null_parameter_controller = null_parameter_controller
        self.alt_parameter_controller = alt_parameter_controller
        self.parameter_controllers =  [self.null_parameter_controller, self.alt_parameter_controller]
    
    def simplify(self, null_result, alt_result):
        return (null_result.getLogLikelihood(), alt_result.getLogLikelihood())
    
    def getObservedlnL(self):
        return self.observed
    
    def getSamplelnL(self):
        return self.results
    
    def getSampleLRList(self):
        LR = [2 * (alt_lnL - null_lnL) for (null_lnL, alt_lnL) in self.results]
        LR.sort()
        LR.reverse()
        return LR
    
    def getObservedLR(self):
        return 2 * (self.observed[1] - self.observed[0])
    
    def getEstimatedProb(self):
        """Return the estimated probability.
        
        Calculated as the number of sample LR's >= observed LR
        divided by the number of replicates.
        """
        
        observed_LR = self.getObservedLR()
        sample_LRs = self.getSampleLRList()
        
        for (count, value) in enumerate(sample_LRs):
            if value <= observed_LR:
                return float(count) / len(sample_LRs)
        return 1.0
    

class EstimateConfidenceIntervals(ParametricBootstrapCore):
    """Estimate confidence interval(s) for one or many statistics
    by parametric bootstrapping."""
    
    def __init__(self, parameter_controller, func_calcstats, alignment):
        # func_calcstats takes a param dict and returns the statistic of interest
        ParametricBootstrapCore.__init__(self)
        self.alignment = alignment
        self.parameter_controller = parameter_controller
        self.parameter_controllers =  [parameter_controller]
        self.func_calcstats = func_calcstats
    
    def simplify(self, result):
        return (result.getLogLikelihood(), self.func_calcstats(result)) #.getStatisticsAsDict())
    
    def getObservedStats(self):
        return self.observed[1]
    
    def getSampleStats(self):
        return [s for (lnL, s) in self.results]
    
    def getSamplelnL(self):
        return [lnL for (lnL, s) in self.results]
    
    def getObservedlnL(self):
        return self.observed[0]
    

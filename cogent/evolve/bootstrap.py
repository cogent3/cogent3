#!/usr/bin/env python
"""
Provides services for parametric bootstrapping. These include
the ability to estimate probabilities or estimate confidence intervals.

The estimation of probabilities is done by the EstimateProbability class.
Functions that provide ParameterController objects for the 'null' and
'alternative' cases are provided to the constructor. Numerous aspects of the
bootstrapping can be controlled such as the choice of numerical optimiser, and
the number of samples from which to estimate the probability. This class
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
from __future__ import with_statement, division
from cogent.util import parallel
from cogent.util import progress_display as UI

import random

__author__ = "Gavin Huttley, Andrew Butterfield and Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley","Andrew Butterfield", "Matthew Wakefield",
                    "Edward Lang", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
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
    
    @UI.display_wrap
    def run(self, ui, **opt_args):
        # Sets self.observed and self.results (a list _numreplicates long) to
        # whatever is returned from self.simplify([LF result from each PC]).
        # self.simplify() is used as the entire LF result might not be picklable
        # for MPI. Subclass must provide self.alignment and
        # self.parameter_controllers
        if 'random_series' not in opt_args and not opt_args.get('local', None):
            opt_args['random_series'] = random.Random()
        
        null_pc = self.parameter_controllers[0]
        pcs = len(self.parameter_controllers)
        if pcs == 1:
            model_label = ['']
        elif pcs == 2:
            model_label = ['null', 'alt ']
        else:
            model_label = ['null'] + ['alt%s'%i for i in range(1,pcs)]
        
        @UI.display_wrap
        def each_model(alignment, ui):
            def one_model(pc):
                pc.setAlignment(alignment)
                return pc.optimise(return_calculator=True, **opt_args)
            # This is not done in parallel because we depend on the side-
            # effect of changing the parameter_controller current values 
            memos = ui.eager_map(one_model, self.parameter_controllers, 
                    labels=model_label, pure=False)
            concise_result = self.simplify(*self.parameter_controllers)
            return (memos, concise_result)
        
        #optimisations = pcs * (self._numreplicates + 1)
        init_work = pcs / (self._numreplicates + pcs)
        ui.display('Original data', 0.0, init_work)
        (starting_points, self.observed) = each_model(self.alignment)
        
        ui.display('Randomness', init_work, 0.0)
        alignment_random_state = random.Random(self.seed).getstate()
        if self.seed is None:
            comm  = parallel.getCommunicator()
            alignment_random_state = comm.bcast(alignment_random_state, 0)
        
        def one_replicate(i):
            for (pc, start_point) in zip(self.parameter_controllers, starting_points):
                # may have fewer CPUs per replicate than for original
                pc.setupParallelContext()
                # using a calculator as a memo object to reset the params
                pc.updateFromCalculator(start_point)
            aln_rnd = random.Random(0)
            aln_rnd.setstate(alignment_random_state)
            aln_rnd.jumpahead(i*10**9)
            simalign = null_pc.simulateAlignment(random_series=aln_rnd)
            (dummy, result) = each_model(simalign)
            return result
        
        ui.display('Bootstrap', init_work)
        self.results = ui.eager_map(
                one_replicate, range(self._numreplicates), noun='replicate',
                start=init_work)


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
        return (result.getLogLikelihood(), self.func_calcstats(result))
    
    def getObservedStats(self):
        return self.observed[1]
    
    def getSampleStats(self):
        return [s for (lnL, s) in self.results]
    
    def getSamplelnL(self):
        return [lnL for (lnL, s) in self.results]
    
    def getObservedlnL(self):
        return self.observed[0]
    

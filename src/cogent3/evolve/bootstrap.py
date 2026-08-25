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

import random
from typing import Any

from scinexus.progress import Progress, get_progress


class ParametricBootstrapCore:
    """Core parametric bootstrap services."""

    def __init__(self) -> None:
        """Constructor for core parametric bootstrap services class."""
        self._numreplicates = 10
        self.seed = None
        self.results = []

    def set_num_replicates(self, num) -> None:
        self._numreplicates = num

    def set_seed(self, seed) -> None:
        self.seed = seed

    def run(
        self,
        show_progress: bool | Progress | dict[str, Any] = False,
        **opt_args,  # type: ignore[type-arg]
    ) -> None:
        # Sets self.observed and self.results (a list _numreplicates long) to
        # whatever is returned from self.simplify([LF result from each PC]).
        # self.simplify() is used as the entire LF result might not be picklable
        # for MPI. Subclass must provide self.alignment and
        # self.parameter_controllers
        if "random_series" not in opt_args and not opt_args.get("local"):
            opt_args["random_series"] = random.Random()

        null_pc = self.parameter_controllers[0]

        def each_model(alignment):
            def one_model(pc):
                pc.set_alignment(alignment)
                return pc.optimise(return_calculator=True, **opt_args)

            # This is not done in parallel because we depend on the side-
            # effect of changing the parameter_controller current values
            memos = [one_model(pc) for pc in self.parameter_controllers]
            concise_result = self.simplify(*self.parameter_controllers)
            return (memos, concise_result)

        progress = (
            get_progress(show_progress=True, **show_progress)
            if isinstance(show_progress, dict)
            else get_progress(show_progress)
        )
        ctx = progress.context(msg="Bootstrap")

        ctx.update(progress=0.0, msg="Original data")
        (starting_points, self.observed) = each_model(self.alignment)

        init_work = len(self.parameter_controllers) / (
            self._numreplicates + len(self.parameter_controllers)
        )
        ctx.update(progress=init_work, msg="Randomness")
        alignment_random_state = random.Random(self.seed).getstate()

        def one_replicate(i):
            for pc, start_point in zip(
                self.parameter_controllers,
                starting_points,
                strict=False,
            ):
                # may have fewer CPUs per replicate than for original
                # using a calculator as a memo object to reset the params
                pc.update_from_calculator(start_point)
            aln_rnd = random.Random(0)
            aln_rnd.setstate(alignment_random_state)
            # TODO jumpahead was deprecated, we need to consider an alternate
            # approach here. Commenting out for now.
            # `aln_rnd.jumpahead(i*10**9)`
            simalign = null_pc.simulate_alignment(random_series=aln_rnd)
            (dummy, result) = each_model(simalign)
            return result

        ctx.update(progress=init_work, msg="Bootstrap")
        replicates = range(self._numreplicates)
        self.results = [
            one_replicate(i)
            for i in progress(replicates, total=self._numreplicates, msg="replicates")
        ]
        ctx.close()


class EstimateProbability(ParametricBootstrapCore):
    # 2 parameter controllers, LR

    def __init__(
        self,
        null_parameter_controller,
        alt_parameter_controller,
        alignment,
    ) -> None:
        ParametricBootstrapCore.__init__(self)
        self.alignment = alignment
        self.null_parameter_controller = null_parameter_controller
        self.alt_parameter_controller = alt_parameter_controller
        self.parameter_controllers = [
            self.null_parameter_controller,
            self.alt_parameter_controller,
        ]

    def simplify(self, null_result, alt_result):
        return (null_result.get_log_likelihood(), alt_result.get_log_likelihood())

    def get_observed_lnL(self):
        return self.observed

    def get_sample_lnL(self):
        return self.results

    def get_sample_LR_list(self):
        LR = [2 * (alt_lnL - null_lnL) for (null_lnL, alt_lnL) in self.results]
        LR.sort()
        LR.reverse()
        return LR

    def get_observed_LR(self):
        return 2 * (self.observed[1] - self.observed[0])

    def get_estimated_prob(self):
        """Return the estimated probability.

        Calculated as the number of sample LR's >= observed LR
        divided by the number of replicates.
        """

        observed_LR = self.get_observed_LR()
        sample_LRs = self.get_sample_LR_list()

        for count, value in enumerate(sample_LRs):
            if value <= observed_LR:
                return float(count) / len(sample_LRs)
        return 1.0


class EstimateConfidenceIntervals(ParametricBootstrapCore):
    """Estimate confidence interval(s) for one or many statistics
    by parametric bootstrapping."""

    def __init__(self, parameter_controller, func_calcstats, alignment) -> None:
        # func_calcstats takes a param dict and returns the statistic of
        # interest
        ParametricBootstrapCore.__init__(self)
        self.alignment = alignment
        self.parameter_controller = parameter_controller
        self.parameter_controllers = [parameter_controller]
        self.func_calcstats = func_calcstats

    def simplify(self, result):
        return (result.get_log_likelihood(), self.func_calcstats(result))

    def get_observed_stats(self):
        return self.observed[1]

    def getSampleStats(self):
        return [s for (lnL, s) in self.results]

    def get_sample_lnL(self):
        return [lnL for (lnL, s) in self.results]

    def get_observed_lnL(self):
        return self.observed[0]

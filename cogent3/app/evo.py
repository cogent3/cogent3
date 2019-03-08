import os
import json
from tqdm import trange, tqdm
from cogent3 import LoadTree
from cogent3.evolve.models import get_model
from .composable import ComposableHypothesis, ComposableModel, ErrorResult
from .result import hypothesis_result, model_result, bootstrap_result
from cogent3.util import parallel
from cogent3.util import progress_display as UI


class model(ComposableModel):
    def __init__(self, sm, tree=None, name=None, sm_args=None,
                 lf_args=None, time_heterogeneity=None, param_rules=None,
                 opt_args=None, split_codons=False, show_progress=False):
        """represents a substitution model + tree"""
        super(model, self).__init__(input_type='aligned',
                                    output_type=(
                                        'model_result', 'serialisable'))
        self._formatted_params()
        sm_args = sm_args or {}
        if type(sm) == str:
            sm = get_model(sm, **sm_args)
        self._sm = sm
        if len(sm.get_motifs()[0]) > 1:
            split_codons = False
        self._tree = tree
        self._lf_args = lf_args or {}
        if not name:
            name = sm.name or 'unnamed model'
        self.name = name
        self._opt_args = opt_args or dict(max_restarts=5,
                                          show_progress=show_progress)
        self._opt_args['show_progress'] = self._opt_args.get('show_progress',
                                                             show_progress)
        param_rules = param_rules or {}
        if param_rules:
            for rule in param_rules:
                rule['upper'] = rule.get('upper', 50)  # default upper bound
        self._param_rules = param_rules
        self._time_heterogeneity = time_heterogeneity
        self._split_codons = split_codons
        self.func = self.fit

    def _configure_lf(self, aln, identifier, initialise=None):
        lf = self._sm.make_likelihood_function(self._tree, **self._lf_args)
        # lf.setup_parallel_context()

        lf.set_alignment(aln)
        verbose = self._opt_args.get('show_progress', False)
        if self._param_rules:
            lf.apply_param_rules(self._param_rules)
        elif self._time_heterogeneity == 'max':
            if not initialise:
                if verbose:
                    print("Time homogeneous fit..")

                # we opt with a time-homogeneous process first
                opt_args = self._opt_args.copy()
                opt_args.update(dict(max_restart=1, tolerance=1e-3))
                lf.optimise(**self._opt_args)
                if verbose:
                    print(lf)
            lf.set_time_heterogeneity(is_independent=True, upper=50)
        else:
            rules = lf.get_param_rules()
            for rule in rules:
                if rule['par_name'] != 'mprobs':
                    rule['upper'] = rule.get('upper', 50)

            lf.apply_param_rules([rule])

        if initialise:
            initialise(lf, identifier)

        self._lf = lf

    def _fit_aln(self, aln, identifier=None, initialise=None, construct=True,
                 **opt_args):
        if construct:
            self._configure_lf(aln=aln, identifier=identifier,
                               initialise=initialise)
        lf = self._lf
        kwargs = self._opt_args.copy()
        kwargs.update(opt_args)

        verbose = kwargs.get('show_progress', False)
        if verbose:
            print("Fit...")

        calc = lf.optimise(return_calculator=True, **kwargs)
        lf.calculator = calc

        if verbose:
            print(lf)

        return lf

    def fit(self, aln, initialise=None, construct=True, **opt_args):
        evaluation_limit = opt_args.get('max_evaluations', None)
        if self._tree is None:
            assert len(aln.names) == 3
            self._tree = LoadTree(tip_names=aln.names)

        result = model_result(name=self.name, stat=sum,
                              source=aln.info.source,
                              evaluation_limit=evaluation_limit)
        if not self._split_codons:
            lf = self._fit_aln(aln, initialise=initialise, construct=construct,
                               **opt_args)
            result[self.name] = lf
            result.num_evaluations = lf.calculator.evaluations
            result.elapsed_time = lf.calculator.elapsed_time
        else:
            num_evals = 0
            elapsed_time = 0
            for i in range(3):
                codon_pos = aln[i::3]
                lf = self._fit_aln(codon_pos, identifier=str(i + 1),
                                   initialise=initialise, construct=construct,
                                   **opt_args)
                result[i + 1] = lf
                num_evals += lf.calculator.evaluations
                elapsed_time += lf.calculator.elapsed_time

            result.num_evaluations = num_evals
            result.elapsed_time = elapsed_time

        return result


class hypothesis(ComposableHypothesis):
    def __init__(self, null, *alternates):
        super(hypothesis, self).__init__(input_type='aligned',
                                         output_type=('result', 'serialisable'))
        self._formatted_params()
        self.null = null
        self._alts = alternates
        self.func = self.test_hypothesis

    def _initialised_alt_from_null(self, null, aln):
        def init(alt):
            try:
                alt.initialise_from_nested(null)
            except:
                pass
            return alt

        results = []
        for alt in self._alts:
            result = alt(aln, initialise=init)
            results.append(result)
        return results

    def test_hypothesis(self, aln):
        try:
            null = self.null(aln)
        except ValueError as err:
            msg = "Hypothesis null had bounds error %s" % aln.info.source
            raise ValueError(msg)
        try:
            alts = [alt for alt in self._initialised_alt_from_null(null, aln)]
        except ValueError as err:
            msg = "Hypothesis alt had bounds error %s" % aln.info.source
            raise ValueError(msg)
        results = {alt.name: alt for alt in alts}
        results.update({null.name: null})

        result = hypothesis_result(name_of_null=null.name,
                                   source=aln.info.source)
        result.update(results)
        return result


class bootstrap(ComposableHypothesis):
    def __init__(self, hyp, num_reps, verbose=False):
        super(bootstrap, self).__init__(input_type='aligned',
                                        output_type=('result', 'serialisable'))
        self._formatted_params()
        self._hyp = hyp
        self._num_reps = num_reps
        self._verbose = verbose
        self.func = self.run

    def _fit_sim(self, rep_num):
        sim_aln = self._null.simulate_alignment()
        sim_aln.info.source = "%s - simalign %d" % (self._inpath, rep_num)

        try:
            sym_result = self._hyp(sim_aln)
        except ValueError:
            sym_result = None
        return sym_result

    def run(self, aln):
        result = bootstrap_result(aln.info.source)
        try:
            obs = self._hyp(aln)
        except ValueError as err:
            result = ErrorResult('ERROR', str(self._hyp), err.args[0])
            return result
        result.observed = obs
        self._null = obs.null
        self._inpath = aln.info.source

        sym_results = [r for r in
                       parallel.imap(self._fit_sim, range(self._num_reps)) if r]
        for sym_result in sym_results:
            if not sym_result:
                continue

            result.add_to_null(sym_result)

        return result

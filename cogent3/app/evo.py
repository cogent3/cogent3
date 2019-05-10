from tqdm import tqdm
from cogent3 import LoadTree
from cogent3.evolve.models import get_model
from .composable import ComposableHypothesis, ComposableModel, NotCompletedResult
from .result import hypothesis_result, model_result, bootstrap_result
from cogent3.util import parallel


class model(ComposableModel):
    """represents a substitution model + tree"""

    def __init__(self, sm, tree=None, name=None, sm_args=None,
                 lf_args=None, time_het=None, param_rules=None,
                 opt_args=None, split_codons=False, show_progress=False,
                 verbose=False):
        """

        Parameters
        ----------
        sm : str or instance
            substitution model if string must be available via get_model()
        tree : str, Tree instance, or None
            if None, assumes a star phylogeny (only valid for 3 taxa)
        name
            name of the model
        sm_args
            arguments to be passed to the substitution model constructor, e.g.
            dict(optimise_motif_probs=True)
        lf_args
            arguments to be passed to the likelihood function constructor
        time_het
            'max' or a list of dicts corresponding to edge_sets, e.g.
            [dict(edges=['Human', 'Chimp'], is_independent=False, upper=10)].
            Passed to the likelihood function .set_time_heterogeneity()
            method.
        param_rules
            other parameter rules, passed to the likelihood function
            set_param_rule() method
        opt_args
            arguments for the numerical optimiser, e.g.
            dict(max_restarts=5, tolerance=1e-6, max_evaluations=1000,
            limit_action='ignore')
        split_codons : bool
            if True, incoming alignments are split into the 3 frames and each
            frame is fit separately
        show_progress : bool
            show progress bars during numerical optimisation
        verbose : bool
            prints intermediate states to screen during fitting

        Returns
        -------
        Calling an instance with an alignment returns a model_result instance
        with the optimised likelihood function. In the case of split_codons,
        the result object has a separate entry for each.
        """
        super(model, self).__init__(input_type='aligned',
                                    output_type=(
                                        'model_result', 'serialisable'))
        self._verbose = verbose
        self._formatted_params()
        sm_args = sm_args or {}
        if type(sm) == str:
            sm = get_model(sm, **sm_args)
        self._sm = sm
        if len(sm.get_motifs()[0]) > 1:
            split_codons = False

        if type(tree) == str:
            tree = LoadTree(treestring=tree)

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
        self._time_het = time_het
        self._split_codons = split_codons
        self.func = self.fit

    def _configure_lf(self, aln, identifier, initialise=None):
        lf = self._sm.make_likelihood_function(self._tree, **self._lf_args)

        lf.set_alignment(aln)
        if self._param_rules:
            lf.apply_param_rules(self._param_rules)
        elif self._time_het:
            if not initialise:
                if self._verbose:
                    print("Time homogeneous fit..")

                # we opt with a time-homogeneous process first
                opt_args = self._opt_args.copy()
                opt_args.update(dict(max_restart=1, tolerance=1e-3))
                lf.optimise(**self._opt_args)
                if self._verbose:
                    print(lf)
            if self._time_het == 'max':
                lf.set_time_heterogeneity(is_independent=True, upper=50)
            else:
                lf.set_time_heterogeneity(edge_sets=self._time_het)
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

        if self._verbose:
            print("Fit...")

        calc = lf.optimise(return_calculator=True, **kwargs)
        lf.calculator = calc

        if identifier:
            lf.set_name(f'LF id: {identifier}')

        if self._verbose:
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
    def __init__(self, null, *alternates, init_alt=None):
        # todo document! init_alt needs to be able to take null, alt and *args
        super(hypothesis, self).__init__(input_type='aligned',
                                         output_type=('result', 'serialisable'))
        self._formatted_params()
        self.null = null
        names = {a.name for a in alternates}
        names.add(null.name)
        if len(names) != len(alternates) + 1:
            msg = f"{names} model names not unique"
            raise ValueError(msg)

        self._alts = alternates
        self.func = self.test_hypothesis
        self._init_alt = init_alt

    def _initialised_alt_from_null(self, null, aln):
        def init(alt, *args, **kwargs):
            try:
                alt.initialise_from_nested(null.lf)
            except:
                pass
            return alt

        if callable(self._init_alt):
            init_func = self._init_alt(null)
        else:
            init_func = init

        results = []
        for alt in self._alts:
            result = alt(aln, initialise=init_func)
            results.append(result)
        return results

    def test_hypothesis(self, aln):
        try:
            null = self.null(aln)
        except ValueError as err:
            msg = f"Hypothesis null had bounds error {aln.info.source}"
            return NotCompletedResult('ERROR', self, msg, source=aln)
        try:
            alts = [
                alt for alt in self._initialised_alt_from_null(null, aln)]
        except ValueError as err:
            msg = f"Hypothesis alt had bounds error {aln.info.source}"
            return NotCompletedResult('ERROR', self, msg, source=aln)
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
            result = NotCompletedResult('ERROR', str(self._hyp), err.args[0])
            return result
        result.observed = obs
        self._null = obs.null
        self._inpath = aln.info.source

        sym_results = [r for r in
                       parallel.imap(self._fit_sim, range(self._num_reps)) if r]
        for sym_result in tqdm(sym_results):
            if not sym_result:
                continue

            result.add_to_null(sym_result)

        return result

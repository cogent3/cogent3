import os

from tqdm import tqdm

from cogent3 import LoadTree
from cogent3.evolve.models import get_model
from cogent3.util import parallel

from .composable import (
    ComposableHypothesis,
    ComposableModel,
    ComposableTabular,
    NotCompleted,
)
from .result import (
    bootstrap_result,
    hypothesis_result,
    model_result,
    tabular_result,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class model(ComposableModel):
    """Define a substitution model + tree for maximum likelihood evaluation.
    Returns model_result."""

    _input_type = frozenset(["aligned"])
    _output_type = frozenset(["result", "model_result", "serialisable"])
    _data_types = frozenset(["ArrayAlignment", "Alignment"])

    def __init__(
        self,
        sm,
        tree=None,
        name=None,
        sm_args=None,
        lf_args=None,
        time_het=None,
        param_rules=None,
        opt_args=None,
        split_codons=False,
        show_progress=False,
        verbose=False,
    ):
        """

        Parameters
        ----------
        sm : str or instance
            substitution model if string must be available via get_model()
        tree
            if None, assumes a star phylogeny (only valid for 3 taxa). Can be a
            newick formatted tree, a path to a file containing one, or a Tree
            instance.
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
        super(model, self).__init__()
        self._verbose = verbose
        self._formatted_params()
        sm_args = sm_args or {}
        if type(sm) == str:
            sm = get_model(sm, **sm_args)
        self._sm = sm
        if len(sm.get_motifs()[0]) > 1:
            split_codons = False

        if type(tree) == str:
            if os.path.exists(tree):
                kwargs = dict(filename=tree, underscore_unmunge=True)
            else:
                kwargs = dict(treestring=tree, underscore_unmunge=True)
            tree = LoadTree(**kwargs)

        self._tree = tree
        self._lf_args = lf_args or {}
        if not name:
            name = sm.name or "unnamed model"
        self.name = name
        self._opt_args = opt_args or dict(max_restarts=5, show_progress=show_progress)
        self._opt_args["show_progress"] = self._opt_args.get(
            "show_progress", show_progress
        )
        param_rules = param_rules or {}
        if param_rules:
            for rule in param_rules:
                if rule.get("is_constant"):
                    continue
                rule["upper"] = rule.get("upper", 50)  # default upper bound
        self._param_rules = param_rules
        self._time_het = time_het
        self._split_codons = split_codons
        self.func = self.fit

    def _configure_lf(self, aln, identifier, initialise=None):
        lf = self._sm.make_likelihood_function(self._tree, **self._lf_args)

        lf.set_alignment(aln)
        if self._param_rules:
            lf.apply_param_rules(self._param_rules)

        if self._time_het:
            if not initialise:
                if self._verbose:
                    print("Time homogeneous fit..")

                # we opt with a time-homogeneous process first
                opt_args = self._opt_args.copy()
                opt_args.update(dict(max_restart=1, tolerance=1e-3))
                lf.optimise(**self._opt_args)
                if self._verbose:
                    print(lf)
            if self._time_het == "max":
                lf.set_time_heterogeneity(is_independent=True, upper=50)
            else:
                lf.set_time_heterogeneity(edge_sets=self._time_het)
        else:
            rules = lf.get_param_rules()
            for rule in rules:
                if rule["par_name"] not in ("mprobs", "psubs"):
                    rule["upper"] = rule.get("upper", 50)

            lf.apply_param_rules([rule])

        if initialise:
            initialise(lf, identifier)

        self._lf = lf

    def _fit_aln(
        self, aln, identifier=None, initialise=None, construct=True, **opt_args
    ):
        if construct:
            self._configure_lf(aln=aln, identifier=identifier, initialise=initialise)
        lf = self._lf
        kwargs = self._opt_args.copy()
        kwargs.update(opt_args)

        if self._verbose:
            print("Fit...")

        calc = lf.optimise(return_calculator=True, **kwargs)
        lf.calculator = calc

        if identifier:
            lf.set_name(f"LF id: {identifier}")

        if self._verbose:
            print(lf)

        return lf

    def fit(self, aln, initialise=None, construct=True, **opt_args):
        evaluation_limit = opt_args.get("max_evaluations", None)
        if self._tree is None:
            assert len(aln.names) == 3
            self._tree = LoadTree(tip_names=aln.names)

        result = model_result(
            name=self.name,
            stat=sum,
            source=aln.info.source,
            evaluation_limit=evaluation_limit,
        )
        if not self._split_codons:
            lf = self._fit_aln(
                aln, initialise=initialise, construct=construct, **opt_args
            )
            result[self.name] = lf
            result.num_evaluations = lf.calculator.evaluations
            result.elapsed_time = lf.calculator.elapsed_time
        else:
            num_evals = 0
            elapsed_time = 0
            for i in range(3):
                codon_pos = aln[i::3]
                lf = self._fit_aln(
                    codon_pos,
                    identifier=str(i + 1),
                    initialise=initialise,
                    construct=construct,
                    **opt_args,
                )
                result[i + 1] = lf
                num_evals += lf.calculator.evaluations
                elapsed_time += lf.calculator.elapsed_time

            result.num_evaluations = num_evals
            result.elapsed_time = elapsed_time

        return result


class hypothesis(ComposableHypothesis):
    """Specify a hypothesis through defining two models. Returns a
    hypothesis_result."""

    _input_type = frozenset(["aligned"])
    _output_type = frozenset(["result", "hypothesis_result", "serialisable"])
    _data_types = frozenset(["ArrayAlignment", "Alignment"])

    def __init__(self, null, *alternates, init_alt=None):
        # todo document! init_alt needs to be able to take null, alt and *args
        super(hypothesis, self).__init__()
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
            return NotCompleted("ERROR", self, msg, source=aln)
        try:
            alts = [alt for alt in self._initialised_alt_from_null(null, aln)]
        except ValueError as err:
            msg = f"Hypothesis alt had bounds error {aln.info.source}"
            return NotCompleted("ERROR", self, msg, source=aln)
        results = {alt.name: alt for alt in alts}
        results.update({null.name: null})

        result = hypothesis_result(name_of_null=null.name, source=aln.info.source)
        result.update(results)
        return result


class bootstrap(ComposableHypothesis):
    """Parametric bootstrap for a provided hypothesis. Returns a bootstrap_result."""

    _input_type = frozenset(["aligned"])
    _output_type = frozenset(["result", "bootstrap_result", "serialisable"])
    _data_types = frozenset(["ArrayAlignment", "Alignment"])

    def __init__(self, hyp, num_reps, verbose=False):
        super(bootstrap, self).__init__()
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
            result = NotCompleted("ERROR", str(self._hyp), err.args[0])
            return result
        result.observed = obs
        self._null = obs.null
        self._inpath = aln.info.source

        sym_results = [
            r for r in parallel.imap(self._fit_sim, range(self._num_reps)) if r
        ]
        for sym_result in tqdm(sym_results):
            if not sym_result:
                continue

            result.add_to_null(sym_result)

        return result


class ancestral_states(ComposableTabular):
    _input_type = frozenset(["model_result"])
    _output_type = frozenset(["result", "tabular_result", "serialisable"])
    _data_types = frozenset(["model_result"])

    def __init__(self):
        super(ancestral_states, self).__init__()
        self._formatted_params()
        self.func = self.recon_ancestor

    def recon_ancestor(self, result):
        """returns a tabular_result of posterior probabilities of ancestral states"""
        anc = result.lf.reconstruct_ancestral_seqs()
        fl = result.lf.get_full_length_likelihoods()
        template = None
        tab = tabular_result(source=result.source)
        for name in anc:
            state_p = anc[name]
            if template is None:
                template = state_p.template
            pp = (state_p.array.T / fl).T
            tab[name] = template.wrap(pp)
        return tab


class tabulate_stats(ComposableTabular):
    """Extracts all model statistics from model_result as Table."""

    _input_type = frozenset(["model_result"])
    _output_type = frozenset(["result", "tabular_result", "serialisable"])
    _data_types = frozenset(["model_result"])

    def __init__(self):
        super(tabulate_stats, self).__init__()
        self._formatted_params()
        self.func = self.extract_stats

    def extract_stats(self, result):
        """returns Table for all statistics returned by likelihood function
        get_statistics"""
        stats = result.lf.get_statistics(with_titles=True, with_motif_probs=True)
        tab = tabular_result(source=result.source)
        for table in stats:
            tab[table.title] = table
        return tab

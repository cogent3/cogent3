#!/usr/bin/env python

import random
from collections import defaultdict
import numpy

from cogent3.core.alignment import ArrayAlignment
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.evolve.simulate import AlignmentEvolver, random_sequence
from cogent3.util import parallel, table
from cogent3.util.misc import adjusted_gt_minprob
from cogent3.recalculation.definition import ParameterController
from cogent3.maths.matrix_logarithm import is_generator_unique
from cogent3.maths.matrix_exponential_integration import expected_number_subs
from cogent3.evolve import substitution_model

from cogent3.util.warning import discontinued, deprecated

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Andrew Butterfield", "Peter Maxwell",
               "Matthew Wakefield", "Rob Knight", "Brett Easton",
               "Ben Kaehler", "Ananias Iliadis"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# cogent3.evolve.parameter_controller.LikelihoodParameterController tells the
# recalculation framework to use this subclass rather than the generic
# recalculation Calculator.  It adds methods which are useful for examining
# the parameter, psub, mprob and likelihood values after the optimisation is
# complete.


def _get_param_mapping(rich, simple):
    """returns {simple_param_name: rich_param_name, ...}, the mapping of simple
    to rich parameters based on matrix coordinates
    """
    assert len(rich) >= len(simple)
    simple_to_rich = defaultdict(set)
    rich_to_simple = defaultdict(set)
    for simple_param in simple:
        simple_coords = simple[simple_param]
        for rich_param in rich:
            rich_coords = rich[rich_param]
            if rich_coords <= simple_coords:
                simple_to_rich[simple_param].add(rich_param)
                rich_to_simple[rich_param].add(simple_param)

    for rich_param in rich_to_simple:
        simple_counterparts = rich_to_simple[rich_param]
        if len(simple_counterparts) == 1:
            continue

        sized_simple = [(len(simple[param]), param) for param in simple_counterparts]
        sized_simple.sort()
        if sized_simple[0][0] == sized_simple[1][0]:
            msg = "%s and %s tied for matrix space" % (sized_simple[0][1],
                                                       sized_simple[1][1])
            raise ValueError(msg)

        _, chosen = sized_simple.pop(0)
        rich_to_simple[rich_param] = [chosen]
        for _, simple_param in sized_simple:
            simple_to_rich[simple_param].pop(rich_param)

    return simple_to_rich


class _ParamProjection:
    """projects parameter names, values between nested models"""

    def __init__(self, simple_model, rich_model, motif_probs, same=True):
        # construct following by calling the functions we wrote
        self._rich_coords = rich_model.get_param_matrix_coords(include_ref_cell=True)
        self._simple_coords = simple_model.get_param_matrix_coords(include_ref_cell=True)
        self._param_map = _get_param_mapping(self._rich_coords, self._simple_coords)
        self._same = same
        # end of constructing attributes
        self._motif_probs = motif_probs
        self._ref_val = self._set_ref_val(same)
        self.projected_rate = {False: self._rate_not_same}.get(same,
                                                               self._rate_same)

    def _set_ref_val(self, same):
        """returns the motif prob corresponding to the model reference cell"""
        if same:
            return 1
        else:
            i, j = list(self._rich_coords['ref_cell'])[0]
            return self._motif_probs[j]

    def _rate_not_same(self, simple_param, mle):
        """returns {rich_param: val, ...} from simple_param: val"""
        ref_val = self._ref_val
        new_terms = {}
        for rich_param in self._param_map[simple_param]:
            if rich_param == 'ref_cell':
                continue
            for i, j in self._rich_coords[rich_param]:
                new_terms[rich_param] = self._motif_probs[j] * mle / ref_val

        return new_terms

    def _rate_same(self, simple_param, mle):
        new_terms = {}
        for rich_param in self._param_map[simple_param]:
            if rich_param == 'ref_cell':
                continue
            for i, j in self._rich_coords[rich_param]:
                new_terms[rich_param] = mle
        return new_terms

    def update_param_rules(self, rules):
        new_rules = []
        if not self._same:
            rules = rules[:] + [dict(par_name='ref_cell', init=1.0, edges=None)]

        for rule in rules:
            # get the param name, mle, call self.projected_rate
            name = rule['par_name']
            if name in ('mprobs', 'length'):
                new_rules.append(rule)
                continue

            if rule.get('is_constant', False):
                par_val_key = 'value'
            else:
                par_val_key = 'init'

            mle = rule[par_val_key]

            proj_rate = self.projected_rate(name, mle)
            for new_name, new_mle in proj_rate.items():
                rule_dict = rule.copy()
                rule_dict['par_name'] = new_name
                # update it with the new parname and mle and append to new rules
                rule_dict['init'] = new_mle
                new_rules.append(rule_dict)

        return new_rules


def compatible_likelihood_functions(lf1, lf2):
    """returns True if all attributes of the two likelihood functions are compatible
    for mapping parameters, else raises ValueError or an AssertionError"""
    # tree's must have the same topology AND be oriented the same way
    # plus have the same edge names
    if len(lf1.bin_names) != 1 or len(lf1.bin_names) != len(lf2.bin_names):
        raise NotImplementedError("Too many bins")
    if len(lf1.locus_names) != 1 or len(lf1.locus_names) != len(lf2.locus_names):
        raise NotImplementedError("Too many loci")
    if lf1.model.get_motifs() != lf2.model.get_motifs():
        raise AssertionError("Motifs don't match")
    if lf1.tree.get_newick(with_node_names=True) != lf2.tree.get_newick(with_node_names=True):
        raise AssertionError("Topology, Orientation or node names don't match")
    return True


class LikelihoodFunction(ParameterController):

    def get_log_likelihood(self):
        return self.get_final_result()

    def get_psub_for_edge(self, name, **kw):
        """returns the substitution probability matrix for the named edge"""
        try:
            # For PartialyDiscretePsubsDefn
            array = self.get_param_value('dpsubs', edge=name, **kw)
        except KeyError:
            array = self.get_param_value('psubs', edge=name, **kw)
        return DictArrayTemplate(self._motifs, self._motifs).wrap(array)

    def get_rate_matrix_for_edge(self, name, calibrated=True, **kw):
        """returns the rate matrix (Q) for the named edge

        If calibrated=False, expm(Q) will give the same result as
        get_psub_for_edge(name)"""
        try:
            array = self.get_param_value('Q', edge=name, **kw)
            array = array.copy()
            if not calibrated:
                length = self.get_param_value('length', edge=name, **kw)
                array *= length
        except KeyError as err:
            if err[0] == 'Q' and name != 'Q':
                raise RuntimeError('rate matrix not known by this model')
            else:
                raise
        return DictArrayTemplate(self._motifs, self._motifs).wrap(array)

    def _getLikelihoodValuesSummedAcrossAnyBins(self, locus=None):
        if self.bin_names and len(self.bin_names) > 1:
            root_lhs = [self.get_param_value('lh', locus=locus, bin=bin) for
                        bin in self.bin_names]
            bprobs = self.get_param_value('bprobs')
            root_lh = bprobs.dot(root_lhs)
        else:
            root_lh = self.get_param_value('lh', locus=locus)
        return root_lh

    def get_full_length_likelihoods(self, locus=None):
        """Array of [site, motif] likelihoods from the root of the tree"""
        root_lh = self._getLikelihoodValuesSummedAcrossAnyBins(locus=locus)
        root_lht = self.get_param_value('root', locus=locus)
        return root_lht.get_full_length_likelihoods(root_lh)

    def get_G_statistic(self, return_table=False, locus=None):
        """Goodness-of-fit statistic derived from the unambiguous columns"""
        root_lh = self._getLikelihoodValuesSummedAcrossAnyBins(locus=locus)
        root_lht = self.get_param_value('root', locus=locus)
        return root_lht.calc_G_statistic(root_lh, return_table)

    def reconstruct_ancestral_seqs(self, locus=None):
        """returns a dict of DictArray objects containing probabilities
        of each alphabet state for each node in the tree.

        Arguments:
            - locus: a named locus"""
        result = {}
        array_template = None
        for restricted_edge in self._tree.get_edge_vector():
            if restricted_edge.istip():
                continue
            try:
                r = []
                for motif in range(len(self._motifs)):
                    self.set_param_rule('fixed_motif', value=motif,
                                      edge=restricted_edge.name, locus=locus,
                                      is_constant=True)
                    likelihoods = self.get_full_length_likelihoods(locus=locus)
                    r.append(likelihoods)
                    if array_template is None:
                        array_template = DictArrayTemplate(
                            likelihoods.shape[0], self._motifs)
            finally:
                self.set_param_rule('fixed_motif', value=-1,
                                  edge=restricted_edge.name, locus=locus,
                                  is_constant=True)
            # dict of site x motif arrays
            result[restricted_edge.name] = array_template.wrap(
                numpy.transpose(numpy.asarray(r)))
        return result

    def likely_ancestral_seqs(self, locus=None):
        """Returns the most likely reconstructed ancestral sequences as an
        alignment.

        Arguments:
            - locus: a named locus"""
        prob_array = self.reconstruct_ancestral_seqs(locus=locus)
        seqs = []
        for edge, probs in list(prob_array.items()):
            seq = []
            for row in probs:
                by_p = [(p, state) for state, p in list(row.items())]
                seq.append(max(by_p)[1])
            seqs += [(edge, self.model.moltype.make_seq("".join(seq)))]
        return ArrayAlignment(data=seqs, moltype=self.model.moltype)

    def get_bin_probs(self, locus=None):
        hmm = self.get_param_value('bindex', locus=locus)
        lhs = [self.get_param_value('lh', locus=locus, bin=bin)
               for bin in self.bin_names]
        array = hmm.get_posterior_probs(*lhs)
        return DictArrayTemplate(self.bin_names, array.shape[1]).wrap(array)

    def _valuesForDimension(self, dim):
        # in support of __str__
        if dim == 'edge':
            result = [e.name for e in self._tree.get_edge_vector()]
        elif dim == 'bin':
            result = self.bin_names[:]
        elif dim == 'locus':
            result = self.locus_names[:]
        elif dim.startswith('motif'):
            result = self._mprob_motifs
        elif dim == 'position':
            result = self.posn_names[:]
        else:
            raise KeyError(dim)
        return result

    def _valuesForDimensions(self, dims):
        # in support of __str__
        result = [[]]
        for dim in dims:
            new_result = []
            for r in result:
                for cat in self._valuesForDimension(dim):
                    new_result.append(r + [cat])
            result = new_result
        return result

    def _for_display(self):
        """processes statistics tables for display"""
        title = self._name if self._name else 'Likelihood function statistics'
        result = []
        result += self.get_statistics(with_motif_probs=True, with_titles=True)
        for i, table in enumerate(result):
            if 'motif' in table.title and table.shape[1] == 2 and\
               table.shape[0] >= 60: # just sort codon motif probs, then truncate
                table = table.sorted(columns="motif")
                data = table.tolist()
                data = data[:5] + [["...", "..."]] + data[-5:]
                table = table.__class__(header=table.header, rows=data,
                                        digits=table._digits, title=table.title)
                result[i] = table
        return title, result

    def _repr_html_(self):
        """for jupyter notebook display"""
        def row_cell_func(val, col, row):
            val = '<td style="font-family: monospace, monospace">%s</td>' % val
            return val

        try:
            lnL = "<p>log-likelihood = %.4f</p>" % self.get_log_likelihood()
        except ValueError:
            # alignment probably not yet set
            lnL = ""

        nfp = "<p>number of free parameters = %d</p>" % \
            self.get_num_free_params()
        title, results = self._for_display()
        for i, table in enumerate(results):
            table.title = table.title.capitalize()
            results[i] = table.to_rich_html(row_cell_func=row_cell_func)
        results = ["<h4>%s</h4>" % title, lnL, nfp] + results
        return "\n".join(results)

    def __str__(self):
        title, results = self._for_display()

        try:
            lnL = "log-likelihood = %.4f" % self.get_log_likelihood()
        except ValueError:
            # alignment probably not yet set
            lnL = None

        nfp = "number of free parameters = %d" % self.get_num_free_params()
        for table in results:
            table.title = ""
        
        if lnL:
            results = [title, lnL, nfp] + results
        else:
            results = [title, nfp] + results

        return '\n'.join(map(str, results))

    def get_annotated_tree(self, length_as_ens=False):
        """returns tree with model attributes on node.params
        
        length_as_ens : bool
            replaces 'length' param with expected number of substition,
            which will be different to standard length if the substition
            model is non-stationary
        """
        d = self.get_param_value_dict(['edge'])
        if length_as_ens:
            lengths = self.get_lengths_as_ens()
        tree = self._tree.deepcopy()
        for edge in tree.get_edge_vector():
            if edge.name == 'root':
                continue
            for par in d:
                val = d[par][edge.name]
                if length_as_ens:
                    val = lengths[edge.name]
                edge.params[par] = val

        return tree

    def get_motif_probs(self, edge=None, bin=None, locus=None):
        motif_probs_array = self.get_param_value(
            'mprobs', edge=edge, bin=bin, locus=locus)
        # these can fall below the minimum allowed value due to
        # rounding errors, so I adjust these
        motif_probs_array = adjusted_gt_minprob(motif_probs_array,
                                                minprob=1e-6)
        return DictArrayTemplate(self._mprob_motifs).wrap(motif_probs_array)

    def get_bin_prior_probs(self, locus=None):
        bin_probs_array = self.get_param_value('bprobs', locus=locus)
        return DictArrayTemplate(self.bin_names).wrap(bin_probs_array)

    def get_scaled_lengths(self, predicate, bin=None, locus=None):
        """A dictionary of {scale:{edge:length}}"""
        if not hasattr(self._model, 'get_scaled_lengthsFromQ'):
            return {}

        def valueOf(param, **kw):
            return self.get_param_value(param, locus=locus, **kw)

        if bin is None:
            bin_names = self.bin_names
        else:
            bin_names = [bin]

        if len(bin_names) == 1:
            bprobs = [1.0]
        else:
            bprobs = valueOf('bprobs')

        mprobs = [valueOf('mprobs', bin=b) for b in bin_names]

        scaled_lengths = {}
        for edge in self._tree.get_edge_vector():
            if edge.isroot():
                continue
            Qs = [valueOf('Qd', bin=b, edge=edge.name).Q for b in bin_names]
            length = valueOf('length', edge=edge.name)
            scaled_lengths[edge.name] = length * self._model.get_scale_from_Qs(
                Qs, bprobs, mprobs, predicate)
        return scaled_lengths

    def get_lengths_as_ens(self):
        """returns {edge: ens, ...} where ens is the expected number of substitutions
        
        for a stationary Markov process, this is just branch length"""
        node_names = self.tree.get_node_names()
        node_names.remove('root')
        lengths = {e: self.get_param_value('length', edge=e) for e in node_names}
        if not isinstance(self.model, substitution_model.Stationary):
            ens = {}
            for e in node_names:
                mprobs = self.get_motif_probs(edge=e)
                Q = self.get_rate_matrix_for_edge(e)
                length = expected_number_subs(mprobs, Q, lengths[e])
                ens[e] = length

            lengths = ens

        return lengths

    def get_param_rules(self):
        """returns the [{rule}, ..] that would allow reconstruction"""
        def mprob_rule(defn, index, edges):
            uniq = defn.uniq[index]
            value = uniq.value
            # default min prob is 1e-6,
            # values can drop below due to precision
            # so we adjust
            value = adjusted_gt_minprob(value, minprob=1e-6)
            val = dict(zip(defn.bin_names, value))
            rule = dict(par_name=defn.name, edges=edges)
            if uniq.is_constant:
                rule.update(dict(is_constant=True, value=val))
            else:
                rule.update(dict(init=val))
            return rule

        rules = []
        
        # markov model rate terms
        param_names = self.get_param_names()
        rate_names = self.model.get_param_list()
        for param_name in param_names:
            defn = self.defn_for[param_name]
            scoped = defaultdict(list)
            for key, index in defn.index.items():
                edge_name = key[0]
                scoped[index].append(edge_name)
                
            if len(scoped) == 1:  # we have a global
                if param_name == 'mprobs':
                    rule = mprob_rule(defn, 0, None)
                else:
                    val = defn.values[0]
                    rule = dict(par_name=param_name, edges=None)
                    if defn.uniq[0].is_constant:
                        rule.update(dict(is_constant=True, value=val))
                    else:
                        rule.update(dict(init=val, upper=defn.upper,
                                         lower=defn.lower))

                rules.append(rule)

                continue
            
            for index in scoped:
                edges = scoped[index]
                if param_name == 'mprobs':
                    rule = mprob_rule(defn, index, edges)
                else:
                    uniq = defn.uniq[index]
                    rule = dict(par_name=param_name, edges=edges)
                    if defn.uniq[index].is_constant:
                        rule.update(dict(is_constant=True, value=uniq.value))
                    else:
                        rule.update(dict(init=uniq.value, lower=uniq.lower,
                                         upper=uniq.upper))
                rules.append(rule)

        return rules


    def get_statistics(self, with_motif_probs=True, with_titles=True):
        """returns the parameter values as tables/dict

        Arguments:
            - with_motif_probs: include the motif probability table
            - with_titles: include a title for each table based on it's
              dimension"""
        result = []
        group = {}
        param_names = self.get_param_names()

        mprob_name = [n for n in param_names if 'mprob' in n]
        if mprob_name:
            mprob_name = mprob_name[0]
        else:
            mprob_name = ''

        if not with_motif_probs:
            param_names.remove(mprob_name)

        for param in param_names:
            dims = tuple(self.get_used_dimensions(param))
            if dims not in group:
                group[dims] = []
            group[dims].append(param)
        table_order = list(group.keys())
        table_order.sort()
        for table_dims in table_order:
            raw_table = self.get_param_value_dict(
                dimensions=table_dims, params=group[table_dims])
            param_names = group[table_dims]
            param_names.sort()
            if table_dims == ('edge',):
                if 'length' in param_names:
                    param_names.remove('length')
                    param_names.insert(0, 'length')
                raw_table['parent'] = dict([(e.name, e.parent.name)
                                            for e in self._tree.get_edge_vector()
                                            if not e.isroot()])
                param_names.insert(0, 'parent')
            list_table = []
            heading_names = list(table_dims) + param_names
            row_order = self._valuesForDimensions(table_dims)
            for scope in row_order:
                row = {}
                row_used = False
                for param in param_names:
                    d = raw_table[param]
                    try:
                        for part in scope:
                            d = d[part]
                    except KeyError:
                        d = 'NA'
                    else:
                        row_used = True
                    row[param] = d
                if row_used:
                    row.update(dict(list(zip(table_dims, scope))))
                    row = [row[k] for k in heading_names]
                    list_table.append(row)
            if table_dims:
                title = ['', '%s params' % ' '.join(table_dims)][with_titles]
                row_ids = True
            else:
                row_ids = False
                title = ['', 'global params'][with_titles]
            result.append(table.Table(
                heading_names, list_table,
                max_width=80, row_ids=row_ids,
                title=title, **self._format))
        return result

    # For tests.  Compat with old LF interface
    def set_name(self, name):
        self._name = name

    def get_name(self):
        return self._name or 'unnamed'

    def set_tables_format(self, space=4, digits=4):
        """sets display properties for statistics tables. This affects results
        of str(lf) too."""
        space = [space, 4][type(space) != int]
        digits = [digits, 4][type(digits) != int]
        self._format = dict(space=space, digits=digits)

    def get_motif_probs_by_node(self, edges=None, bin=None, locus=None):
        kw = dict(bin=bin, locus=locus)
        mprobs = self.get_param_value('mprobs', **kw)
        mprobs = self._model.calc_word_probs(mprobs)
        result = self._nodeMotifProbs(self._tree, mprobs, kw)
        if edges is None:
            edges = [name for (name, m) in result]
        result = dict(result)
        values = [result[name] for name in edges]
        return DictArrayTemplate(edges, self._mprob_motifs).wrap(values)

    def _nodeMotifProbs(self, tree, mprobs, kw):
        result = [(tree.name, mprobs)]
        for child in tree.children:
            psub = self.get_psub_for_edge(child.name, **kw)
            child_mprobs = numpy.dot(mprobs, psub)
            result.extend(self._nodeMotifProbs(child, child_mprobs, kw))
        return result

    def simulate_alignment(self, sequence_length=None, random_series=None,
                          exclude_internal=True, locus=None, seed=None, root_sequence=None):
        """
        Returns an alignment of simulated sequences with key's corresponding to
        names from the current attached alignment.

        Arguments:
            - sequence_length: the legnth of the alignment to be simulated,
              default is the length of the attached alignment.
            - random_series: a random number generator.
            - exclude_internal: if True, only sequences for tips are returned.
            - root_sequence: a sequence from which all others evolve.
        """

        if sequence_length is None:
            lht = self.get_param_value('lht', locus=locus)
            sequence_length = len(lht.index)
            leaves = self.get_param_value('leaf_likelihoods', locus=locus)
            orig_ambig = {}
            for (seq_name, leaf) in list(leaves.items()):
                orig_ambig[seq_name] = leaf.get_ambiguous_positions()
        else:
            orig_ambig = {}

        if random_series is None:
            random_series = random.Random()
            random_series.seed(seed)
            parallel.sync_random(random_series)

        def psub_for(edge, bin):
            return self.get_psub_for_edge(edge, bin=bin, locus=locus)

        if len(self.bin_names) > 1:
            hmm = self.get_param_value('bdist', locus=locus)
            site_bins = hmm.emit(sequence_length, random_series)
        else:
            site_bins = numpy.zeros([sequence_length], int)

        evolver = AlignmentEvolver(random_series, orig_ambig, exclude_internal,
                                   self.bin_names, site_bins, psub_for, self._motifs)

        if root_sequence is not None:  # we convert to a vector of motifs
            if isinstance(root_sequence, str):
                root_sequence = self._model.moltype.make_seq(root_sequence)
            motif_len = self._model.get_alphabet().get_motif_len()
            root_sequence = root_sequence.get_in_motif_size(motif_len)
        else:
            mprobs = self.get_param_value('mprobs', locus=locus, edge='root')
            mprobs = self._model.calc_word_probs(mprobs)
            mprobs = dict((m, p) for (m, p) in zip(self._motifs, mprobs))
            root_sequence = random_sequence(
                random_series, mprobs, sequence_length)

        simulated_sequences = evolver(self._tree, root_sequence)

        return ArrayAlignment(
            data=simulated_sequences,
            moltype=self._model.moltype)

    def all_psubs_DLC(self):
        """Returns True if every Psub matrix is Diagonal Largest in Column"""
        for edge in self.tree.get_edge_vector(include_root=False):
            P = self.get_psub_for_edge(edge.name).asarray()
            if (P.diagonal() < P).any():
                return False
        return True

    def all_rate_matrices_unique(self):
        """Returns True if every rate matrix is unique for its Psub matrix"""
        for edge in self.tree.get_edge_vector(include_root=False):
            Q = self.get_rate_matrix_for_edge(edge.name).asarray()
            t = self.get_param_value('length', edge=edge.name)
            if not is_generator_unique(Q * t):
                return False
        return True

    def initialise_from_nested(self, nested_lf, length_kwargs=None, param_kwargs=None):
        from cogent3.evolve.substitution_model import TimeReversible
        assert self.get_num_free_params() > nested_lf.get_num_free_params(), "wrong order for likelihood functions"
        compatible_likelihood_functions(self, nested_lf)

        same = (isinstance(self.model, TimeReversible) and isinstance(nested_lf.model, TimeReversible))\
            or (not isinstance(self.model, TimeReversible) and not isinstance(nested_lf.model, TimeReversible))

        if length_kwargs is None:
            length_kwargs = {}
        if param_kwargs is None:
            param_kwargs = {}

        mprobs = nested_lf.get_motif_probs()
        edge_names = self.tree.get_node_names()
        edge_names.remove('root')
        lengths = {name: nested_lf.get_param_value('length', edge=name) for name in edge_names}
        param_proj = _ParamProjection(nested_lf.model, self.model, mprobs, same=same)
        param_rules = nested_lf.get_param_rules()
        param_rules = param_proj.update_param_rules(param_rules)

        with self.updates_postponed():
            for rule in param_rules:
                rule.update(param_kwargs)
                self.set_param_rule(**rule)

        return self

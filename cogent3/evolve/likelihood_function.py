#!/usr/bin/env python

import random
import numpy

from cogent3.core.alignment import Alignment
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.evolve.simulate import AlignmentEvolver, randomSequence
from cogent3.util import parallel, table
from cogent3.recalculation.definition import ParameterController
from cogent3.maths.matrix_logarithm import is_generator_unique

from cogent3.util.warning import discontinued, deprecated

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Andrew Butterfield", "Peter Maxwell",
               "Matthew Wakefield", "Rob Knight", "Brett Easton",
               "Ben Kaehler"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# cogent3.evolve.parameter_controller.LikelihoodParameterController tells the
# recalculation framework to use this subclass rather than the generic
# recalculation Calculator.  It adds methods which are useful for examining
# the parameter, psub, mprob and likelihood values after the optimisation is
# complete.


class LikelihoodFunction(ParameterController):

    def setpar(self, param_name, value, edge=None, **scope):
        deprecated('method', 'setpar', 'setParamRule', '1.6')
        return self.setParamRule(param_name, edge=edge, value=value, is_constant=True, **scope)

    def get_log_likelihood(self):
        return self.getFinalResult()

    def get_psub_for_edge(self, name, **kw):
        """returns the substitution probability matrix for the named edge"""
        try:
            # For PartialyDiscretePsubsDefn
            array = self.get_param_value('dpsubs', edge=name, **kw)
        except KeyError:
            array = self.get_param_value('psubs', edge=name, **kw)
        return DictArrayTemplate(self._motifs, self._motifs).wrap(array)

    def get_rate_matrix_for_edge(self, name, **kw):
        """returns the rate matrix (Q) for the named edge

        Note: expm(Q) will give the same result as get_psub_for_edge(name)"""
        try:
            array = self.get_param_value('Q', edge=name, **kw)
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
        return root_lht.calcGStatistic(root_lh, return_table)

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
                    self.setParamRule('fixed_motif', value=motif,
                                      edge=restricted_edge.name, locus=locus,
                                      is_constant=True)
                    likelihoods = self.get_full_length_likelihoods(locus=locus)
                    r.append(likelihoods)
                    if array_template is None:
                        array_template = DictArrayTemplate(
                            likelihoods.shape[0], self._motifs)
            finally:
                self.setParamRule('fixed_motif', value=-1,
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
            seqs += [(edge, self.model.moltype.make_sequence("".join(seq)))]
        return Alignment(data=seqs, moltype=self.model.moltype)

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

    def __str__(self):
        if not self._name:
            title = 'Likelihood Function Table'
        else:
            title = self._name
        result = [title]
        result += self.getStatistics(with_motif_probs=True, with_titles=False)
        return '\n'.join(map(str, result))

    def getAnnotatedTree(self):
        d = self.get_param_valueDict(['edge'])
        tree = self._tree.deepcopy()
        for edge in tree.get_edge_vector():
            if edge.name == 'root':
                continue
            for par in d:
                edge.params[par] = d[par][edge.name]
        return tree

    def getMotifProbs(self, edge=None, bin=None, locus=None):
        motif_probs_array = self.get_param_value(
            'mprobs', edge=edge, bin=bin, locus=locus)
        return DictArrayTemplate(self._mprob_motifs).wrap(motif_probs_array)
        # return dict(zip(self._motifs, motif_probs_array))

    def getBinPriorProbs(self, locus=None):
        bin_probs_array = self.get_param_value('bprobs', locus=locus)
        return DictArrayTemplate(self.bin_names).wrap(bin_probs_array)

    def getScaledLengths(self, predicate, bin=None, locus=None):
        """A dictionary of {scale:{edge:length}}"""
        if not hasattr(self._model, 'getScaledLengthsFromQ'):
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
            scaled_lengths[edge.name] = length * self._model.getScaleFromQs(
                Qs, bprobs, mprobs, predicate)
        return scaled_lengths

    def getStatistics(self, with_motif_probs=True, with_titles=True):
        """returns the parameter values as tables/dict

        Arguments:
            - with_motif_probs: include the motif probability table
            - with_titles: include a title for each table based on it's
              dimension"""
        result = []
        group = {}
        param_names = self.getParamNames()

        mprob_name = [n for n in param_names if 'mprob' in n]
        if mprob_name:
            mprob_name = mprob_name[0]
        else:
            mprob_name = ''

        if not with_motif_probs:
            param_names.remove(mprob_name)

        for param in param_names:
            dims = tuple(self.getUsedDimensions(param))
            if dims not in group:
                group[dims] = []
            group[dims].append(param)
        table_order = list(group.keys())
        table_order.sort()
        for table_dims in table_order:
            raw_table = self.get_param_valueDict(
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

    def getStatisticsAsDict(self, with_parent_names=True,
                            with_edge_names=False):
        """Returns a dictionary containing the statistics for each edge of the
        tree, and any other information provided by the substitution model. The
        dictionary is keyed at the top-level by parameter name, and then by
        edge.name.

        Arguments:
            - with_edge_names: if True, an ordered list of edge names is
              included under the top-level key 'edge.names'. Default is
              False.
        """

        discontinued('method', "'getStatisticsAsDict' "
                     "use 'get_param_valueDict(['edge'])' is nearly equivalent",
                     '1.6')

        stats_dict = self.get_param_valueDict(['edge'])

        if hasattr(self.model, 'scale_masks'):
            for predicate in self.model.scale_masks:
                stats_dict[predicate] = self.getScaledLengths(predicate)

        edge_vector = [e for e in self._tree.get_edge_vector() if not e.isroot()]

        # do the edge names
        if with_parent_names:
            parents = {}
            for edge in edge_vector:
                if edge.parent.isroot():
                    parents[edge.name] = "root"
                else:
                    parents[edge.name] = str(edge.parent.name)
            stats_dict["edge.parent"] = parents

        if with_edge_names:
            stats_dict['edge.name'] = (
                [e.name for e in edge_vector if e.istip()] +
                [e.name for e in edge_vector if not e.istip()])

        return stats_dict

    # For tests.  Compat with old LF interface
    def setName(self, name):
        self._name = name

    def get_name(self):
        return self._name or 'unnamed'

    def setTablesFormat(self, space=4, digits=4):
        """sets display properties for statistics tables. This affects results
        of str(lf) too."""
        space = [space, 4][type(space) != int]
        digits = [digits, 4][type(digits) != int]
        self._format = dict(space=space, digits=digits)

    def getMotifProbsByNode(self, edges=None, bin=None, locus=None):
        kw = dict(bin=bin, locus=locus)
        mprobs = self.get_param_value('mprobs', **kw)
        mprobs = self._model.calcWordProbs(mprobs)
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

    def simulateAlignment(self, sequence_length=None, random_series=None,
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
            orig_ambig = {}  # alignment.getPerSequenceAmbiguousPositions()
            for (seq_name, leaf) in list(leaves.items()):
                orig_ambig[seq_name] = leaf.getAmbiguousPositions()
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
                root_sequence = self._model.moltype.make_sequence(root_sequence)
            motif_len = self._model.getAlphabet().get_motif_len()
            root_sequence = root_sequence.get_in_motif_size(motif_len)
        else:
            mprobs = self.get_param_value('mprobs', locus=locus, edge='root')
            mprobs = self._model.calcWordProbs(mprobs)
            mprobs = dict((m, p) for (m, p) in zip(self._motifs, mprobs))
            root_sequence = randomSequence(
                random_series, mprobs, sequence_length)

        simulated_sequences = evolver(self._tree, root_sequence)

        return Alignment(
            data=simulated_sequences,
            moltype=self._model.moltype)

    def allPsubsDLC(self):
        """Returns True if every Psub matrix is Diagonal Largest in Column"""
        for edge in self.tree.get_edge_vector(include_root=False):
            P = self.get_psub_for_edge(edge.name).asarray()
            if (P.diagonal() < P).any():
                return False
        return True

    def allRateMatricesUnique(self):
        """Returns True if every rate matrix is unique for its Psub matrix"""
        for edge in self.tree.get_edge_vector(include_root=False):
            Q = self.get_rate_matrix_for_edge(edge.name).asarray()
            t = self.get_param_value('length', edge=edge.name)
            if not is_generator_unique(Q * t):
                return False
        return True

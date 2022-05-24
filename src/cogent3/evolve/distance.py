#!/usr/bin/env python
"""Estimating pairwise distances between sequences.
"""
from itertools import combinations
from warnings import warn

from cogent3 import make_tree
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.maths.stats.number import NumberCounter
from cogent3.util import progress_display as UI
from cogent3.util import table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def get_name_combinations(names, group_size):
    """returns combinations of names"""
    combined = list(tuple(sorted(p)) for p in combinations(names, group_size))
    combined.sort()
    return combined


def get_pairwise_distance_from_triad(data, summary_function="mean"):
    """returns pairwise distances from lengths estimated from triads

    Parameters
    ----------
    data
        a dict keyed as {(a,b,c)
    summary_function
        a string naming the function used for
        estimating param from threeway distances. Valid values are 'mean'
        (default) and 'median'.

    """
    summary_func = summary_function.lower()
    pairwise_stats = {}
    lengths = {}
    for key in data:
        a, b, c = key
        for x, y in [(a, b), (a, c), (b, c)]:
            length = data[key]["length"][x] + data[key]["length"][y]
            try:
                lengths[(x, y)].append(length)
            except KeyError:
                lengths[(x, y)] = [length]

    # get all the distances involving this pair
    for pair in lengths:
        values = NumberCounter(lengths[pair])
        pairwise_stats[pair] = getattr(values, summary_func)

    return pairwise_stats


class EstimateDistances(object):
    """base class used for estimating pairwise distances between sequences.
    Can also estimate other parameters from pairs."""

    def __init__(
        self,
        seqs,
        submodel,
        threeway=False,
        motif_probs=None,
        do_pair_align=False,
        rigorous_align=False,
        est_params=None,
        modify_lf=None,
    ):
        """Arguments:
            - seqs: an Alignment or SeqCollection instance with > 1 sequence
            - submodel: substitution model object Predefined models can
              be imported from cogent3.evolve.models
            - threeway: a boolean flag for using threeway comparisons to
              estimate distances. default False. Ignored if do_pair_align is
              True.
            - do_pair_align: if the input sequences are to be pairwise aligned
              first and then the distance will be estimated. A pair HMM based
              on the submodel will be used.
            - rigorous_align: if True the pairwise alignments are actually
              numerically optimised, otherwise the current substitution model
              settings are used. This slows down estimation considerably.
            - est_params: substitution model parameters to save estimates from
              in addition to length (distance)
            - modify_lf: a callback function for that takes a likelihood
              function (with alignment set) and modifies it. Can be used to
              configure local_params, set bounds, optimise using a restriction
              for faster performance.

        Note: Unless you know a priori your alignment will be flush ended
        (meaning no sequence has terminal gaps) it is advisable to construct a
        substitution model that recodes gaps. Otherwise the terminal gaps will
        significantly bias the estimation of branch lengths when using
        do_pair_align.
        """

        if do_pair_align:
            self._threeway = False
        else:
            # whether pairwise is to be estimated from 3-way
            self._threeway = [threeway, False][do_pair_align]

        self._seq_collection = seqs
        self._seqnames = seqs.names[:]
        self._motif_probs = motif_probs

        # the following may be pairs or three way combinations
        self._combination_aligns = None
        self._do_pair_align = do_pair_align
        self._rigorous_align = rigorous_align

        # substitution model stuff
        self._sm = submodel
        self._modify_lf = modify_lf

        # store for the results
        self._param_ests = {}
        self._est_params = list(est_params or [])

        self._run = False  # a flag indicating whether estimation completed

    def __str__(self):
        return str(self.get_table())

    def _make_pair_alignment(self, seqs, opt_kwargs):
        lf = self._sm.make_likelihood_function(
            make_tree(tip_names=seqs.names), aligned=False
        )
        lf.set_sequences(seqs.named_seqs)

        # allow user to modify the lf config
        if self._modify_lf:
            lf = self._modify_lf(lf)

        if self._rigorous_align:
            lf.optimise(**opt_kwargs)
        lnL = lf.get_log_likelihood()
        return lnL.edge.get_viterbi_path().get_alignment()

    @UI.display_wrap
    def _doset(self, sequence_names, dist_opt_args, aln_opt_args, ui):
        # slice the alignment
        seqs = self._seq_collection.take_seqs(sequence_names)
        if self._do_pair_align:
            ui.display("Aligning", progress=0.0)
            align = self._make_pair_alignment(seqs, aln_opt_args)
            ui.display("", progress=0.5)

        else:
            align = seqs
            ui.display("", progress=0.0)
        # note that we may want to consider removing the redundant gaps

        # create the tree object
        tree = make_tree(tip_names=sequence_names)

        # make the parameter controller
        lf = self._sm.make_likelihood_function(tree)
        if not self._threeway:
            lf.set_param_rule("length", is_independent=False)

        if self._motif_probs:
            lf.set_motif_probs(self._motif_probs)

        lf.set_alignment(align)

        # allow user modification of lf using the modify_lf
        if self._modify_lf:
            lf = self._modify_lf(lf)

        lf.optimise(**dist_opt_args)

        # get the statistics
        stats_dict = lf.get_param_value_dict(
            ["edge"], params=["length"] + self._est_params
        )

        # if two-way, grab first distance only
        if not self._threeway:
            result = {"length": list(stats_dict["length"].values())[0] * 2.0}
        else:
            result = {"length": stats_dict["length"]}

        # include any other params requested
        for param in self._est_params:
            result[param] = list(stats_dict[param].values())[0]

        return result

    @UI.display_wrap
    def run(self, dist_opt_args=None, aln_opt_args=None, ui=None, **kwargs):
        """Start estimating the distances between sequences. Distance estimation
        is done using the Powell local optimiser. This can be changed using the
        dist_opt_args and aln_opt_args.

        Parameters
        ----------
        show_progress
            whether to display progress. More detailed progress
            information from individual optimisation is controlled by the
            ..opt_args.
        dist_opt_args, aln_opt_args
            arguments for the optimise method for
            the distance estimation and alignment estimation respectively.

        """

        if "local" in kwargs:
            warn(
                "local argument ignored, provide it to dist_opt_args or"
                " aln_opt_args",
                DeprecationWarning,
                stacklevel=2,
            )

        ui.display("Distances")
        dist_opt_args = dist_opt_args or {}
        aln_opt_args = aln_opt_args or {}
        # set the optimiser defaults
        dist_opt_args["local"] = dist_opt_args.get("local", True)
        aln_opt_args["local"] = aln_opt_args.get("local", True)
        # generate the list of unique sequence sets (pairs or triples) to be
        # analysed
        if self._threeway:
            combination_aligns = get_name_combinations(self._seq_collection.names, 3)
            desc = "triplet "
        else:
            combination_aligns = get_name_combinations(self._seq_collection.names, 2)
            desc = "pair "
        labels = [desc + ",".join(names) for names in combination_aligns]

        def _one_alignment(comp):
            result = self._doset(comp, dist_opt_args, aln_opt_args)
            return (comp, result)

        for (comp, value) in ui.imap(_one_alignment, combination_aligns, labels=labels):
            self._param_ests[comp] = value

    def get_pairwise_param(self, param, summary_function="mean"):
        """Return the pairwise statistic estimates as a dictionary keyed by
        (seq1, seq2)

        Parameters
        ----------
        param
            name of a parameter in est_params or 'length'
        summary_function
            a string naming the function used for
            estimating param from threeway distances. Valid values are 'mean'
            (default) and 'median'.

        """
        pairwise_stats = {}
        assert param in self._est_params + ["length"], f"unrecognised param {param}"
        if not self._param_ests:
            return None

        if self._threeway and param == "length":
            pairwise_stats = get_pairwise_distance_from_triad(
                self._param_ests, summary_function=summary_function
            )
        else:
            # no additional processing of the distances is required
            for comp_names, param_vals in list(self._param_ests.items()):
                pairwise_stats[comp_names] = param_vals[param]

        return pairwise_stats

    def get_pairwise_distances(self, summary_function="mean", **kwargs):
        """Return the pairwise distances as a dictionary keyed by (seq1, seq2).
        Convenience interface to get_pairwise_param.

        Parameters
        ----------
        summary_function
            a string naming the function used for
            estimating param from threeway distances. Valid values are 'mean'
            (default) and 'median'.

        """
        dists = self.get_pairwise_param(
            "length", summary_function=summary_function, **kwargs
        )
        return None if not dists else DistanceMatrix(dists)

    def get_param_values(self, param, **kwargs):
        """Returns a Numbers object with all estimated values of param.

        Parameters
        ----------
        param
            name of a parameter in est_params or 'length'
        **kwargs
            arguments passed to get_pairwise_param

        """
        ests = self.get_pairwise_param(param, **kwargs)
        return NumberCounter(list(ests.values()))

    def get_all_param_values(self):
        """returns raw estimated parameter dictionary"""
        return self._param_ests.copy()

    def get_table(self, summary_function="mean", **kwargs):
        """returns a Table instance of the distance matrix.

        Parameters
        ----------
        summary_function
            a string naming the function used for
            estimating param from threeway distances. Valid values are 'mean'
            (default) and 'median'.

        """
        d = self.get_pairwise_distances(summary_function=summary_function, **kwargs)
        if not d:
            d = {}
            for s1 in self._seqnames:
                for s2 in self._seqnames:
                    if s1 == s2:
                        continue
                    else:
                        d[(s1, s2)] = "Not Done"
        twoD = []
        for s1 in self._seqnames:
            row = [s1]
            for s2 in self._seqnames:
                if s1 == s2:
                    row.append("")
                    continue
                try:
                    row.append(d[(s1, s2)])
                except KeyError:
                    row.append(d[(s2, s1)])
            twoD.append(row)
        return table.Table(
            [r"Seq1 \ Seq2"] + self._seqnames,
            twoD,
            index_name=r"Seq1 \ Seq2",
            missing_data="*",
        )

    def get_newick_trees(self):
        """Returns a list of Newick format trees for supertree methods."""
        trees = []
        for comp_names, param_vals in list(self._param_ests.items()):
            tips = []
            for name in comp_names:
                tips.append(repr(name) + f":{param_vals[name]}")
            trees.append("(" + ",".join(tips) + ");")

        return trees

    def write(self, filename, summary_function="mean", format="phylip", **kwargs):
        """Save the pairwise distances to a file using phylip format. Other
        formats can be obtained by getting to a Table.

        Parameters
        ----------
        filename
            where distances will be written, required.
        summary_function
            a string naming the function used for
            estimating param from threeway distances. Valid values are 'mean'
            (default) and 'median'.
        format
            output format of distance matrix

        """

        table = self.get_table(summary_function=summary_function, **kwargs)
        table.write(filename, format=format)

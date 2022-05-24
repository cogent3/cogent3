from numpy import exp, log

from cogent3.util.io import atomic_write

from . import consensus


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"


class _UserList(list):
    def __getitem__(self, index):
        # Helpful to keep type after truncation like [self[:10]],
        # but not self[0] or self[x,y,-1]
        result = list.__getitem__(self, index)
        if isinstance(index, slice) and index.step is None:
            result = type(self)(result)
        return result


class ScoredTreeCollection(_UserList):
    """An ordered list of (score, tree) tuples"""

    def write(self, filename):
        with atomic_write(filename, mode="wt") as f:
            for (score, tree) in self:
                f.writelines(
                    self.scored_tree_format(
                        tree.get_newick(with_distances=True), str(score)
                    )
                )

    def scored_tree_format(self, tree, score):
        return [tree, "\t[", score, "]\n"]

    def get_consensus_tree(self, strict=None, method="unrooted"):
        ctrees = self.get_consensus_trees(strict, method=method)
        assert len(ctrees) == 1, len(ctrees)
        return ctrees[0]

    def get_consensus_trees(self, strict=True, method="unrooted"):
        if strict is None:
            strict = True
        return consensus.weighted_majority_rule(self, strict, method=method)


class UsefullyScoredTreeCollection(ScoredTreeCollection):
    def scored_tree_format(self, tree, score):
        return [score, "\t", tree, "\n"]


class WeightedTreeCollection(UsefullyScoredTreeCollection):
    """An ordered list of (weight, tree) tuples"""

    def get_consensus_trees(self, strict=False, method="unrooted"):
        if strict is None:
            strict = False
        return consensus.weighted_majority_rule(self, strict, method=method)


class LogLikelihoodScoredTreeCollection(UsefullyScoredTreeCollection):
    """An ordered list of (log likelihood, tree) tuples"""

    def __init__(self, trees):
        list.__init__(self, trees)
        # Quick and very dirty check of order
        assert self[0][0] >= self[-1][0]

    def get_consensus_tree(self, cutoff=None, strict=False, alpha=0.05):
        """See documentation for get_consensus_trees"""
        return self.get_consensus_trees(cutoff, strict, alpha)[0]

    def get_consensus_trees(self, cutoff=None, strict=False, alpha=0.05):
        """Returns a weighted consensus tree as described in Holland (2006).
        Weights transformed according to the class IV transformation in Jermiin
        (1997).

        Args:
            cutoff: Discard trees with weights that sum to cutoff.
            strict: Discard splits with consensus weight <= 0.5.
            alpha: Significance level, see Jermiiin (1997).

        Returns:
            Single consensus tree, in a list for backward compatibility.

        Holland, B. R., Jermiin, L. S., Moulton, V., and
        Investigators, S. T.-N. Y. (2006). Proceedings of the SMBE
        Tri-National Young Investigators' Workshop 2005. Improved consensus
        network techniques for genome-scale phylogeny. Mol Biol Evol,
        23(5), 848-855. doi:10.1093/molbev/msj061

        Jermiin, L. S., Olsen, G. J., Mengerson, K., and Easteal, S. (1997).
        Majority-rule consensus of phylogenetic trees obtained by
        maximum-likelihood analysis. Molecular Biology and Evolution,
        14(12):1296.
        """
        return self.get_weighted_trees(cutoff, alpha).get_consensus_trees(strict)

    def get_weighted_trees(self, cutoff=None, alpha=0.05):
        if cutoff is None:
            cutoff = 0.99
        assert 0 <= cutoff <= 1.0
        max_lnL = self[0][0]
        forgotten = log(alpha) / (self[-1][0] - max_lnL)
        weights = [exp(forgotten * (lnL - max_lnL)) for (lnL, t) in self]
        # add from smallest end to avoid rounding errors
        weights.reverse()
        tail = (1.0 - cutoff) * sum(weights)
        dropped = 0.0
        for (index, weight) in enumerate(weights):
            dropped += weight
            if dropped > tail:
                weights = weights[index:]
                break
        denominator = sum(weights)
        weights.reverse()
        return WeightedTreeCollection(
            (weight / denominator, tree) for (weight, (lnL, tree)) in zip(weights, self)
        )


def make_trees(filename):
    """Parse a file of (score, tree) lines. Scores can be positive probabilities
    or negative log likelihoods."""
    from cogent3 import make_tree

    infile = open(filename, "r")
    trees = []
    klass = list
    # expect score, tree
    for line in infile:
        line = line.split(None, 1)
        lnL = float(line[0])
        if lnL > 1:
            raise ValueError(f"likelihoods expected, not {lnL}")
        elif lnL > 0:
            assert klass in [list, WeightedTreeCollection]
            klass = WeightedTreeCollection
        else:
            assert klass in [list, LogLikelihoodScoredTreeCollection]
            klass = LogLikelihoodScoredTreeCollection
        tree = make_tree(treestring=line[1])
        trees.append((lnL, tree))
    trees.sort(reverse=True)
    infile.close()
    return klass(trees)

import itertools

from cogent3 import get_moltype, make_tree
from cogent3.evolve.fast_distance import DistanceMatrix, get_calculator
from cogent3.evolve.models import get_model

from .composable import ComposableDistance


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.30a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_calculators = {
    "fast_calculators": ["hamming"],
    "slow_calculators": ["gtr"],
    "fast_and_slow_calculators": ["tn93"],
}


class fast_slow_dist(ComposableDistance):
    """Pairwise distance calculation. Uses fast (but less
    numerically robust) approach where possible, slow (robust)
    approach when not. Returns a DistanceMatrix."""

    def __init__(self, distance=None, moltype="dna", fast_calc=None, slow_calc=None):
        super(fast_slow_dist, self).__init__(input_types="aligned",
                                             output_types=("pairwise_distances", "serialisable"),
                                             data_types=("ArrayAlignment", "Alignment"))
        self._moltype = get_moltype(moltype)
        self._sm = None

        if distance is None and fast_calc is None:
            fast_calc = "hamming"

        if (fast_calc or slow_calc) and distance:
            raise ValueError("cannot combine distance and fast/slow")

        if (
            fast_calc
            and fast_calc.lower()
            not in _calculators["fast_calculators"]
            + _calculators["fast_and_slow_calculators"]
        ):
            raise ValueError("%s is not a fast calculator" % fast_calc)

        if (
            slow_calc
            and slow_calc.lower()
            not in _calculators["slow_calculators"]
            + _calculators["fast_and_slow_calculators"]
        ):
            raise ValueError("%s is not a slow calculator" % slow_calc)

        try:
            fast_calc = get_calculator(distance, moltype)
        except (ValueError, AttributeError):
            self.fast_calc = None

        try:
            self._sm = get_model(distance)
        except ValueError:
            self._sm = None
            self.slow_calc = None

        if type(fast_calc) is str:
            fast_calc = get_calculator(fast_calc, moltype)
        self.fast_calc = fast_calc if fast_calc else get_calculator("hamming", moltype)
        self.slow_calc = slow_calc

    def _est_dist_pair_slow(self, aln):
        """returns distance between seq pairs in aln"""
        assert len(aln.names) == 2
        tree = make_tree(tip_names=aln.names)
        lf = self._sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        lf.set_param_rule("length", is_independent=False)
        lf.optimise(max_restarts=0, show_progress=False)
        dist = 2 * lf.get_param_value("length", edge=aln.names[0])
        return dist

    def __call__(self, aln):
        aln = aln.to_moltype(self._moltype)
        if self.fast_calc:
            self.fast_calc(aln, show_progress=False)
            dists = self.fast_calc.get_pairwise_distances()
        else:
            empty = {p: 0 for p in itertools.product(aln.names, aln.names)}
            dists = DistanceMatrix(empty)
        if self.slow_calc:
            for a in dists.template.names[0]:
                for b in dists.template.names[1]:
                    if not dists[a, b] and a != b:
                        subset = aln.take_seqs([a, b])
                        dist = self._est_dist_pair_slow(subset)
                        dists[a, b] = dists[b, a] = dist
        return dists


def get_fast_slow_calc(distance, **kwargs):
    """returns FastSlow instance for a given distance name"""
    return fast_slow_dist(distance, **kwargs)

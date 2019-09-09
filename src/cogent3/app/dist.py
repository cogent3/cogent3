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


class fast_slow_dist(ComposableDistance):
    """Pairwise distance calculation. Uses fast (but less
    numerically robust) approach where possible, slow (robust)
    approach when not. Returns a DistanceMatrix."""

    def __init__(self, distance=None, moltype=None, fast_calc=None, slow_calc=None):
        super(fast_slow_dist, self).__init__(
            input_types="aligned",
            output_types=("pairwise_distances", "serialisable"),
            data_types=("ArrayAlignment", "Alignment"),
        )
        self._moltype = moltype if moltype is None else get_moltype(moltype)
        self._sm = None

        if (fast_calc or slow_calc) and distance:
            raise ValueError("cannot combine distance and fast/slow")

        if distance:
            fast_calc = distance
            slow_calc = distance

        d = set(["hamming", "paralinear_continuous_time", "logdet"]) & set([slow_calc, fast_calc])
        if d and not self._moltype:
            raise ValueError(f"you must provide a moltype for {d}")

        try:
            fast_calc = get_calculator(fast_calc, moltype=self._moltype)
        except (ValueError, AttributeError):
            fast_calc = None

        try:
            slow_calc = get_model(slow_calc)
        except ValueError:
            slow_calc = None

        if not (fast_calc or slow_calc):
            raise ValueError(f"invalid values for {slow_calc} or {fast_calc}")

        self.fast_calc = fast_calc
        if fast_calc and self._moltype and fast_calc.moltype != self._moltype:
            raise ValueError(
                f"{self._moltype} incompatible moltype with fast calculator {fast_calc.moltype}"
            )
        elif fast_calc:
            self._moltype = fast_calc.moltype

        if slow_calc and self._moltype and slow_calc.moltype != self._moltype:
            raise ValueError("incompatible moltype with slow calculator")
        elif slow_calc:
            self._moltype = slow_calc.moltype
        self._sm = slow_calc

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
        if self._sm:
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

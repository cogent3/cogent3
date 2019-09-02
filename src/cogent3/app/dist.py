from cogent3 import get_moltype, make_tree
from cogent3.evolve.fast_distance import get_calculator
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

    _input_type = frozenset(["aligned"])
    _output_type = frozenset(["pairwise_distances", "serialisable"])
    _data_types = frozenset(["DistanceMatrix"])

    def __init__(self, distance="TN93", moltype="dna", fast_calc=None, slow_calc=None):
        super(fast_slow_dist, self).__init__()
        self._moltype = get_moltype(moltype)
        self._sm = None
        if distance:
            if fast_calc and slow_calc:
                distance = None
                raise UserWarning('distance="%s" being ignored' % distance)
            elif fast_calc or slow_calc:
                raise ValueError("must specify both fast and slow calculator")
            else:
                fast_calc = get_calculator(distance, moltype)
                self._sm = get_model(distance)
                slow_calc = None

        self._distance = distance
        self.slow_calc = slow_calc
        self.fast_calc = fast_calc

    def _est_dist_pair(self, aln):
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
        self.fast_calc(aln, show_progress=False)
        dists = self.fast_calc.get_pairwise_distances()
        for a in dists.template.names[0]:
            for b in dists.template.names[1]:
                if not dists[a, b] and a != b:
                    subset = aln.take_seqs([a, b])
                    dist = self._est_dist_pair(subset)
                    dists[a, b] = dists[b, a] = dist
        return dists


def get_fast_slow_calc(distance, **kwargs):
    """returns FastSlow instance for a given distance name"""
    return fast_slow_dist(distance, **kwargs)

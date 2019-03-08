from cogent3 import LoadTree, get_moltype
from cogent3.evolve.pairwise_distance import get_calculator
from cogent3.evolve.models import get_model

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


# todo tests for dist

class FastSlowDist:  # todo make this a composable type
    def __init__(self, distance='TN93', moltype='dna', fast_calc=None,
                 slow_calc=None):
        self._moltype = get_moltype(moltype)
        self._sm = None
        if distance:
            if fast_calc and slow_calc:
                distance = None
                raise UserWarning('distance="%s" being ignored' % distance)
            elif fast_calc or slow_calc:
                raise ValueError("most specify both fast and slow calculator")
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
        tree = LoadTree(tip_names=aln.names)
        lf = self._sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        lf.set_param_rule('length', is_independent=False)
        lf.optimise(max_restarts=0, show_progress=False)
        dist = 2 * lf.get_param_value('length', edge=aln.names[0])
        return dist

    def __call__(self, aln):
        aln = aln.to_moltype(self._moltype)
        self.fast_calc(aln)
        dists = self.fast_calc.get_pairwise_distances()
        for (a, b) in dists:
            if not dists[(a, b)]:
                subset = aln.take_seqs([a, b])
                dist = self._est_dist_pair(subset)
                dists[(a, b)] = dists[(b, a)] = dist
        return dists


def get_fast_slow_calc(distance, **kwargs):
    """returns FastSlow instance for a given distance name"""
    return FastSlowDist(distance, **kwargs)

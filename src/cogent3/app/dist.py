import itertools

from cogent3 import get_moltype, make_tree
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    get_distance_calculator,
)
from cogent3.evolve.models import get_model

from .composable import (
    ALIGNED_TYPE,
    PAIRWISE_DISTANCE_TYPE,
    SERIALISABLE_TYPE,
    TABULAR_TYPE,
    ComposableDistance,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class fast_slow_dist(ComposableDistance):
    """Pairwise distance calculation for aligned sequences.

    Uses fast (but less numerically robust) approach where possible, slow (robust)
    approach when not. Returns a DistanceMatrix.
    """

    _input_types = ALIGNED_TYPE
    _output_types = (PAIRWISE_DISTANCE_TYPE, TABULAR_TYPE, SERIALISABLE_TYPE)
    _data_types = ("ArrayAlignment", "Alignment")

    def __init__(self, distance=None, moltype=None, fast_calc=None, slow_calc=None):
        """
        Parameters
        ----------
        moltype : str
            cogent3 moltype
        distance : str
            Name of a distance method available as both fast and slow calculator.
        fast_calc
            Name of a fast distance calculator. See cogent3.available_distances().
        slow_calc
            Name of a slow distance calculator. See cogent3.available_models().

        Notes
        -----
        If you provide fast_calc or slow_calc, you must specify the moltype.
        """
        super(fast_slow_dist, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
        self._moltype = moltype if moltype is None else get_moltype(moltype)
        self._sm = None

        if (fast_calc or slow_calc) and distance:
            raise ValueError("cannot combine distance and fast/slow")

        if distance:
            fast_calc = distance
            slow_calc = distance

        d = {"hamming", "percent", "paralinear", "logdet"} & {slow_calc, fast_calc}
        if d and not self._moltype:
            raise ValueError(f"you must provide a moltype for {d}")

        try:
            fast_calc = get_distance_calculator(fast_calc, moltype=self._moltype)
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
        self.func = self.calc_distance

    def _est_dist_pair_slow(self, aln):
        """returns distance between seq pairs in aln"""
        assert len(aln.names) == 2
        tree = make_tree(tip_names=aln.names)
        lf = self._sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        lf.set_param_rule("length", is_independent=False)
        lf.optimise(max_restarts=0, show_progress=False)
        return 2 * lf.get_param_value("length", edge=aln.names[0])

    def calc_distance(self, aln):
        if self._moltype and self._moltype != aln.moltype:
            aln = aln.to_moltype(self._moltype)

        if self.fast_calc:
            self.fast_calc(aln, show_progress=False)
            dists = self.fast_calc.get_pairwise_distances()
        else:
            empty = {p: 0 for p in itertools.product(aln.names, aln.names)}
            dists = DistanceMatrix(empty)
        dists.source = aln.info.source
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

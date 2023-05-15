import itertools

from typing import Union

from numpy import log, polyval, triu_indices, zeros_like

from cogent3 import get_moltype, make_tree
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    get_distance_calculator,
)
from cogent3.evolve.models import get_model
from cogent3.maths.distance_transform import jaccard

from .composable import define_app
from .typing import (
    AlignedSeqsType,
    PairwiseDistanceType,
    SerialisableType,
    UnalignedSeqsType,
)


# The following coefficients are derived from a polynomial fit between Jaccard distance
# and the proportion of different sites for mammalian DNA sequences. NOTE: the Jaccard
# distance used kmers where k=10.
JACCARD_PDIST_POLY_COEFFS = [
    2271.7714153914335,
    -11998.34362001251,
    27525.142573445955,
    -35922.0159776342,
    29337.5102940838,
    -15536.693064681693,
    5346.929667208838,
    -1165.616998965176,
    151.8581396241204,
    -10.489082251524346,
    0.3334853259953467,
    0.0,
]


@define_app
class fast_slow_dist:
    """Pairwise distance calculation for aligned sequences.

    Uses fast (but less numerically robust) approach where possible, slow (robust)
    approach when not.
    """

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

    def _est_dist_pair_slow(self, aln):
        """returns distance between seq pairs in aln"""
        assert len(aln.names) == 2
        tree = make_tree(tip_names=aln.names)
        lf = self._sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        lf.set_param_rule("length", is_independent=False)
        lf.optimise(max_restarts=0, show_progress=False)
        return 2 * lf.get_param_value("length", edge=aln.names[0])

    def main(
        self, aln: AlignedSeqsType
    ) -> Union[SerialisableType, PairwiseDistanceType]:
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


@define_app
def jaccard_dist(seq_coll: UnalignedSeqsType, k: int = 10) -> PairwiseDistanceType:
    """Calculates the pairwise Jaccard distance between the sets of kmers generated from
    sequences in the collection. A measure of distance for unaligned sequences.

    Parameters
    ----------
    seq_coll: UnalignedSeqsType
        a collection of unaligned sequences
    k: int
        size of kmer to use for

    Returns
    -------
    Pairwise Jaccard distance between sequences in the collection.
    """

    kmers = {name: set(seq.get_kmers(k)) for name, seq in seq_coll.named_seqs.items()}
    seq_names = sorted(kmers.keys())
    num_seqs = len(seq_names)

    jaccard_dists = {}

    for i in range(num_seqs):
        for j in range(i):
            name1, name2 = seq_names[i], seq_names[j]
            dist = jaccard(kmers[name1], kmers[name2])
            jaccard_dists[(name1, name2)] = dist
            jaccard_dists[(name2, name1)] = dist

    return DistanceMatrix(jaccard_dists)


@define_app()
def approx_pdist(jaccard_dists: PairwiseDistanceType) -> PairwiseDistanceType:
    """approximate the proportion sites different from Jaccard distances

    Parameters
    ----------
    jaccard_dists
        The pairwise Jaccard distance matrix

    Returns
    -------
    DistanceMatrix of approximated proportion sites different

    Notes
    -----
    Coefficients were derived from a polynomial fit between Jaccard distance
    of kmers with k=10 and the proportion of sites different using mammalian
    106 protein coding gene DNA sequence alignments.
    """
    upper_indices = triu_indices(n=jaccard_dists.shape[0], k=1)
    result = deepcopy(jaccard_dists)  # so the original matrix not modified
    arr = result.array

    # Convert only the upper indices from Jaccard distance to approximate PDist
    upper_vals = polyval(JACCARD_PDIST_POLY_COEFFS, j_dists[upper_indices])
    p_dists[upper_indices] = upper_vals

    # Reflect the upper triangle to the lower triangle
    lower_indices = (upper_indices[1], upper_indices[0])
    p_dists[lower_indices] = upper_vals

    # add dists to dictionary where {(seq, seq) : dist}
    names = jaccard_dists.names
    named_dists = {
        (names[i], names[j]): p_dists[i, j]
        for i, j in itertools.combinations(range(len(names)), 2)
    }

    return DistanceMatrix(named_dists)


@define_app
def approx_jc69(
    pdist_predicted: PairwiseDistanceType,
) -> PairwiseDistanceType:
    """converts p-distances and returns pairwise JC69 distances

    Parameters
    ----------
    pdists
        The pairwise PDist matrix
    num_states
        Number of sequence states, default is for the traditional
        JC69 modelling of DNA substitutions

    Returns
    -------
    DistanceMatrix of pairwise JC69 distances

    """
    pdists = pdist_predicted.array
    jc_dists = zeros_like(pdists)

    # calculate approx jc dist from approx pdist for upper triangle of matrix
    upper_indices = triu_indices(n=pdists.shape[0], k=1)
    upper_vals = _jc69_from_pdist(pdists[upper_indices])
    jc_dists[upper_indices] = upper_vals

    # reflect into lower triangle of matrix
    lower_indices = (upper_indices[1], upper_indices[0])
    jc_dists[lower_indices] = upper_vals

    # add dists to dictionary where {(seq, seq) : dist}, so can be wrapped in DistanceMatrix
    names = pdist_predicted.names
    named_dists = {
        (names[i], names[j]): jc_dists[i, j]
        for i, j in itertools.combinations(range(len(names)), 2)
    }

    return DistanceMatrix(named_dists)


def _jc69_from_pdist(p):
    """convert proportion of sites different to Jukes Cantor distance"""
    return -3.0 * log(1 - (4 / 3) * p) / 4

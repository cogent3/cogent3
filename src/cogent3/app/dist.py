import itertools

from copy import deepcopy
from itertools import combinations
from typing import Callable, Union

import numpy

from numpy import polyval, triu_indices

from cogent3 import get_moltype, make_tree
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    get_distance_calculator,
)
from cogent3.evolve.models import get_model
from cogent3.maths.distance_transform import jaccard
from cogent3.util.misc import get_true_spans

from .composable import NotCompleted, define_app
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

    def __init__(
        self,
        distance: str = None,
        moltype: str = None,
        fast_calc: str = None,
        slow_calc: str = None,
    ):
        """
        Parameters
        ----------
        moltype
            cogent3 moltype
        distance
            Name of a distance method available as both fast and slow calculator.
        fast_calc
            Name of a fast distance calculator. See cogent3.available_distances().
        slow_calc
            Name of a slow distance calculator. See cogent3.available_models().

        Notes
        -----
        If you provide fast_calc or slow_calc, you must specify the moltype.

        Examples
        --------

        Create a sample alignment and app to calculate the pairwise hamming
        distance.

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> seqs = {
        ...     "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        ...     "Bandicoot": "---NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
        ...     "Rhesus": "GCCAGCTCATTACAGCATGAGAAC---AGTTTGTTACTCACT",
        ...     "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT"
        ... }
        >>> aln = make_aligned_seqs(seqs, moltype="dna")
        >>> app = get_app("fast_slow_dist", distance="hamming", moltype="dna")
        >>> result = app(aln)
        >>> result.to_dict()
        {('Bandicoot', 'FlyingFox'): 11.0, ('Bandicoot', 'Human'): 11.0, ('Bandicoot', 'Rhesus'): 12.0,...

        Create an app with ``fast_calc="TN93"``.

        >>> app = get_app("fast_slow_dist", fast_calc="TN93", moltype="dna")
        >>> result = app(aln)
        >>> result.to_dict()
        {('Bandicoot', 'FlyingFox'): 0.4494289084991177,...

        """
        self._moltype = moltype if moltype is None else get_moltype(moltype)
        self._sm = None

        if (fast_calc or slow_calc) and distance:
            raise ValueError("cannot combine distance and fast/slow")

        if distance:
            fast_calc = distance
            slow_calc = distance

        d = {"hamming", "pdist", "paralinear", "logdet"} & {slow_calc, fast_calc}
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
        size of kmer. Default is 10.

    Returns
    -------
    PairwiseDistanceType
        Pairwise Jaccard distance between sequences in the collection.

    Examples
    --------

    Create an sample sequence collection and an app that calculates the
    pairwise Jaccard distance with kmer size ``k=2``.

    >>> from cogent3 import make_unaligned_seqs, get_app
    >>> seqs = make_unaligned_seqs({"s1": "ACGTA", "s2": "ACGTC"}, moltype="dna")
    >>> app = get_app("jaccard_dist", k=2)
    >>> result = app(seqs)
    >>> print(result.to_dict())
    {('s1', 's2'): 0.4, ('s2', 's1'): 0.4}
    """

    kmers = {
        name: set(seq.iter_kmers(k, strict=True))
        for name, seq in seq_coll.named_seqs.items()
    }
    seq_names = sorted(kmers.keys())
    num_seqs = len(seq_names)

    dists = {}

    for i in range(num_seqs):
        for j in range(i):
            name1, name2 = seq_names[i], seq_names[j]
            dist = jaccard(kmers[name1], kmers[name2])
            dists[(name1, name2)] = dist
            dists[(name2, name1)] = dist
    if not dists:
        names = ", ".join(f"{n!r}" for n in seq_coll.names)
        return NotCompleted(
            "ERROR",
            jaccard_dist.__name__,
            f"could not compute distances between {names}",
            source=seq_coll,
        )
    return DistanceMatrix(dists)


@define_app
def approx_pdist(jaccard_dists: PairwiseDistanceType) -> PairwiseDistanceType:
    """Calculates an approximation of the p-distance between sequences based
    on Jaccard distances (see Notes for details).

    Parameters
    ----------
    jaccard_dists : DistanceMatrix
        A DistanceMatrix containing pairwise Jaccard distances between sequences.

    Returns
    -------
    DistanceMatrix
        Pairwise approximated p-distance between sequences.

    Notes
    -----
    This approximation assumes the input DistanceMatrix contains Jaccard
    distances. See the ``jaccard_dist`` app for calculating Jaccard distances.

    The true p-distance measures the proportion of differing sites between two
    aligned sequences. For the true p-distance, 0 indicates identical sequences,
    and 1 indicates completely different sequences.

    This approximation does not guarantee bounds of 0 and 1.

    The method approximates the p-distance using coefficients derived from a
    polynomial fit between Jaccard distance and p-distance. The coefficients
    were fitted using data from 106 DNA sequences of mammalian protein coding
    genes, with kmers of size k=10.

    Examples
    --------

    Create a sample sequence collection and apps for calculating the Jaccard
    distances, and an approximation of the p-distance.

    >>> from cogent3 import make_unaligned_seqs, get_app
    >>> seqs = make_unaligned_seqs({"s1": "ACGTA", "s2": "ACGTC"}, moltype="dna")
    >>> jaccard_dist = get_app("jaccard_dist", k=2)
    >>> pdist = get_app("approx_pdist")

    Create a composable app to return the p-distances.

    >>> app = jaccard_dist + pdist
    >>> result = app(seqs)
    >>> print(result.to_dict())
    {('s1', 's2'): 0.02597466474649073, ('s2', 's1'): 0.02597466474649073}
    """
    upper_indices = triu_indices(n=jaccard_dists.shape[0], k=1)
    result = deepcopy(jaccard_dists)  # so the original matrix not modified
    arr = result.array

    # Convert only the upper indices from Jaccard distance to approximate PDist
    arr[upper_indices] = polyval(JACCARD_PDIST_POLY_COEFFS, arr[upper_indices])
    # Reflect the upper triangle to the lower triangle
    arr.T[upper_indices] = arr[upper_indices]
    result.array = arr
    return result


@define_app
def approx_jc69(
    pdists: PairwiseDistanceType, num_states: int = 4
) -> PairwiseDistanceType:
    """Converts p-distances and returns pairwise JC69 distances

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

    Examples
    --------

    Create a sample sequence collection and apps for calculating the Jaccard
    distances, and an approximation of the pairwise JC69 distances.

    >>> from cogent3 import make_unaligned_seqs, get_app
    >>> seqs = make_unaligned_seqs({"s1": "ACGTA", "s2": "ACGTC"}, moltype="dna")
    >>> jaccard_dist = get_app("jaccard_dist", k=2)
    >>> jc69 = get_app("approx_jc69")

    Create a composable app to return the JC69 distances.

    >>> app = jaccard_dist + jc69
    >>> result = app(seqs)
    >>> print(result.to_dict())
    {('s1', 's2'): 0.5716050390351726, ('s2', 's1'): 0.5716050390351726}
    """
    upper_indices = triu_indices(n=pdists.shape[0], k=1)
    result = deepcopy(pdists)  # so the original matrix not modified
    arr = result.array
    n_1 = num_states - 1
    arr[upper_indices] = (
        -n_1 / num_states * numpy.log(1 - (num_states / n_1 * arr[upper_indices]))
    )
    arr.T[upper_indices] = arr[upper_indices]
    result.array = arr
    return result


def get_approx_dist_calc(dist: str, num_states: int = 4) -> Callable:
    """Return the corresponding callable approximate distance calculator

    Parameters
    ----------
    dist : str
        The distance measure for the calculator, either "pdist" or "jc69"
    num_states
        Number of sequence states, default is for the traditional
        JC69 modelling of DNA substitutions

    Returns
    -------
    cogent3.app for calculating DistanceMatrix for a SequenceCollection
    """
    if dist == "pdist":
        dist_calc_app = jaccard_dist() + approx_pdist()
    elif dist == "jc69":
        dist_calc_app = (
            jaccard_dist() + approx_pdist() + approx_jc69(num_states=num_states)
        )
    else:
        raise ValueError(f"No support for calc={dist}. Use either 'pdist' or 'jc69'")
    return dist_calc_app


@define_app
class gap_dist:
    """compute the pairwise difference in gaps using affine gap score"""

    def __init__(self, gap_insert: float = 12.0, gap_extend: float = 1.0):
        """
        Parameters
        ----------
        gap_insert
            gap insertion penalty
        gap_extend
            gap extension penalty

        Examples
        --------

        Create a sample alignment and app to calculate the pairwise difference in gaps.

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> seqs = {
        ...     "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        ...     "Bandicoot": "---NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
        ...     "Rhesus": "GCCAGCTCATTACAGCATGAGAAC---AGTTTGTTACTCACT",
        ...     "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
        ... }
        >>> aln = make_aligned_seqs(seqs, moltype="dna")
        >>> app = get_app("gap_dist")
        >>> result = app(aln)
        >>> result.to_dict()
        {('Bandicoot', 'FlyingFox'): 30.0, ('Bandicoot', 'Human'): 15.0, ('Bandicoot', 'Rhesus'): 30.0,...

        Calculate the pairwise difference in gaps with different ``gap_insert``
        and ``gap_extend`` penalties.

        >>> app = get_app("gap_dist", gap_insert=10, gap_extend=2.8)
        >>> result = app(aln)
        >>> result.to_dict()
        {('Bandicoot', 'FlyingFox'): 36.8, ('Bandicoot', 'Human'): 18.4, ('Bandicoot', 'Rhesus'): 36.8,...

        """
        self._insert = gap_insert
        self._extend = gap_extend

    def main(self, aln: AlignedSeqsType) -> PairwiseDistanceType:
        gap_diffs = {}
        gap_data = dict(zip(aln.names, aln.get_gap_array()))
        # convert each gap vector to gap spans
        for k, gap_vector in gap_data.items():
            gap_data[k] = set(tuple(pr) for pr in get_true_spans(gap_vector).tolist())

        for i, j in combinations(range(aln.num_seqs), 2):
            n1, n2 = aln.names[i], aln.names[j]

            # unique gaps
            unique = gap_data[n1] ^ gap_data[n2]

            # compute the score as
            # number of gaps * gap_insert + total length of gaps * gap_extend
            score = self._insert * len(unique) + self._extend * sum(
                l for p, l in unique
            )
            gap_diffs[(n1, n2)] = gap_diffs[(n2, n1)] = score

        return DistanceMatrix(gap_diffs)

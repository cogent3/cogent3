import warnings

from bisect import bisect_left
from itertools import combinations
from typing import Union

from numpy import array

from cogent3 import make_tree
from cogent3.align import (
    global_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.align.progressive import tree_align
from cogent3.app import dist
from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.location import gap_coords_to_map
from cogent3.core.moltype import get_moltype
from cogent3.evolve.fast_distance import get_distance_calculator
from cogent3.evolve.models import get_model
from cogent3.maths.util import safe_log

from .composable import NotCompleted, define_app
from .tree import quick_tree, scale_branches
from .typing import AlignedSeqsType, SerialisableType, UnalignedSeqsType


class _GapOffset:
    """computes sum of gap lengths preceding a position. Acts like a dict
    for getting the offset for an integer key with the __getitem__ returning
    the offset.
    If your coordinate is an alignment position, set invert=True.
    Examples
    --------
    From sequence coordinate to an alignment coordinate

    >>> seq2aln = _GapOffset({1:3, 7:1})
    >>> seq_pos = 2
    >>> aln_pos = seq_pos + seq2aln[seq_pos]
    >>> aln_pos
    5

    From alignment coordinate to a sequence coordinate

    >>> aln2seq = _GapOffset({1:3, 7:1}, invert=True)
    >>> seq_pos = aln_pos - aln2seq[seq_pos]
    >>> seq_pos
    2
    """

    def __init__(self, gaps_lengths, invert=False):
        """
        Parameters
        ----------
        gaps_lengths : dict
            {pos: length, ...} where pos is a gap insert position and length
            how long the gap is.
        invert : bool
            if True, query keys are taken as being in alignment coordinates
        """
        offset = 0
        min_val = None
        result = {}
        k = -1
        for k, l in sorted(gaps_lengths.items()):
            if invert:
                result[k + offset + l] = offset + l
                result[k + offset] = offset
            else:
                result[k] = offset

            offset += l
            if min_val is None:
                min_val = k

        self._store = result
        self.min_pos = min_val
        self.max_pos = k + offset if invert else k
        self.total = offset
        self._ordered = None
        self._invert = invert

    def __repr__(self):
        return repr(self._store)

    def __str__(self):
        return str(self._store)

    def __getitem__(self, k):
        if not self._store:
            return 0

        if k in self._store:
            return self._store[k]

        if k < self.min_pos:
            return 0

        if k > self.max_pos:
            return self.total

        if self._ordered is None:
            self._ordered = sorted(self._store)

        # k is definitely bounded by min and max positions here
        i = bisect_left(self._ordered, k)
        pos = self._ordered[i]
        if self._invert:
            pos = pos if pos in [k, 0] else self._ordered[i - 1]
        return self._store[pos]


def _gap_union(seqs) -> dict:
    """returns the union of all gaps in seqs"""
    seq_name = None
    all_gaps = {}
    for seq in seqs:
        if not isinstance(seq, Aligned):
            raise TypeError(f"must be Aligned instances, not {type(seq)}")
        if seq_name is None:
            seq_name = seq.name
        if seq.name != seq_name:
            raise ValueError("all sequences must have the same name")

        gaps_lengths = dict(seq.map.get_gap_coordinates())
        all_gaps = _merged_gaps(all_gaps, gaps_lengths)
    return all_gaps


def _gap_difference(seq_gaps: dict, union_gaps: dict) -> tuple:
    """

    Parameters
    ----------
    seq_gaps
        {gap pos: length, } of a sequence used to generate the gap union
    union_gaps
        {gap pos: maximum length, ...} derived from the same seq aligned
        to different sequences

    Returns
    -------
    gaps missing from seq_gaps, seq_gaps that overlap with union gaps
    """
    missing = {}
    overlapping = {}
    for position, length in union_gaps.items():
        if position not in seq_gaps:
            missing[position] = length
        elif seq_gaps[position] != length:
            overlapping[position] = length - seq_gaps[position]
    return missing, overlapping


def _merged_gaps(a_gaps: dict, b_gaps: dict) -> dict:
    """merges gaps that occupy same position

    Parameters
    ----------
    a_gaps, b_gaps
        [(gap position, length),...]

    Returns
    -------
    Merged places as {gap position: length, ...}

    Notes
    -----
    If a_gaps and b_gaps are from the same underlying sequence, set
    function to 'max'. Use 'sum' when the gaps derive from different
    sequences.
    """

    if not a_gaps:
        return b_gaps

    if not b_gaps:
        return a_gaps

    places = set(a_gaps) | set(b_gaps)
    return {
        place: max(
            a_gaps.get(place, 0),
            b_gaps.get(place, 0),
        )
        for place in places
    }


def _subset_gaps_to_align_coords(
    subset_gaps: dict, orig_gaps: dict, seq_2_aln: _GapOffset
) -> dict:
    """compute alignment coords of subset gaps

    Parameters
    ----------
    subset_gaps : dict
        {position: length delta} lengths are the adjusted gap lengths
    orig_gaps : dict
        {position: orig length} the original gap lengths are from a pairwise
        alignment
    seq_2_aln : dict
        {seq position: alignment position, ...}

    Returns
    -------
    dict
        {alignment position + orig length: length delta, ...}

    Notes
    -----
    """
    result = {}
    for p in subset_gaps:
        offset = seq_2_aln[p]
        result[offset + p + orig_gaps[p]] = subset_gaps[p]

    return result


def _combined_refseq_gaps(seq_gaps: dict, union_gaps: dict) -> dict:
    # takes union gaps and refseq gaps, converts into diffs and
    # subset diffs
    seq2aln = _GapOffset(seq_gaps)
    diff_gaps, subset_gaps = _gap_difference(seq_gaps, union_gaps)
    align_coord_gaps = _subset_gaps_to_align_coords(subset_gaps, seq_gaps, seq2aln)
    align_coord_gaps.update({p + seq2aln[p]: diff_gaps[p] for p in diff_gaps})
    return align_coord_gaps


def _gaps_for_injection(other_seq_gaps: dict, refseq_gaps: dict, seqlen: int) -> dict:
    """projects refseq aligned gaps into otherseq

    Parameters
    ----------
    other_seq_gaps : dict
        {gap in other seq position: gap length}
    refseq_gaps : dict
        {gap as alignment position: gap length}
    seqlen : int
        length of sequence being injected into

    Returns
    -------
    dict
        {gap in other seq position: gap length}
    """
    aln2seq = _GapOffset(other_seq_gaps, invert=True)
    # to inject a gap means to convert it from alignment coordinates into
    # sequence coordinates
    # we probably need to include the refseq gap union because we need to
    # establish whether a refseq gap overlaps with a gap in other seq
    # and
    all_gaps = {}
    all_gaps.update(other_seq_gaps)
    for gap_pos, gap_length in sorted(refseq_gaps.items()):
        offset = aln2seq[gap_pos]
        gap_pos -= offset
        gap_pos = min(seqlen, gap_pos)
        if gap_pos < 0:
            raise ValueError(
                f"computed gap_pos {gap_pos} < 0, correct reference sequence?"
            )
        if gap_pos in all_gaps:
            gap_length += all_gaps[gap_pos]

        all_gaps[gap_pos] = gap_length

    return all_gaps


def pairwise_to_multiple(pwise, ref_seq, moltype, info=None):
    """
    turns pairwise alignments to a reference into a multiple alignment

    Parameters
    ----------
    pwise
        Series of pairwise alignments to ref_seq as
        [(non-refseq name, aligned pair), ...]
    ref_seq
        The sequence common in all pairwise alignments
    moltype
        molecular type for the returned alignment
    info
        info object

    Returns
    -------
    ArrayAlign
    """
    if not hasattr(ref_seq, "name"):
        raise TypeError(f"ref_seq must be a cogent3 sequence, not {type(ref_seq)}")

    refseqs = [s for _, aln in pwise for s in aln.seqs if s.name == ref_seq.name]
    ref_gaps = _gap_union(refseqs)

    m = gap_coords_to_map(ref_gaps, len(ref_seq))
    aligned = [Aligned(m, ref_seq)]
    for other_name, aln in pwise:
        curr_ref = aln.named_seqs[ref_seq.name]
        curr_ref_gaps = dict(curr_ref.map.get_gap_coordinates())
        other_seq = aln.named_seqs[other_name]
        other_gaps = dict(other_seq.map.get_gap_coordinates())
        diff_gaps = _combined_refseq_gaps(curr_ref_gaps, ref_gaps)
        inject = _gaps_for_injection(other_gaps, diff_gaps, len(other_seq.data))
        if inject:
            m = gap_coords_to_map(inject, len(other_seq.data))
            other_seq = Aligned(m, other_seq.data)

        aligned.append(other_seq)
    # default to ArrayAlign
    return Alignment(aligned, moltype=moltype, info=info).to_type(
        array_align=True, moltype=moltype
    )


@define_app
class align_to_ref:
    """Aligns sequences to a nominated reference in the unaligned collection.

    This is much faster, and requires much less memory, than progressive_align
    but the quality will likely be lower. Alignment quality will be strongly
    affected by choice of reference.
    """

    def __init__(
        self,
        ref_seq="longest",
        score_matrix=None,
        insertion_penalty=20,
        extension_penalty=2,
        moltype="dna",
    ):
        """
        Parameters
        ----------
        ref_seq : str
            either a name to be found in the data, or 'longest'.
            If latter, the longest sequence will be chosen as the reference
        score_matrix
            scoring dict for DNA, defaults to `make_dna_scoring_dict(10, -1, -8)`
        insertion_penalty
            penalty for gap insertion
        extension_penalty
            penalty for gap extension
        moltype : str
            molecular type
        """
        assert moltype
        moltype = get_moltype(moltype)
        self._moltype = moltype
        S = score_matrix or (
            make_dna_scoring_dict(10, -1, -8)
            if self._moltype.label == "dna"
            else make_generic_scoring_dict(10, self._moltype)
        )
        self._kwargs = dict(
            S=S, d=insertion_penalty, e=extension_penalty, return_score=False
        )
        if ref_seq.lower() == "longest":
            self._func = self.align_to_longest
        else:
            self._func = self.align_to_named_seq
            self._ref_name = ref_seq

        self._gap_state = None  # can be character or int, depends on aligner

    def align_to_longest(self, seqs):
        """returns alignment to longest seq"""
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        lengths = seqs.get_lengths()
        lengths = [(l, n) for n, l in lengths.items()]
        _, ref_seq_name = max(lengths)
        self._ref_name = ref_seq_name
        return self.align_to_named_seq(seqs)

    def align_to_named_seq(self, seqs):
        """returns alignment to named seq"""
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        ref_seq = seqs.get_seq(self._ref_name)
        kwargs = self._kwargs.copy()
        pwise = []
        for seq in seqs.seqs:
            if seq.name == self._ref_name:
                continue

            aln = global_pairwise(ref_seq, seq, **kwargs).to_type(array_align=False)
            pwise.append((seq.name, aln))

        return pairwise_to_multiple(pwise, ref_seq, self._moltype, info=seqs.info)

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, seqs: UnalignedSeqsType) -> T:
        """return aligned sequences"""
        return self._func(seqs)


@define_app
class progressive_align:
    """Progressive multiple sequence alignment via any cogent3 model."""

    def __init__(
        self,
        model,
        gc=None,
        param_vals=None,
        guide_tree=None,
        unique_guides=False,
        indel_length=1e-1,
        indel_rate=1e-10,
        distance="pdist",
    ):
        """
        Parameters
        ----------
        model
            substitution model instance or name. If 'codon'
            (uses MG94HKY), 'nucleotide' (uses HKY85), 'protein'
            (uses WG01). These choices provide also provide default
            settings for param_vals.
        gc : int or string
            the genetic code for a codon alignment, defaults to the standard
            genetic code
        param_vals : dict
            param name, values for parameters in model. Overrides
            default choices.
        guide_tree
            newick string, tree instance (must have branch lengths), or a
            callable that will build a tree from unaligned collection. If not
            provided, estimated ONCE via constructing a crude alignment. In the
            case of callable, or not provided, the computed guide tree is stored
            in the returned alignment.info['guide_tree'].
        unique_guides : bool
            whether each alignment requires a new guide tree
        indel_rate : float
            probability of gap insertion
        indel_length : float
            probability of gap extension
        distance : string
            the distance measure for building a guide tree. Default is 'percent',
            the proportion of differences. This is applicable for any moltype,
            and sequences with very high percent identity. For more diverged
            sequences we recommend 'paralinear'.
        """
        self._param_vals = {
            "codon": dict(omega=0.4, kappa=3),
            "nucleotide": dict(kappa=3),
        }.get(model, param_vals)
        sm = {"codon": "MG94HKY", "nucleotide": "HKY85", "protein": "JTT92"}.get(
            model, model
        )
        kwargs = {} if gc is None else dict(gc=gc)
        sm = get_model(sm, **kwargs)
        moltype = sm.alphabet.moltype
        self._model = sm
        self._scalar = sm.word_length
        self._indel_length = indel_length
        self._indel_rate = indel_rate
        self._moltype = moltype
        self._unique_guides = unique_guides
        self._distance = distance
        if callable(guide_tree):
            self._make_tree = guide_tree
            guide_tree = None  # callback takes precedence
        else:
            al_to_ref = align_to_ref(moltype=self._moltype)
            dist_calc = dist.fast_slow_dist(
                distance=self._distance, moltype=self._moltype
            )
            est_tree = quick_tree()
            self._make_tree = al_to_ref + dist_calc + est_tree

        if guide_tree is not None:
            if type(guide_tree) == str:
                guide_tree = make_tree(treestring=guide_tree, underscore_unmunge=True)
                if guide_tree.children[0].length is None:
                    raise ValueError("Guide tree must have branch lengths")
            # make sure no zero lengths
            guide_tree = scale_branches()(guide_tree)

        self._guide_tree = guide_tree
        self._kwargs = dict(
            indel_length=self._indel_length,
            indel_rate=self._indel_rate,
            tree=self._guide_tree,
            param_vals=self._param_vals,
            show_progress=False,
        )

    def _build_guide(self, seqs):
        tree = self._make_tree(seqs)
        if self._scalar != 1:
            scaler = scale_branches(scalar=self._scalar)
            tree = scaler(tree)
        return tree

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, seqs: UnalignedSeqsType) -> T:
        """returned progressively aligned sequences"""
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        if self._guide_tree is None or self._unique_guides:
            self._guide_tree = self._build_guide(seqs)
            self._kwargs["tree"] = self._guide_tree
            diff = set(self._guide_tree.get_tip_names()) ^ set(seqs.names)
            if diff:
                numtips = len(set(self._guide_tree.get_tip_names()))
                print(f"numseqs={len(seqs.names)} not equal to numtips={numtips}")
                print(f"These were different: {diff}")
                seqs = seqs.take_seqs(self._guide_tree.get_tip_names())

        kwargs = self._kwargs.copy()

        with warnings.catch_warnings(record=False):
            warnings.simplefilter("ignore")
            try:
                result, tree = tree_align(self._model, seqs, **kwargs)
                if self._moltype and self._moltype != result.moltype:
                    result = result.to_moltype(self._moltype)

                result.info.update(seqs.info)
            except ValueError as err:
                # probably an internal stop
                result = NotCompleted("ERROR", self, err.args[0], source=seqs)
                return result
        return result


@define_app
class ic_score:
    """compute the Information Content alignment quality score

    Returns
    -------
    The score or 0.0 if it cannot be computed.

    Notes
    -----
    Based on eq. (2) in noted reference.

    Hertz, G. Z. & Stormo, G. D. Identifying DNA and protein patterns with
    statistically significant alignments of multiple sequences.
    Bioinformatics vol. 15 563–577 (Oxford University Press, 1999)
    """

    def __init__(self, equifreq_mprobs=True):
        """
        Parameters
        ----------
        equifreq_mprobs : bool
            If true, specifies equally frequent motif probabilities.
        """
        self._equi_frequent = equifreq_mprobs

    def main(self, aln: AlignedSeqsType) -> float:
        counts = aln.counts_per_pos(include_ambiguity=False, allow_gap=False)
        if counts.array.max() == 0 or aln.num_seqs == 1:
            return 0.0

        motif_probs = aln.get_motif_probs(include_ambiguity=False, allow_gap=False)

        if self._equi_frequent:
            # we reduce motif_probs to observed states
            motif_probs = {m: v for m, v in motif_probs.items() if v > 0}
            num_motifs = len(motif_probs)
            motif_probs = {m: 1 / num_motifs for m in motif_probs}

        p = array([motif_probs.get(b, 0.0) for b in counts.motifs])

        cols = p != 0
        p = p[cols]
        counts = counts.array[:, cols]
        frequency = counts / aln.num_seqs
        log_f = safe_log(frequency / p)
        I_seq = log_f * frequency
        return I_seq.sum()


@define_app
def cogent3_score(aln: AlignedSeqsType) -> float:
    """returns the cogent3 log-likelihood as an alignment quality score

    Parameters
    ----------
    aln
        An alignment instance.

    Returns
    -------
    returns the log-likelihood, or 0.0 if the alignment does not have the
    score

    Notes
    -----
    The cogent3 tree_align algorithm is a progressive-HMM. It records
    the log-likelihood of the alignment (in addition to all other
    the parameter values, including the guide tree) in
    ``alignment.info['align_params']``. This can be used as a measure of
    alignment quality.

    The instance must have been aligned using cogent3 tree_align. In addition,
    if the alignment has been saved, it has must have been serialised
    using a format that preserves the score.
    """
    if aln.num_seqs == 1 or len(aln) == 0:
        return 0.0

    align_params = aln.info.get("align_params", {})
    return align_params.get("lnL", 0.0)


@define_app
class sp_score:
    """Compute a variant of Sum of Pairs alignment quality score

    Returns
    -------
    The score or 0.0 if it cannot be computed.

    Notes
    -----
    We first compute pairwise genetic distances by the user specified
    method. For each pair, this genetic distance is converted into an
    estimate of the number of unchanged positions (with the minimum being
    zero). The result is a sequence similarity matrix.

    We compute pairwise affine gap distance using the user specified
    values. The first matrix minus the second gives the statistic for
    individual pairs and their sum is the Sum of Pairs score.

    Carrillo, H. & Lipman, D. The Multiple Sequence Alignment Problem in
    Biology. Siam J Appl Math 48, 1073–1082 (1988).
    """

    def __init__(
        self, calc: str = "JC69", gap_insert: float = 12.0, gap_extend: float = 1.0
    ):
        """
        Parameters
        ----------
        calc
            name of the pairwise distance calculator
        gap_insert
            gap insertion penalty
        gap_extend
            gap extension penalty

        Notes
        -----
        see available_distances() for the list of available methods.
        """
        self._insert = gap_insert
        self._extend = gap_extend
        self._calc = get_distance_calculator(calc)
        self._gap_calc = dist.gap_dist(gap_insert=gap_insert, gap_extend=gap_extend)

    def main(self, aln: AlignedSeqsType) -> float:
        # this is a distance matrix, we want a similarity matrix
        # because genetic distances can exceed 1, we need to adjust by
        # multiplying by the length of the alignment to get the estimated
        # number of changes and subtract that from the alignment length
        if aln.num_seqs == 1 or len(aln) == 0:
            return 0.0
        self._calc(aln, show_progress=False)
        dmat = self._calc.dists
        lengths = self._calc.lengths
        for i, j in combinations(range(aln.num_seqs), 2):
            n1, n2 = dmat.names[i], dmat.names[j]
            length = lengths[n1, n2]
            unchanged = length - dmat[i, j] * length
            # for very divergent cases, unchanged could be < 0,
            # limit this to 0
            dmat[i, j] = dmat[j, i] = max(unchanged, 0.0)
        # compute the pairwise gap as number of gaps * gap_insert + length
        # of gaps * gap_extend
        gap_scores = self._gap_calc(aln)
        # make sure the two dmats in same name order
        gap_scores = gap_scores.take_dists(dmat.names)
        dmat.array -= gap_scores
        return dmat.array.sum() / 2

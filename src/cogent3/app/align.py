import warnings
from bisect import bisect_left
from itertools import combinations
from typing import Optional, Union

from numpy import array, isnan

from cogent3.align import (
    classic_align_pairwise,
    global_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.align.progressive import tree_align
from cogent3.app import dist
from cogent3.app.tree import interpret_tree_arg
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

    >>> seq2aln = _GapOffset({1: 3, 7: 1})
    >>> seq_pos = 2
    >>> aln_pos = seq_pos + seq2aln[seq_pos]
    >>> aln_pos
    5

    From alignment coordinate to a sequence coordinate

    >>> aln2seq = _GapOffset({1: 3, 7: 1}, invert=True)
    >>> seq_pos = aln_pos - aln2seq[aln_pos]
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
        cum_gap_length = 0
        result = {}
        gap_pos = -1
        for gap_pos, gap_length in sorted(gaps_lengths.items()):
            if invert:
                result[gap_pos + cum_gap_length] = cum_gap_length
                result[gap_pos + cum_gap_length + gap_length] = (
                    cum_gap_length + gap_length
                )
            else:
                result[gap_pos] = cum_gap_length

            cum_gap_length += gap_length

        self._store = result
        self.min_pos = min(gaps_lengths) if gaps_lengths else None
        self.max_pos = gap_pos + cum_gap_length if invert else gap_pos
        self.total = cum_gap_length
        self._ordered = None
        self._invert = invert

    def __repr__(self):
        return repr(self._store)

    def __str__(self):
        return str(self._store)

    def __getitem__(self, index):
        if not self._store:
            return 0

        if index in self._store:
            return self._store[index]

        if index < self.min_pos:
            return 0

        if index > self.max_pos:
            return self.total

        if self._ordered is None:
            self._ordered = sorted(self._store)

        # index is definitely bounded by min and max positions here
        i = bisect_left(self._ordered, index)
        pos = self._ordered[i]
        if self._invert:
            pos = pos if pos in [index, 0] else self._ordered[i - 1]
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
    # todo convert to using IndelMap functions
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
    # todo convert these functions to using IndelMap and the numpy set
    #  operation functions
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
        if inject := _gaps_for_injection(other_gaps, diff_gaps, len(other_seq.data)):
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
        iters: Optional[int] = None,
        approx_dists: bool = True,
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
            the distance measure for building a guide tree. Default is 'pdist',
            the proportion of differences. This is applicable for any moltype,
            and sequences with very high percent identity. For more diverged
            sequences we recommend 'paralinear'.
        iters
            the number of times the alignment process is repeated. The guide tree
            is updated on each iteration from pairwise distances computed from the
            alignment produced by the previous iteration. If None, does not do any
            iterations.
        approx_dists
            if no guide tree, and model is for DNA / Codons, estimates pairwise
            distances using an approximation and JC69. Otherwise, estimates
            genetic distances from pairwise alignments (which is slower).

        Examples
        --------

        Create an unaligned sequence collection of BRCA1 genes from 4 species,
        and an app for alignment with nucleotide model ``model="HKY85"``.

        >>> from cogent3 import make_unaligned_seqs, get_app
        >>> aln = make_unaligned_seqs(
        ...     {
        ...         "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        ...         "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
        ...         "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
        ...         "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
        ...     },
        ...     moltype="dna",
        ... )
        >>> app = get_app("progressive_align", model="HKY85")
        >>> result = app(aln)
        >>> print(
        ...     result.to_pretty(
        ...         name_order=["Human", "Bandicoot", "Rhesus", "FlyingFox"]
        ...     )
        ... )
            Human    GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT...

        Optionally, a pre-computed guide tree can be provided.

        >>> newick = "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus:0.06," "Human:0.0):0.04);"
        >>> app_guided = get_app("progressive_align", model="HKY85", guide_tree=newick)
        >>> result = app_guided(aln)
        >>> print(
        ...     result.to_pretty(
        ...         name_order=["Human", "Bandicoot", "Rhesus", "FlyingFox"]
        ...     )
        ... )
            Human    GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT
        Bandicoot    NA.TCA.T.A.G.TTG.AACC.G...---......GTC..AC
           Rhesus    ..........................---...G.........
        FlyingFox    ........T.................---.......TA....
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
        moltype = sm.moltype
        self._model = sm
        self._scalar = sm.word_length
        self._indel_length = indel_length
        self._indel_rate = indel_rate
        self._moltype = moltype
        self._unique_guides = unique_guides
        self._distance = distance
        self._iters = iters
        if callable(guide_tree):
            self._make_tree = guide_tree
            guide_tree = None  # callback takes precedence
        elif approx_dists and len(moltype.alphabet) == 4:
            dist_app = dist.get_approx_dist_calc(dist="jc69", num_states=4)
            est_tree = quick_tree()
            self._make_tree = dist_app + est_tree
        else:
            al_to_ref = align_to_ref(moltype=self._moltype)
            dist_calc = dist.fast_slow_dist(
                distance=self._distance, moltype=self._moltype
            )
            est_tree = quick_tree()
            self._make_tree = al_to_ref + dist_calc + est_tree

        if guide_tree is not None:
            guide_tree = interpret_tree_arg(guide_tree)
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
            iters=self._iters,
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
            if not self._guide_tree:
                return self._guide_tree
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
class smith_waterman:
    """Local alignment of two sequences using the Smith-Waterman algorithm

    Notes
    -----
    It records the score of the alignment (in addition to all other the
    parameter values) as "sw_score" in ``alignment.info['align_params']``.
    """

    def __init__(
        self,
        score_matrix: dict = None,
        insertion_penalty: int = 20,
        extension_penalty: int = 2,
        moltype: str = "dna",
    ):
        """
        Parameters
        ----------
        score_matrix
            scoring dict, defaults to `make_dna_scoring_dict(10, -1, -8)` for
            DNA and `make_generic_scoring_dict(10, moltype)` for other moltype.
        insertion_penalty
            penalty for gap insertion
        extension_penalty
            penalty for gap extension
        moltype
            molecular type of sequences, defaults to "dna"

        Note
        ----
        If the provided molecular type differs from the moltype of the
        SequenceCollection to be aligned, the sequences are converted to
        the provided moltype.

        Examples
        --------

        Create and align two sequences using the Smith-Waterman algorithm.

        >>> from cogent3 import make_unaligned_seqs, get_app
        >>> seqs = make_unaligned_seqs(
        ...     {
        ...         "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        ...         "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
        ...     },
        ...     moltype="dna",
        ... )
        >>> app = get_app("smith_waterman")
        >>> result = app(seqs)
        >>> print(result.to_pretty())
            Human    AGCTCATTACAGCATGAGAACAGCAGTTTATTACTCA
        Bandicoot    NA.......AT..T...A.C............GTC..

        Align with a less stringent scoring matrix. i.e. for more divergent
        sequences.

        >>> from cogent3.align import make_dna_scoring_dict
        >>> app = get_app(
        ...     "smith_waterman", score_matrix=make_dna_scoring_dict(5, -2, -5)
        ... )
        >>> result = app(seqs)
        >>> print(result.to_pretty())
            Human    AGCTCATTACAGCATGAGAACAGCAGTTTATTACTCA
        Bandicoot    NA.......AT..T...A.C............GTC..
        """
        self.moltype = get_moltype(moltype)
        self._score_matrix = score_matrix or (
            make_dna_scoring_dict(10, -1, -8)
            if self.moltype.label == "dna"
            else make_generic_scoring_dict(10, self.moltype)
        )
        self._insertion_penalty = insertion_penalty
        self._extension_penalty = extension_penalty

    def main(self, seqs: UnalignedSeqsType) -> AlignedSeqsType:
        if seqs.num_seqs > 2:
            return NotCompleted(
                "ERROR",
                self,
                message="maximum number of two seqs per collection",
                source=seqs,
            )
        seqs = seqs.to_moltype(self.moltype)
        seq1, seq2 = seqs.seqs
        aln, score = classic_align_pairwise(
            seq1,
            seq2,
            self._score_matrix,
            self._insertion_penalty,
            self._extension_penalty,
            True,
            return_score=True,
        )

        aln.info["align_params"] = dict(
            score_matrix=self._score_matrix,
            insertion_penalty=self._insertion_penalty,
            extension_penalty=self._extension_penalty,
            sw_score=score,
        )
        return aln.to_moltype(self.moltype)


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

        Examples
        --------

        Create a sample alignment and compute the Information Content alignment
        quality score. The default is equally frequent motif probabilities.

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> aln = make_aligned_seqs(
        ...     {"s1": "AATTGA", "s2": "AGGTCC", "s3": "AGGATG", "s4": "AGGCGT"}
        ... )
        >>> app = get_app("ic_score")
        >>> result = app(aln)
        >>> print(result)
        5.377443751081734
        """
        self._equi_frequent = equifreq_mprobs

    def main(self, aln: AlignedSeqsType) -> float:
        counts = aln.counts_per_pos(include_ambiguity=False, allow_gap=False)
        if counts.array.max() == 0 or aln.num_seqs == 1:
            msg = "zero length" if len(aln) == 0 else "single sequence"
            return NotCompleted(
                "FAIL",
                self,
                f"cannot compute alignment quality because {msg}",
                source=aln.info,
            )

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

    Examples
    --------

    Prepare an alignment of BRCA1 genes from 4 species. Create an unaligned
    sequence collection, guide tree, and an app for alignment using cogent3
    ``progressive_align()``.

    >>> from cogent3 import make_unaligned_seqs, get_app
    >>> aln = make_unaligned_seqs(
    ...     {
    ...         "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    ...         "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    ...         "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    ...         "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
    ...     },
    ...     moltype="dna",
    ... )
    >>> newick = "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus:0.06," "Human:0.0):0.04);"
    >>> aligner = get_app("progressive_align", model="HKY85")

    Create a composable app that aligns the sequences and returns the
    cogent3 log-likelihood alignment score.

    >>> scorer = get_app("cogent3_score")
    >>> app = aligner + score
    >>> result = app(aln)
    >>> print(result)
    -130.47375615734916
    """
    if aln.num_seqs == 1 or len(aln) == 0:
        msg = "zero length" if len(aln) == 0 else "single sequence"
        return NotCompleted(
            "FAIL",
            "cogent3_score",
            f"cannot compute alignment quality because {msg}",
            source=aln.info,
        )

    align_params = aln.info.get("align_params", {})
    return align_params.get(
        "lnL",
        NotCompleted(
            "FAIL", "cogent3_score", "no alignment quality score", source=aln.info
        ),
    )


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

        Examples
        --------

        Create a sample alignment and an app to calculate the Sum of Pairs
        alignment score with ``calc="pdist"`` and no gap penalties.

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> aln = make_aligned_seqs({"s1": "AAGAA-A", "s2": "-ATAATG", "s3": "C-TGG-G"})
        >>> app = get_app("sp_score", calc="pdist", gap_extend=0, gap_insert=0)
        >>> result = app(aln)
        >>> print(result)
        5.0

        Penalise gap extensions with ``gap_extend=1`` and insertions with
        ``gap_insert=2``.

        >>> app_gap_penalty = get_app(
        ...     "sp_score", calc="pdist", gap_extend=1, gap_insert=2
        ... )
        >>> result = app_gap_penalty(aln)
        >>> print(result)
        -13.0
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
            msg = "zero length" if len(aln) == 0 else "single sequence"
            return NotCompleted(
                "FAIL",
                self,
                f"cannot compute alignment quality because {msg}",
                source=aln.info,
            )

        self._calc(aln, show_progress=False)
        dmat = self._calc.dists

        if isnan(dmat.array).sum():
            nans = isnan(dmat.array).sum(axis=0) > 0
            names = ", ".join(f"{n!r}" for n in array(dmat.names)[nans])
            return NotCompleted(
                "ERROR",
                self,
                f"Some genetic distances involving {names} were NaN's. "
                "Using calc='pdist' will prevent this issue.",
                source=aln,
            )

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

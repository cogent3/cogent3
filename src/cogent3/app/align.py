import warnings

from bisect import bisect_left

from cogent3 import make_tree
from cogent3.align import (
    global_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.align.progressive import TreeAlign
from cogent3.app import dist
from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.location import (
    _gap_insertion_data,
    _gap_pos_to_map,
    _interconvert_seq_aln_coords,
    _merged_gaps,
)
from cogent3.core.moltype import get_moltype
from cogent3.evolve.models import get_model

from .composable import (
    ALIGNED_TYPE,
    SEQUENCE_TYPE,
    SERIALISABLE_TYPE,
    ComposableSeq,
    NotCompleted,
)
from .tree import quick_tree, scale_branches


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2021.5.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


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
        l = 0
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


def _map_ref_gaps_to_seq(all_ref_seq, curr_ref, other):
    """inject ref gaps missing from curr_ref into other map

    Returns
    -------
    Map
    """
    all_gap_pos, _ = _gap_insertion_data(all_ref_seq)
    c_gpos, c_offsets = _gap_insertion_data(curr_ref)
    o_gpos, o_offsets = _gap_insertion_data(other)
    unique_ref_gaps = set(all_gap_pos) - set(c_gpos)
    if not unique_ref_gaps:  # other seq map unchanged
        return other.map

    # we add the offset from the full ref up to the new gap
    new_gaps = []
    for unique in unique_ref_gaps:
        # convert ref seq coord into current alignment coord
        aln_pos = _interconvert_seq_aln_coords(
            c_gpos, c_offsets, unique[0], seq_pos=True
        )
        # convert current alignment coord into query seq coord
        seq_pos = _interconvert_seq_aln_coords(
            o_gpos, o_offsets, aln_pos, seq_pos=False
        )
        new_gaps.append((seq_pos, unique[1]))

    # because we combine gap positions from two different sequences,
    # the reference and the other, we need to merge overlapping gaps by
    # summation, rather than max
    o_gpos = _merged_gaps(o_gpos, new_gaps, function="sum")
    gap_pos, gap_lengths = list(zip(*o_gpos))
    seq = other.data
    return _gap_pos_to_map(gap_pos, gap_lengths, len(seq))


class align_to_ref(ComposableSeq):
    """Aligns sequences to a nominated reference in the unaligned collection.
    This is much faster, and requires much less memory, than progressive_align
    but the quality will likely be lower. Alignment quality will be strongly
    affected by choice of reference.

    Returns
    -------
    ArrayAlignment.
    """

    _input_types = SEQUENCE_TYPE
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = "SequenceCollection"

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
        super(align_to_ref, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )
        self._formatted_params()
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
            self.func = self.align_to_longest
        else:
            self.func = self.align_to_named_seq
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
        ref_gaps = []
        kwargs = self._kwargs.copy()
        pwise = []
        for seq in seqs.seqs:
            if seq.name == self._ref_name:
                continue

            aln = global_pairwise(ref_seq, seq, **kwargs).to_type(array_align=False)
            pwise.append(((seq.name, aln)))
            gap_pos, _ = _gap_insertion_data(aln.named_seqs[ref_seq.name])
            ref_gaps = _merged_gaps(ref_gaps, gap_pos, function="max")

        gap_pos, gap_lengths = list(zip(*ref_gaps)) if ref_gaps else ([], [])
        m = _gap_pos_to_map(gap_pos, gap_lengths, len(ref_seq))
        aligned = [Aligned(m, ref_seq)]
        for other_name, aln in pwise:
            curr_ref = aln.named_seqs[ref_seq.name]
            other_seq = aln.named_seqs[other_name]
            m = _map_ref_gaps_to_seq(aligned[0], curr_ref, other_seq)
            aligned.append(
                other_seq if m is other_seq.map else Aligned(m, other_seq.data)
            )

        # default to ArrayAlign
        return Alignment(aligned, moltype=self._moltype, info=seqs.info).to_type(
            array_align=True, moltype=self._moltype
        )


class progressive_align(ComposableSeq):
    """Progressive multiple sequence alignment via any cogent3 model.
    Returns an Alignment object."""

    _input_types = SEQUENCE_TYPE
    _output_types = (ALIGNED_TYPE, SERIALISABLE_TYPE)
    _data_types = "SequenceCollection"

    def __init__(
        self,
        model,
        gc=None,
        param_vals=None,
        guide_tree=None,
        unique_guides=False,
        indel_length=1e-1,
        indel_rate=1e-10,
        distance="percent",
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
        super(progressive_align, self).__init__(
            input_types=self._input_types,
            output_types=self._output_types,
            data_types=self._data_types,
        )

        self._param_vals = {
            "codon": dict(omega=0.4, kappa=3),
            "nucleotide": dict(kappa=3),
        }.get(model, param_vals)
        sm = {"codon": "MG94HKY", "nucleotide": "HKY85", "protein": "JTT92"}.get(
            model, model
        )
        self._formatted_params()
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

        self.func = self.multiple_align

    def _build_guide(self, seqs):
        tree = self._make_tree(seqs)
        if self._scalar != 1:
            scaler = scale_branches(scalar=self._scalar)
            tree = scaler(tree)
        return tree

    def multiple_align(self, seqs):
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        if self._guide_tree is None or self._unique_guides:
            self._guide_tree = self._build_guide(seqs)
            self._kwargs["tree"] = self._guide_tree
            diff = set(self._guide_tree.get_tip_names()) ^ set(seqs.names)
            if diff:
                numtips = len(set(self._guide_tree.get_tip_names()))
                print(f"numseqs={len(seqs.names)} not equal " f"to numtips={numtips}")
                print(f"These were different: {diff}")
                seqs = seqs.take_seqs(self._guide_tree.get_tip_names())

        kwargs = self._kwargs.copy()

        with warnings.catch_warnings(record=False):
            warnings.simplefilter("ignore")
            try:
                result, tree = TreeAlign(self._model, seqs, **kwargs)
                if self._moltype and self._moltype != result.moltype:
                    result = result.to_moltype(self._moltype)

                result.info.update(seqs.info)
            except ValueError as err:
                # probably an internal stop
                result = NotCompleted("ERROR", self, err.args[0], source=seqs)
                return result
        return result

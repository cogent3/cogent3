import warnings

from cogent3 import make_tree
from cogent3.align import (
    global_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.align.progressive import TreeAlign
from cogent3.app import dist
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


class _GapInRef:
    """assumes first element of series is reference, returns True if that matches
    gap_state"""

    def __init__(self, moltype, gap):
        self.gap_state = moltype.alphabet.to_indices(gap)[0]
        self.func = self._ref_gap if gap == "-" else self._array_ref_gap

    def _ref_gap(self, x):
        return x[0] != self.gap_state

    def _array_ref_gap(self, x):
        return x.flatten()[0] != self.gap_state

    def __call__(self, x):
        return self.func(x)


class align_to_ref(ComposableSeq):
    """Aligns to a reference seq, no gaps in the reference.
    Returns an Alignment object."""

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
        aligned = None
        kwargs = self._kwargs.copy()
        no_ref_gap = None

        for i in range(seqs.num_seqs):
            seq = seqs.seqs[i]
            if seq.name == self._ref_name:
                continue

            result = global_pairwise(ref_seq, seq, **kwargs)
            if no_ref_gap is None:
                no_ref_gap = _GapInRef(result.moltype, seqs.moltype.gap)

            # as we're going to be using a pairwise distance that excludes gaps
            # eliminating positions with deletions in the reference
            result = result.filtered(no_ref_gap)
            aligned = result if aligned is None else aligned.add_from_ref_aln(result)

        # default to ArrayAlign
        new = aligned.to_type(array_align=True, moltype=self._moltype)
        return new


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

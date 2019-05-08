from cogent3 import LoadTree
from cogent3.core.alignment import ArrayAlignment
from cogent3.core.moltype import get_moltype
from cogent3.align import global_pairwise, make_dna_scoring_dict
from cogent3.align.progressive import TreeAlign
from cogent3.evolve.models import get_model, protein_models

from .composable import ComposableSeq

from .tree import quick_tree, scale_branches

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class align_to_ref(ComposableSeq):
    """returns an alignment to a reference seq. No gaps in the reference."""

    def __init__(self, ref_seq='longest', score_matrix=None,
                 insertion_penalty=20, extension_penalty=2, moltype='dna'):
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
            molecular type, currently only DNA or RNA suppported
        """
        super(align_to_ref, self).__init__(input_type='sequences',
                                           output_type='aligned')
        assert moltype
        moltype = get_moltype(moltype)
        self._moltype = moltype
        S = score_matrix or make_dna_scoring_dict(10, -1, -8)
        self._kwargs = dict(S=S, d=insertion_penalty, e=extension_penalty,
                            return_score=False)
        if ref_seq.lower() == 'longest':
            self.func = self.align_to_longest
        else:
            self.func = self.align_to_named_seq
            self._ref_name = ref_seq

    def align_to_longest(self, seqs):
        """returns alignment to longest seq"""
        seqs = seqs.to_moltype(self._moltype)
        lengths = seqs.get_lengths()
        lengths = [(l, n) for n, l in lengths.items()]
        _, ref_seq_name = max(lengths)
        self._ref_name = ref_seq_name
        return self.align_to_named_seq(seqs)

    def align_to_named_seq(self, seqs):
        """returns alignment to named seq"""
        seqs = seqs.to_moltype(self._moltype)
        ref_seq = seqs.get_seq(self._ref_name)
        aligned = None
        kwargs = self._kwargs.copy()

        def gap_in_ref(gap_char):
            gap_char = gap_char[0]

            def array_ref_gap(x):
                r = x.flatten()[0] != gap_char
                return r

            def standard_ref_gap(x):
                r = x[0] != gap_char
                return r

            func = {'-': standard_ref_gap}.get(gap_char, array_ref_gap)
            return func

        no_ref_gap = None

        for i in range(seqs.num_seqs):
            seq = seqs.seqs[i]
            if seq.name == self._ref_name:
                continue

            result = global_pairwise(ref_seq, seq, **kwargs)
            if no_ref_gap is None:
                gap = result.moltype.alphabet.to_indices(seqs.moltype.gap)
                no_ref_gap = gap_in_ref(gap)

            # as we're going to be using a pairwise distance that excludes gaps
            # eliminating positions with deletions in the reference
            result = result.filtered(no_ref_gap)
            if aligned is None:
                aligned = result.to_type(array_align=False)
                continue

            aligned = aligned.add_from_ref_aln(
                result.to_type(array_align=False))

        new = ArrayAlignment(
            aligned.todict(), moltype=seqs.moltype, info=seqs.info)
        return new


class progressive_align(ComposableSeq):
    """returns a multiple sequence alignment."""

    def __init__(self, model, gc=None, param_vals=None, guide_tree=None,
                 unique_guides=False,
                 indel_length=1e-1, indel_rate=1e-10):
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
        """
        super(progressive_align, self).__init__(input_type='sequences',
                                                output_type='aligned')
        if guide_tree is None and model in protein_models + ['protein']:
            raise NotImplementedError('auto-build of guide tree '
                                      'not supported for protein seqs yet')

        self._param_vals = {'codon': dict(omega=0.4, kappa=3),
                            'nucleotide': dict(kappa=3)}.get(model, param_vals)
        sm = {'codon': 'MG94HKY', 'nucleotide': 'HKY85',
              'protein': 'JTT92'}.get(model, model)

        kwargs = {} if gc is None else dict(gc=gc)
        sm = get_model(sm, **kwargs)
        moltype = sm.alphabet.moltype
        self._model = sm
        self._scalar = sm.get_word_length()
        self._indel_length = indel_length
        self._indel_rate = indel_rate
        self._moltype = moltype
        self._unique_guides = unique_guides
        if callable(guide_tree):
            self._make_tree = guide_tree
            guide_tree = None  # callback takes precedence
        else:
            self._make_tree = quick_tree()

        if guide_tree is not None:
            if type(guide_tree) == str:
                guide_tree = LoadTree(treestring=guide_tree)
            # make sure no zero lengths
            guide_tree = scale_branches()(guide_tree)

        self._guide_tree = guide_tree
        self._kwargs = dict(indel_length=self._indel_length,
                            indel_rate=self._indel_rate,
                            tree=self._guide_tree,
                            param_vals=self._param_vals,
                            show_progress=False)

        self.func = self.multiple_align

    def _build_guide(self, seqs):
        crude_aligner = align_to_ref(moltype=self._moltype)
        aln = crude_aligner(seqs)
        tree = self._make_tree(aln)
        if self._scalar != 1:
            scaler = scale_branches(scalar=self._scalar)
            tree = scaler(tree)
        return tree

    def multiple_align(self, seqs):
        seqs = seqs.to_moltype(self._moltype)
        if self._guide_tree is None or self._unique_guides:
            self._guide_tree = self._build_guide(seqs)
            self._kwargs['tree'] = self._guide_tree
            diff = set(self._guide_tree.get_tip_names()) ^ set(seqs.names)
            if diff:
                numtips = len(set(self._guide_tree.get_tip_names()))
                print(f"numseqs={len(seqs.names)} not equal "
                      "to numtips={numtips}")
                print(f"These were different: {diff}")
                seqs = seqs.take_seqs(self._guide_tree.get_tip_names())

        kwargs = self._kwargs.copy()

        try:
            result, tree = TreeAlign(self._model, seqs, **kwargs)
            result = result.to_moltype(self._moltype)
            result.info.update(seqs.info)
        except ValueError:
            # probably an internal stop
            ## TODO need to log this as cause of failure
            result = False
        return result

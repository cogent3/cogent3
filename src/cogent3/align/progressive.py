from typing import Optional, Tuple

from cogent3 import get_app, get_model
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.tree import TreeNode
from cogent3.evolve.distance import EstimateDistances
from cogent3.phylo import nj as NJ
from cogent3.util import progress_display as UI
from cogent3.util import warning as c3warn


@c3warn.deprecated_args(
    "2023.8",
    reason="better name",
    old_new=[("ests_from_pairwise", "params_from_pairwise")],
)
@UI.display_wrap
def tree_align(
    model: str,
    seqs: SequenceCollection,
    tree: Optional[TreeNode] = None,
    indel_rate: float = 0.01,
    indel_length: float = 0.01,
    ui=None,
    params_from_pairwise: bool = True,
    param_vals: dict = None,
    iters: Optional[int] = None,
    approx_dists: bool = True,
) -> Tuple[ArrayAlignment, TreeNode]:
    """Returns a multiple sequence alignment and tree.

    Parameters
    ----------
    model
        a substitution model or the name of one, see available_models()
    seqs
        a sequence collection
    tree
        if None, estimates the guide tree from pairwise distances
    indel_rate, indel_length
        parameters for the progressive pair-HMM
    params_from_pairwise
        if no tree provided and True, the median value
        of the substitution model parameters are used
    param_vals
        named key, value pairs for model parameters. These
        override params_from_pairwise.
    iters
        the number of times the alignment process is repeated. The guide tree
        is updated on each iteration from pairwise distances computed from the
        alignment produced by the previous iteration. If None, does not do any
        iterations.
    approx_dists
        if no guide tree, and model is for DNA / Codons, estimates pairwise
        distances using an approximation and JC69. Otherwise, estimates
        genetic distances from pairwise alignments (which is slower).

    Notes
    -----
    Uses a tree for determining the progressive order. If a tree is not
    provided, a Neighbour Joining tree is constructed from pairwise
    distances. If the model is for DNA, the pairwise distances are from
    SequenceCollection.distance_matrix() for the initial guide tree. For other
    moltypes, and distances are estimated using the provided substitution model
    from pairwise alignments of the sequences.

    Parameters and tree are added to ``<align>.info["align_params"]``.
    """

    _exclude_params = ["mprobs", "rate", "bin_switch"]
    param_vals = dict(param_vals) if param_vals else {}

    model = get_model(model)
    moltype = model.alphabet.moltype

    num_states = len(model.alphabet)

    if isinstance(seqs, SequenceCollection):
        seqs = seqs.to_moltype(moltype)
    else:
        seqs = SequenceCollection(data=seqs, moltype=moltype)

    if tree is not None:
        fix_lengths = get_app("scale_branches", scalar=1.0)
        tip_names = set(tree.get_tip_names())

        seq_names = set(seqs.names)
        assert (
            tip_names == seq_names
        ), f"names don't match between seqs and tree: {tip_names ^ seq_names}"
        tree = tree.bifurcating(name_unnamed=True)
        tree = fix_lengths(tree)
        align = _progressive_hmm(
            indel_length, indel_rate, model, param_vals, seqs, tree
        )

        return align, tree

    if params_from_pairwise:
        est_params = [
            param for param in model.get_param_list() if param not in _exclude_params
        ]
    else:
        est_params = None

    if approx_dists and num_states == 4:
        # we have a nucleic acid alphabet, so we will try the new
        # approximation method
        dmat = seqs.distance_matrix(calc="jc69")
        tree = dmat.quick_tree()
    else:
        # we have to do the pairwise-alignment based approach
        dists, param_vals = _dists_from_pairwise_align(
            est_params, params_from_pairwise, model, param_vals, seqs
        )
        tree = NJ.nj(dists.to_dict())

    tree = tree.bifurcating(name_unnamed=True)
    # makes sure all edges have non-zero length and whether we need to scale
    # the lengths for the codon case
    fix_lengths = get_app("scale_branches", nuc_to_codon=num_states >= 60)
    tree = fix_lengths(tree)

    ui.display("Doing progressive alignment")
    # this is the point at which we do the iterations
    align = _progressive_hmm(
        indel_length, indel_rate, model, {**param_vals}, seqs, tree
    )
    if iters is None:
        return align, tree

    for _ in range(iters):
        dmat = align.distance_matrix(calc="jc69")
        tree = dmat.quick_tree()
        tree = tree.bifurcating(name_unnamed=True)
        tree = fix_lengths(tree)
        align = _progressive_hmm(
            indel_length, indel_rate, model, {**param_vals}, seqs, tree
        )

    return align, tree


def _dists_from_pairwise_align(
    est_params, params_from_pairwise, model, param_vals, seqs
):
    dcalc = EstimateDistances(seqs, model, do_pair_align=True, est_params=est_params)
    dcalc.run()
    if params_from_pairwise and not param_vals:
        # we use the median to avoid the influence of outlier pairs
        param_vals = {}
        for param in est_params:
            numbers = dcalc.get_param_values(param)
            param_vals[param] = numbers.median
    dists = dcalc.get_pairwise_distances()
    return dists, param_vals


def _progressive_hmm(indel_length, indel_rate, model, param_vals, seqs, tree):
    LF = model.make_likelihood_function(tree, aligned=False)
    with LF.updates_postponed():
        for param, val in list(param_vals.items()):
            LF.set_param_rule(param, value=val, is_constant=True)
        LF.set_param_rule("indel_rate", value=indel_rate, is_constant=True)
        LF.set_param_rule("indel_length", value=indel_length, is_constant=True)
        LF.set_sequences(seqs)
    lnL = LF.get_log_likelihood()
    edge = lnL.edge
    try:
        align = edge.get_viterbi_path().get_alignment()
    except ArithmeticError:
        # trying to narrow down conditions for difficult to reproduce exception
        print(
            "###" * 30,
            "",
            tree.get_newick(with_distances=True),
            "",
            "#" * 20,
            "",
            str(LF),
            "",
            "#" * 20,
            "",
            seqs.to_fasta(),
            sep="\n",
        )
        raise

    align = align.to_moltype(model.moltype)
    param_vals.update(
        dict(
            indel_length=indel_length,
            indel_rate=indel_rate,
            guide_tree=tree.get_newick(with_distances=True),
            model=model.name,
            lnL=lnL,
        )
    )
    align.info["align_params"] = param_vals
    return align


def TreeAlign(*args, **kwargs):  # pragma: no cover
    """deprecated, used tree_align()"""
    from cogent3.util.warning import deprecated

    deprecated(
        "function",
        "TreeAlign",
        "tree_align",
        "2023.9",
    )

    return tree_align(*args, **kwargs)

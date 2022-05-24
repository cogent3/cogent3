from cogent3.align import indel_model, pairwise
from cogent3.maths.markov import SiteClassTransitionMatrix
from cogent3.recalculation.definition import (
    CalcDefn,
    CalculationDefn,
    NonParamDefn,
    PartitionDefn,
    ProbabilityParamDefn,
)


__author__ = "Gavin Huttley and Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttleuy"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class IndelParameterDefn(ProbabilityParamDefn):
    # locus means allowed to vary by loci
    valid_dimensions = ("edge", "bin", "locus")
    independent_by_default = False
    default_value = default = 0.4
    lower = 0.0001


def make_indel_model_defn(with_indel_params=True, kn=True):
    if kn:
        klass = indel_model.KnudsenMiyamotoIndelModel
    else:
        klass = indel_model.SimpleIndelModel
    if with_indel_params:
        a = IndelParameterDefn("indel_length")  # P(extend indel)
        r = IndelParameterDefn("indel_rate")  # indels per substitution
        return CalcDefn(klass, name="indels")(r, a)
    else:
        # not optimisable parameter, a constant. Another example is the
        # alignment in an LikFunc
        return NonParamDefn("indel_model")


class FloatWithAttrs(float):
    def __new__(cls, value, **kw):
        return float.__new__(cls, value)

    def __init__(self, value, **kw):
        float.__init__(self)
        for (n, v) in list(kw.items()):
            setattr(self, n, v)


def Edge(seq1, seq2, length, bin_data, switch=1.0, bprobs=None):
    # one sequence pair in, potentialy, a tree
    bins = len(bin_data)
    pair = pairwise.Pair(seq1, seq2)
    EP = pair.make_reversible_emission_probs(
        [(bin.mprobs, bin.Qd) for bin in bin_data], length
    )
    tms = [bin.indel.calc_transition_matrix(length) for bin in bin_data]
    if bins == 1:
        TM = tms[0]
    else:
        assert bprobs
        R = SiteClassTransitionMatrix(switch, bprobs)
        TM = R.nestTransitionMatricies(tms)
    assert min(TM.Matrix.flat) >= 0, bin_data
    return EP.make_pair_HMM(TM)


class BinData(object):
    def __init__(self, mprobs, indel, Qd, rate=1.0):
        self.mprobs = mprobs
        self.indel = indel
        self.Qd = Qd
        self.rate = rate

    def __repr__(self):
        return f"Bin(Pi, Qd, {self.rate}, {vars(self.indel)})"


class AnnotateFloatDefn(CalculationDefn):
    name = "annot"

    def calc(self, value, edge):
        return FloatWithAttrs(value, edge=edge)


class ViterbiPogDefn(CalculationDefn):
    name = "align"

    def calc(self, edge):
        return edge.getaln()


class FwdDefn(CalculationDefn):
    name = "fwd"

    def calc(self, edge):
        return edge.get_forward_score(use_cost_function=False)


class _GetAlign:
    def __init__(self, edge, length1, length2):
        try:
            ratio = length1 / (length1 + length2)
        except (ZeroDivisionError, FloatingPointError):
            ratio = 1.0
        self.edge = edge
        self.ratio = ratio

    def __call__(self):
        return self.edge.get_viterbi_path().get_alignable(self.ratio)


class EdgeSumAndAlignDefn(CalculationDefn):
    name = "pair"

    def calc(self, pog1, pog2, length1, length2, bin):
        edge = Edge(pog1, pog2, length1 + length2, [bin])
        edge.getaln = _GetAlign(edge, length1, length2)
        return edge


class EdgeSumAndAlignDefnWithBins(CalculationDefn):
    name = "pair"

    def calc(self, pog1, pog2, length1, length2, switch, bprobs, *bin_data):
        edge = Edge(pog1, pog2, length1 + length2, bin_data, switch, bprobs)

        def _getaln():
            ratio = length1 / (length1 + length2)
            (vtScore, result) = edge.getViterbiScoreAndAlignable(ratio)
            return result

        edge.getaln = _getaln
        return edge


def _recursive_defns(edge, subst, leaf, edge_defn_constructor, bin_args):
    """A defn which calculates a fwd score with an .edge
    attribute which can provide a viterbi alignment which can be
    provided to a similar defn"""
    scores = []
    args = []
    for child in edge.children:
        if child.istip():
            args.append(leaf.select_from_dimension("edge", child.name))
        else:
            (child_defn, scores2) = _recursive_defns(
                child, subst, leaf, edge_defn_constructor, bin_args
            )
            child_defn = ViterbiPogDefn(child_defn)
            scores.extend(scores2)
            args.append(child_defn)
    child_names = [child.name for child in edge.children]
    assert len(child_names) == 2, child_names
    child_lengths = subst["length"].across_dimension("edge", child_names)
    args.extend(child_lengths)
    args.extend(bin_args)
    edge_defn = edge_defn_constructor(*args)
    # fwd = FwdDefn(edge_defn)
    # scores.append(fwd)
    return (edge_defn, scores)


def make_forward_tree_defn(
    subst_model, tree, bin_names, with_indel_params=True, kn=True
):
    """Pairwise Fwd"""
    indel = make_indel_model_defn(with_indel_params, kn)
    subst = subst_model.make_fundamental_param_controller_defns(bin_names)
    leaf = NonParamDefn("leaf", dimensions=("edge",))

    if len(bin_names) > 1:
        switch = ProbabilityParamDefn("bin_switch", dimensions=["locus"])
        bprobs = PartitionDefn(
            [1.0 / len(bin_names) for bin in bin_names],
            name="bprobs",
            dimensions=["locus"],
            dimension=("bin", bin_names),
        )
        edge_args = [switch, bprobs]
        edge_defn_constructor = EdgeSumAndAlignDefnWithBins
    else:
        edge_args = []
        edge_defn_constructor = EdgeSumAndAlignDefn

    mprobs = subst["word_probs"]
    bin_data = CalcDefn(BinData)(mprobs, indel, subst["Qd"])
    bin_data = bin_data.across_dimension("bin", bin_names)
    edge_args.extend(bin_data)

    (top, scores) = _recursive_defns(
        tree, subst, leaf, edge_defn_constructor, edge_args
    )
    defn = FwdDefn(top)
    # defn = SumDefn(*scores)
    return AnnotateFloatDefn(defn, top)

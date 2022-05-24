#!/usr/bin/env python

import sys  # ,hotshot

from cogent3 import load_aligned_seqs, load_tree
from cogent3.evolve.substitution_model import (
    TimeReversibleCodon,
    TimeReversibleDinucleotide,
    TimeReversibleNucleotide,
)
from cogent3.maths import optimisers
from cogent3.util import parallel


__author__ = "Peter Maxwell and  Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

ALIGNMENT = load_aligned_seqs(filename="data/brca1.fasta")
TREE = load_tree(filename="data/murphy.tree")


def subtree(size):
    names = ALIGNMENT.names[:size]
    assert len(names) == size
    tree = TREE.get_sub_tree(names)  # .balanced()
    return names, tree


def brca_test(subMod, names, tree, length, par_rules, **kw):
    # names = ALIGNMENT.names[:taxa]
    # assert len(names) == taxa
    tree = TREE.get_sub_tree(names)  # .balanced()
    aln = ALIGNMENT.take_seqs(names).omit_gap_pos()[:length]
    assert len(aln) == length, (len(aln), length)
    # the_tree_analysis = LikelihoodFunction(treeobj = tree, submodelobj = subMod, alignobj = aln)
    par_controller = subMod.make_likelihood_function(tree, **kw)
    for par_rule in par_rules:
        par_controller.set_param_rule(**par_rule)
    # lf = par_controller.make_calculator(aln)
    return (par_controller, aln)


def measure_evals_per_sec(pc, aln):
    pc.set_alignment(aln)
    return pc.measure_evals_per_second(time_limit=2.0, wall=False)


def makePC(modelClass, parameterisation, length, taxa, tree, opt_mprobs, **kw):
    modelClass = eval(modelClass)
    if parameterisation is not None:
        predicates = {"silly": silly_predicate}
        par_rules = [{"par_name": "silly", "is_independent": parameterisation}]
    else:
        predicates = {}
        par_rules = []
    subMod = modelClass(
        equal_motif_probs=True,
        optimise_motif_probs=opt_mprobs,
        predicates=predicates,
        recode_gaps=True,
        mprob_model="conditional",
    )
    (pc, aln) = brca_test(subMod, taxa, tree, length, par_rules, **kw)
    return (pc, aln)


def quiet(f, *args, **kw):
    import io
    import sys

    temp = io.StringIO()
    _stdout = sys.stdout
    try:
        sys.stdout = temp
        result = f(*args, **kw)
    finally:
        # pass
        sys.stdout = _stdout
    return result


def evals_per_sec(*args):
    pc, aln = makePC(*args)  # quiet(makeLF, *args)
    speed1 = measure_evals_per_sec(pc, aln)
    speed = str(int(speed1))
    return speed


class CompareImplementations(object):
    def __init__(self, switch):
        self.switch = switch

    def __call__(self, *args):
        self.switch(0)
        (pc, aln) = quiet(makePC, *args)
        speed1 = measure_evals_per_sec(pc, aln)
        self.switch(1)
        (pc, aln) = quiet(makePC, *args)
        speed2 = measure_evals_per_sec(pc, aln)
        if speed1 < speed2:
            speed = f"+{speed2 / speed1:2.1f}"
        else:
            speed = f"-{speed1 / speed2:2.1f}"
        if speed in ["+1.0", "-1.0"]:
            speed = ""
        return speed


def benchmarks(test):
    alphabets = ["Nucleotide", "Dinucleotide", "Codon"]
    sequence_lengths = [18, 2004]
    treesizes = [5, 20]

    for (optimise_motifs, parameterisation) in [
        (False, "global"),
        (False, "local"),
        (True, "global"),
    ]:
        print(parameterisation, ["", "opt motifs"][optimise_motifs])
        print(" " * 14, end=" ")
        wcol = 5 * len(sequence_lengths) + 2
        for alphabet in alphabets:
            print(str(alphabet).ljust(wcol), end=" ")
        print()
        print("%-15s" % "", end=" ")  # "length"
        for alphabet in alphabets:
            for sequence_length in sequence_lengths:
                print("%4s" % sequence_length, end=" ")
            print("  ", end=" ")
        print()
        print(
            " " * 12
            + (
                " | ".join(
                    [""]
                    + ["-" * (len(sequence_lengths) * 5) for alphabet in alphabets]
                    + [""]
                )
            )
        )
        for treesize in treesizes:
            print(("%4s taxa    | " % treesize), end=" ")
            (taxa, tree) = subtree(treesize)
            for alphabet in alphabets:
                for sequence_length in sequence_lengths:
                    speed = test(
                        alphabet,
                        parameterisation == "local",
                        sequence_length,
                        taxa,
                        tree,
                        optimise_motifs,
                    )
                    print("%4s" % speed, end=" ")
                print("| ", end=" ")
            print()
        print()
    print()


def silly_predicate(a, b):
    return a.count("A") > a.count("T") or b.count("A") > b.count("T")


# def asym_predicate((a,b)):
#    print a, b, 'a' in a
#    return 'a' in a
# mA = Codon()
# mA.setPredicates({'asym': asym_predicate})


def exponentiator_switch(switch):
    import cogent3.evolve.substitution_calculation

    cogent3.evolve.substitution_calculation.use_new = switch


if "relative" in sys.argv:
    test = CompareImplementations(exponentiator_switch)
else:
    test = evals_per_sec

parallel.inefficiency_forgiven = True

if parallel.get_rank() > 0:
    # benchmarks(test)
    quiet(benchmarks, test)
else:
    try:
        benchmarks(test)
    except KeyboardInterrupt:
        print(" OK")

#!/usr/bin/env python

"""testing the alignment simulation code. We will first just create a simple
Jukes Cantor model using a four taxon tree with very different branch lengths,
 and a Kimura two (really one) parameter model.

The test is to reestimate the parameter values as accurately as possible."""
from cogent3 import make_tree
from cogent3.evolve import substitution_model


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def _est_simulations():
    # specify the 4 taxon tree, and a 'dummy' alignment
    t = make_tree(treestring="(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1)root;")

    # how long the simulated alignments should be
    # at 1000000 the estimates get nice and close
    length_of_align = 10000

    #########################
    #
    # For a Jukes Cantor model
    #
    #########################

    sm = substitution_model.TimeReversibleNucleotide()
    lf = sm.make_likelihood_function(t)
    lf.set_constant_lengths()
    lf.set_name("True JC model")
    print(lf)
    simulated = lf.simulate_alignment(sequence_length=length_of_align)
    print(simulated)

    new_lf = sm.make_likelihood_function(t)
    new_lf = new_lf.set_alignment(simulated)
    new_lf.optimise(tolerance=1.0)
    new_lf.optimise(local=True)
    new_lf.set_name("True JC model")
    print(new_lf)

    #########################
    #
    # a Kimura model
    #
    #########################

    # has a ts/tv term, different values for every edge
    sm = substitution_model.TimeReversibleNucleotide(predicates={"kappa": "transition"})
    lf = sm.make_likelihood_function(t)
    lf.set_constant_lengths()
    lf.set_param_rule("kappa", is_constant=True, value=4.0, edge_name="a")
    lf.set_param_rule("kappa", is_constant=True, value=0.5, edge_name="b")
    lf.set_param_rule("kappa", is_constant=True, value=0.2, edge_name="c")
    lf.set_param_rule("kappa", is_constant=True, value=3.0, edge_name="d")
    lf.set_param_rule("kappa", is_constant=True, value=2.0, edge_name="edge.0")
    lf.set_name("True Kappa model")
    print(lf)
    simulated = lf.simulate_alignment(sequence_length=length_of_align)
    print(simulated)
    new_lf = sm.make_likelihood_function(t)
    new_lf.set_param_rule("kappa", is_independent=True)
    new_lf.set_alignment(simulated)
    new_lf.optimise(tolerance=1.0)
    new_lf.optimise(local=True)
    new_lf.set_name("Estimated Kappa model")
    print(new_lf)

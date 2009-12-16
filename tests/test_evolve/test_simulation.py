#!/usr/bin/env python

"""testing the alignment simulation code. We will first just create a simple Jukes Cantor model using a four taxon tree with very different branch lengths, and a Kimura two (really one) parameter model.

The test is to reestimate the parameter values as accurately as possible."""
import sys
from cogent.core import alignment, tree
from cogent.evolve import substitution_model

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# specify the 4 taxon tree, and a 'dummy' alignment
t = tree.LoadTree(treestring='(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1)root;')
#al = alignments.LoadSeqs(data={'a':'a','b':'a','c':'a','d':'a'})

# how long the simulated alignments should be
# at 1000000 the estimates get nice and close
length_of_align = 10000

#########################
#
# For a Jukes Cantor model
#
#########################

sm = substitution_model.Nucleotide()
lf = sm.makeLikelihoodFunction(t)
lf.setConstantLengths()
lf.setName('True JC model')
print lf
simulated = lf.simulateAlignment(sequence_length=length_of_align)
print simulated

new_lf = sm.makeLikelihoodFunction(t)
new_lf = new_lf.setAlignment(simulated)
new_lf.optimise(tolerance=1.0)
new_lf.optimise(local=True)
new_lf.setName('True JC model')
print new_lf

#########################
#
# a Kimura model
#
#########################


# has a ts/tv term, different values for every edge
sm = substitution_model.Nucleotide(predicates={'kappa':'transition'})
lf = sm.makeLikelihoodFunction(t)
lf.setConstantLengths()
lf.setParamRule('kappa',is_const = True, value = 4.0, edge_name='a')
lf.setParamRule('kappa',is_const = True, value = 0.5, edge_name='b')
lf.setParamRule('kappa',is_const = True, value = 0.2, edge_name='c')
lf.setParamRule('kappa',is_const = True, value = 3.0, edge_name='d')
lf.setParamRule('kappa',is_const = True, value = 2.0, edge_name='edge.0')
lf.setName('True Kappa model')
print lf
simulated = lf.simulateAlignment(sequence_length=length_of_align)
print simulated
new_lf = sm.makeLikelihoodFunction(t)
new_lf.setParamRule('kappa',is_independent=True)
new_lf.setAlignment(simulated)
new_lf.optimise(tolerance=1.0)
new_lf.optimise(local=True)
new_lf.setName('Estimated Kappa model')
print new_lf



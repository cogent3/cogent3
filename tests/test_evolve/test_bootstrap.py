#!/usr/bin/env python

import os
import unittest

from cogent3 import load_aligned_seqs, load_tree
from cogent3.evolve import bootstrap, substitution_model


__author__ = "Peter Maxwell and  Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Matthew Wakefield",
    "Helen Lindsay",
    "Andrew Butterfield",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, "data")

seqnames = ["Chimpanzee", "Rhesus", "Orangutan", "Human"]

REPLICATES = 2


def float_ge_zero(num, epsilon=1e-6):
    """compare whether a floating point value is >= zero with epsilon
    tolerance."""
    if num >= 0.0:
        return True
    elif abs(num - 0.0) < epsilon:
        return True
    else:
        return False


class BootstrapTests(unittest.TestCase):
    def gettree(self):
        treeobj = load_tree(filename=os.path.join(data_path, "murphy.tree"))

        return treeobj.get_sub_tree(seqnames)

    def getsubmod(self, choice="F81"):
        if choice == "F81":
            return substitution_model.TimeReversibleNucleotide(model_gaps=True)
        else:
            return substitution_model.TimeReversibleNucleotide(
                model_gaps=True, predicates={"kappa": "transition"}
            )

    def getalignmentobj(self):
        moltype = self.getsubmod().moltype
        alignmentobj = load_aligned_seqs(
            os.path.join(data_path, "brca1.fasta"), moltype=moltype
        )
        return alignmentobj.take_seqs(seqnames)[:1000]

    def getcontroller(self, treeobj, submodobj):
        return submodobj.make_likelihood_function(treeobj)

    def create_null_controller(self, alignobj):
        """A null model controller creator.
        We constrain the human chimp branches to be equal."""
        treeobj = self.gettree()
        submodobj = self.getsubmod()
        controller = self.getcontroller(treeobj, submodobj)

        # we are setting a local molecular clock for human/chimp
        controller.set_local_clock("Human", "Chimpanzee")
        return controller

    def create_alt_controller(self, alignobj):
        """An alternative model controller. Chimp/Human
        branches are free to vary."""

        treeobj = self.gettree()
        submodobj = self.getsubmod()
        controller = self.getcontroller(treeobj, submodobj)
        return controller

    def calclength(self, likelihood_function):
        """This extracts the length of the human branch and returns it."""

        return likelihood_function.get_param_value("length", "Human")

    def test_conf_int(self):
        """testing estimation of confidence intervals."""
        alignobj = self.getalignmentobj()

        bstrap = bootstrap.EstimateConfidenceIntervals(
            self.create_null_controller(alignobj), self.calclength, alignobj
        )
        bstrap.set_num_replicates(REPLICATES)
        bstrap.set_seed(1984)
        bstrap.run(local=True)
        samplelnL = bstrap.get_sample_lnL()
        for lnL in samplelnL:
            assert lnL < 0.0, lnL

        observed_stat = bstrap.get_observed_stats()
        assert float_ge_zero(observed_stat)

        samplestats = bstrap.getSampleStats()

        for stat in samplestats:
            assert float_ge_zero(stat)

        self.assertEqual(len(samplelnL), REPLICATES)
        self.assertEqual(len(samplestats), REPLICATES)

    def test_prob(self):
        """testing estimation of probability."""

        alignobj = self.getalignmentobj()
        prob_bstrap = bootstrap.EstimateProbability(
            self.create_null_controller(alignobj),
            self.create_alt_controller(alignobj),
            alignobj,
        )
        prob_bstrap.set_num_replicates(REPLICATES)
        prob_bstrap.set_seed(1984)
        prob_bstrap.run(local=True)

        self.assertEqual(len(prob_bstrap.get_sample_LR_list()), REPLICATES)

        assert float_ge_zero(prob_bstrap.get_observed_LR())

        # check the returned sample LR's for being > 0.0
        for sample_LR in prob_bstrap.get_sample_LR_list():
            # print sample_LR
            assert float_ge_zero(sample_LR), sample_LR

        # check the returned observed lnL fulfill this assertion too, really
        # testing their order
        null, alt = prob_bstrap.get_observed_lnL()
        assert float_ge_zero(2 * (alt - null))

        # now check the structure of the returned sample
        for snull, salt in prob_bstrap.get_sample_lnL():
            # print salt, snull, 2*(salt-snull)
            assert float_ge_zero(2 * (salt - snull))

        # be sure we get something back from getprob if proc rank is 0
        assert float_ge_zero(prob_bstrap.get_estimated_prob())


if __name__ == "__main__":
    unittest.main()

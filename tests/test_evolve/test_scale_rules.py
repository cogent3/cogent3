#!/usr/bin/env python

import unittest

from cogent3 import make_tree
from cogent3.evolve import substitution_model
from cogent3.evolve.predicate import MotifChange, replacement


def a_c(x, y):
    return (x == "A" and y == "C") or (x == "C" and y == "A")


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

a_c = MotifChange("A", "C")
trans = MotifChange("A", "G") | MotifChange("T", "C")

TREE = make_tree(tip_names="ab")


class ScaleRuleTests(unittest.TestCase):
    def _makeModel(self, predicates, scale_rules=None):
        scale_rules = scale_rules or []
        return substitution_model.TimeReversibleNucleotide(
            equal_motif_probs=True,
            model_gaps=False,
            predicates=predicates,
            scales=scale_rules,
        )

    def _get_scaled_lengths(self, model, params):
        LF = model.make_likelihood_function(TREE)
        for param in params:
            LF.set_param_rule(param, value=params[param], is_constant=True)
        result = {}
        for predicate in model.scale_masks:
            result[predicate] = LF.get_scaled_lengths(predicate)["a"]
        return result

    def test_scaled(self):
        """Scale rule requiring matrix entries to have all pars specified"""
        model = self._makeModel({"k": trans}, {"ts": trans, "tv": ~trans})

        self.assertEqual(
            self._get_scaled_lengths(model, {"k": 6.0, "length": 4.0}),
            {"ts": 3.0, "tv": 1.0},
        )

    def test_binned(self):
        model = self._makeModel({"k": trans}, {"ts": trans, "tv": ~trans})

        LF = model.make_likelihood_function(TREE, bins=2)
        LF.set_param_rule("length", value=4.0, is_constant=True)
        LF.set_param_rule("k", value=6.0, bin="bin0", is_constant=True)
        LF.set_param_rule("k", value=1.0, bin="bin1", is_constant=True)

        for (bin, expected) in [("bin0", 3.0), ("bin1", 4.0 / 3), (None, 13.0 / 6)]:
            self.assertEqual(LF.get_scaled_lengths("ts", bin=bin)["a"], expected)

    def test_scaled_or(self):
        """Scale rule where matrix entries can have any of the pars specified"""
        model = self._makeModel(
            {"k": trans, "ac": a_c}, {"or": (trans | a_c), "not": ~(trans | a_c)}
        )

        self.assertEqual(
            self._get_scaled_lengths(model, {"k": 6.0, "length": 6.0, "ac": 3.0}),
            {"or": 5.0, "not": 1.0},
        )

    def test_scaling(self):
        """Testing scaling calculations using Dn and Ds as an example."""
        model = substitution_model.TimeReversibleCodon(
            model_gaps=False,
            recode_gaps=True,
            predicates={"k": trans, "r": replacement},
            motif_probs={
                "TAT": 0.0088813702685557206,
                "TGT": 0.020511736096426307,
                "TCT": 0.024529498836963416,
                "TTT": 0.019454430112074435,
                "TGC": 0.0010573059843518714,
                "TGG": 0.0042292239374074857,
                "TAC": 0.002326073165574117,
                "TTC": 0.0086699090716853451,
                "TCG": 0.0010573059843518714,
                "TTA": 0.020723197293296681,
                "TTG": 0.01036159864664834,
                "TCC": 0.0082469866779445976,
                "TCA": 0.022414886868259674,
                "GCA": 0.015648128568407697,
                "GTA": 0.014590822584055826,
                "GCC": 0.0095157538591668436,
                "GTC": 0.0063438359061112285,
                "GCG": 0.0016916895749629942,
                "GTG": 0.0067667582998519769,
                "CAA": 0.018185662930852189,
                "GTT": 0.021569042080778176,
                "GCT": 0.014167900190315077,
                "ACC": 0.0042292239374074857,
                "GGT": 0.014167900190315077,
                "CGA": 0.0012687671812222456,
                "CGC": 0.0010573059843518714,
                "GAT": 0.030238951152463524,
                "AAG": 0.034891097483611758,
                "CGG": 0.002326073165574117,
                "ACT": 0.028758722774370905,
                "GGG": 0.0071896806935927262,
                "GGA": 0.016282512159018821,
                "GGC": 0.0090928314654260944,
                "GAG": 0.031296257136815393,
                "AAA": 0.05476844998942694,
                "GAC": 0.011207443434129837,
                "CGT": 0.0033833791499259885,
                "GAA": 0.076337492070205112,
                "CTT": 0.010573059843518714,
                "ATG": 0.012687671812222457,
                "ACA": 0.021991964474518927,
                "ACG": 0.00084584478748149711,
                "ATC": 0.0076126030873334746,
                "AAC": 0.022837809262000422,
                "ATA": 0.017762740537111441,
                "AGG": 0.013533516599703954,
                "CCT": 0.025586804821315288,
                "AGC": 0.029393106364982026,
                "AGA": 0.021991964474518927,
                "CAT": 0.021357580883907802,
                "AAT": 0.05772890674561218,
                "ATT": 0.019031507718333687,
                "CTG": 0.012899133009092831,
                "CTA": 0.013744977796574329,
                "CTC": 0.0078240642842038483,
                "CAC": 0.0050750687248889825,
                "CCG": 0.00021146119687037428,
                "AGT": 0.03742863184605625,
                "CAG": 0.024106576443222668,
                "CCA": 0.021357580883907802,
                "CCC": 0.0069782194967223515,
            },
            scales={"dN": replacement, "dS": ~replacement},
            mprob_model="tuple",
        )
        length = 0.1115

        a = self._get_scaled_lengths(
            model, {"k": 3.6491, "r": 0.6317, "length": length}
        )
        b = self._get_scaled_lengths(model, {"k": 3.6491, "r": 1.0, "length": length})
        dN = length * a["dN"] / (3.0 * b["dN"])
        dS = length * a["dS"] / (3.0 * b["dS"])
        # following are results from PAML
        self.assertEqual(f"{dN:.4f}", "0.0325")
        self.assertEqual(f"{dS:.4f}", "0.0514")


if __name__ == "__main__":
    unittest.main()

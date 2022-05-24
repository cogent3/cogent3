#!/usr/bin/env Python
"""Data for molecular weight calculations on proteins and nucleotides."""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

ProteinWeights = {
    "A": 71.09,
    "C": 103.16,
    "D": 115.10,
    "E": 129.13,
    "F": 147.19,
    "G": 57.07,
    "H": 137.16,
    "I": 113.18,
    "K": 128.19,
    "L": 113.18,
    "M": 131.21,
    "N": 114.12,
    "P": 97.13,
    "Q": 128.15,
    "R": 156.20,
    "S": 87.09,
    "T": 101.12,
    "V": 99.15,
    "W": 186.23,
    "Y": 163.19,
    "U": 150.06,
}

RnaWeights = {"A": 313.21, "U": 290.17, "C": 289.19, "G": 329.21}

DnaWeights = {"A": 297.21, "T": 274.17, "C": 273.19, "G": 313.21}

ProteinWeightCorrection = 18.0  # terminal residues not dehydrated
DnaWeightCorrection = 61.96  # assumes 5' monophosphate, 3' OH
RnaWeightCorrection = DnaWeightCorrection


class WeightCalculator(object):
    """Calculates molecular weight of a non-degenerate sequence."""

    def __init__(self, Weights, Correction):
        """Returns a new WeightCalculator object (class, so serializable)."""
        self.Weights = Weights
        self.Correction = Correction

    def __call__(self, seq, correction=None):
        """Returns the molecular weight of a specified sequence."""
        if not seq:
            return 0
        if correction is None:
            correction = self.Correction
        get_mw = self.Weights.get
        return sum([get_mw(i, 0) for i in seq]) + correction


DnaMW = WeightCalculator(DnaWeights, DnaWeightCorrection)
RnaMW = WeightCalculator(RnaWeights, DnaWeightCorrection)
ProteinMW = WeightCalculator(ProteinWeights, ProteinWeightCorrection)

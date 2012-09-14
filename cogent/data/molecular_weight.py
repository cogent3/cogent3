#!/usr/bin/env Python
"""Data for molecular weight calculations on proteins and nucleotides."""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

ProteinWeights = {
    'A': 89.09,
    'C': 121.16,
    'D': 133.10,
    'E': 147.13,
    'F': 165.19,
    'G': 75.07,
    'H': 155.16,
    'I': 131.18,
    'K': 146.19,
    'L': 131.18,
    'M': 149.21,
    'N': 132.12,
    'P': 115.13,
    'Q': 146.15,
    'R': 174.20,
    'S': 105.09,
    'T': 119.12,
    'V': 117.15,
    'W': 204.23,
    'Y': 181.19,
    'U': 168.06,
}

RnaWeights = {
    'A': 313.21,
    'U': 290.17,
    'C': 289.19,
    'G': 329.21,
}

DnaWeights = {
    'A': 297.21,
    'T': 274.17,
    'C': 273.19,
    'G': 313.21,
}

ProteinWeightCorrection = 18.0          #terminal residues not dehydrated
DnaWeightCorrection = 61.96             #assumes 5' monophosphate, 3' OH
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

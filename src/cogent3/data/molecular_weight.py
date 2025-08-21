"""Data for molecular weight calculations on proteins and nucleotides."""

from cogent3.util import warning as c3warn

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


class WeightCalculator:
    """Calculates molecular weight of a non-degenerate sequence."""

    # refactor: array
    @c3warn.deprecated_args(
        "2025.9",
        reason="pep8",
        old_new=[("Weights", "weights"), ("Correction", "correction")],
    )
    def __init__(self, weights: dict[str, float], correction: float) -> None:
        """Returns a new WeightCalculator object (class, so serializable)."""
        self.weights = weights
        self.correction = correction

    def __call__(self, seq: str, correction: float | None = None) -> float:
        """Returns the molecular weight of a specified sequence."""
        if not seq:
            return 0
        if correction is None:
            correction = self.correction
        mws = self.weights
        return sum(mws[i] for i in seq) + correction


DnaMW = WeightCalculator(DnaWeights, DnaWeightCorrection)
RnaMW = WeightCalculator(RnaWeights, DnaWeightCorrection)
ProteinMW = WeightCalculator(ProteinWeights, ProteinWeightCorrection)

#!/usr/bin/env python
"""EnergyParams class and instances of entropy and enthalpy params.

The enthalpy and entropy params for the 10 Watson-Crick nearest-neighbor
interactions, initiation corrections, and symmetry corrections
are taken from SantaLucia, PNAS vol 95 1460-1465

GC_INIT is the initiation parameter for duplexes that contain AT LEAST ONE
GC base pair, while AT_INIT is the initiation parameter for duplexes that
contain ONLY AT base pairs. (quoted from SantaLucia, Allawi, and Seneviratne,
Biochemistry, 1996, 3555-3562.)

The SYM term is added into the entropy or enthalpy calculation IFF the
sequence is self-complementary (seq = Reverse Comp of seq)
"""

__author__ = "Amanda Birmingham"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Amanda Birmingham", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Amanda Birmingham"
__email__ = "amanda.birmingham@thermofisher.com"
__status__ = "Production"


class EnergyParams(object):
    """A data structure to hold parameters used in energy calculations"""

    def __init__(self, nearest_neighbor_vals, gc_init, at_init, sym_correct):
        """Store the input params for later reference

        nearest_neighbor_vals: a dictionary or dictionary-castable object
            keyed by each nearest-neighbor pair and holding a value for each.
        gc_init: a floating-point castable value holding the initiation
            parameter for duplexes that contain AT LEAST ONE GC base pair
        at_init: a floating-point castable value holding the initiation
            parameter for duplexes that contain ONLY AT base pairs
        sym_correct: a floating-point castable value that is added into the
            calculation if a sequence is self-complementary
        """

        self.nearestNeighbors = dict(nearest_neighbor_vals)
        self.gcInit = float(gc_init)
        self.atInit = float(at_init)
        self.symCorrection = float(sym_correct)

    # end __init__


# end EnergyParams

# --------------------------
# Enthalpies are in kcal/mol, assumed to be salt concentration independent
_NN_ENTHALPIES = {
    "AA": -7.9,
    "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5,
    "TG": -8.5,
    "GT": -8.4,
    "AC": -8.4,
    "CT": -7.8,
    "AG": -7.8,
    "GA": -8.2,
    "TC": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0,
    "CC": -8.0,
}

_ENTHALPY_GC_INIT = 0.1
_ENTHALPY_AT_INIT = 2.3
_ENTHALPY_SYM = 0
# --------------------------

# --------------------------
# Entropies are in cal/Kelvin*mol, at 1 M NaCl
_NN_ENTROPIES = {
    "AA": -22.2,
    "TT": -22.2,
    "AT": -20.4,
    "TA": -21.3,
    "CA": -22.7,
    "TG": -22.7,
    "GT": -22.4,
    "AC": -22.4,
    "CT": -21.0,
    "AG": -21.0,
    "GA": -22.2,
    "TC": -22.2,
    "CG": -27.2,
    "GC": -24.4,
    "GG": -19.9,
    "CC": -19.9,
}
_ENTROPY_GC_INIT = -2.8
_ENTROPY_AT_INIT = 4.1
_ENTROPY_SYM = -1.4
# --------------------------

# --------------------------
# Module level public EnergyParams instances (one for entropy, one for energy)
ENTHALPY_PARAMS = EnergyParams(
    _NN_ENTHALPIES, _ENTHALPY_GC_INIT, _ENTHALPY_AT_INIT, _ENTHALPY_SYM
)
ENTROPY_PARAMS = EnergyParams(
    _NN_ENTROPIES, _ENTROPY_GC_INIT, _ENTROPY_AT_INIT, _ENTROPY_SYM
)
# --------------------------

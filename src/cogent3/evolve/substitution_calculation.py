#!/usr/bin/env python
import warnings

from numpy.linalg import LinAlgError

from cogent3.maths.matrix_exponentiation import (
    CheckedExponentiator,
    FastExponentiator,
    PadeExponentiator,
)
from cogent3.recalculation.definition import (
    CalculationDefn,
    PositiveParamDefn,
    RatioParamDefn,
)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Custom subclasses of Defn (see cogent3.recalulation) for use by
# substitution models.


class AlignmentAdaptDefn(CalculationDefn):
    name = "leaf_likelihoods"

    def calc(self, model, alignment):
        return model.convert_alignment(alignment)


class LengthDefn(PositiveParamDefn):
    name = "length"
    valid_dimensions = ("edge",)
    independent_by_default = True
    upper = 10.0


class RateDefn(RatioParamDefn):
    name = "rate"
    valid_dimensions = ("bin", "locus")
    independent_by_default = True
    lower = 1e-3
    upper = 1e3


class SubstitutionParameterDefn(RatioParamDefn):
    valid_dimensions = ("edge", "bin", "locus")
    independent_by_default = False


class _EigenPade:
    """class that tries expm via eig first, then Pade if that fails"""

    def __init__(self, eigen):
        self.eigen = eigen
        self.given_expm_warning = False

    def __call__(self, Q):
        try:
            return self.eigen(Q)
        except (ArithmeticError, LinAlgError) as detail:
            if not self.given_expm_warning:
                warnings.warn(f"using slow exponentiator because '{str(detail)}'")
                self.given_expm_warning = True
            return PadeExponentiator(Q)


class ExpDefn(CalculationDefn):
    name = "exp"

    def calc(self, expm):
        (allow_eigen, check_eigen, allow_pade) = {
            "eigen": (True, False, False),
            "checked": (True, True, False),
            "pade": (False, False, True),
            "either": (True, True, True),
        }[str(expm)]

        if not allow_eigen:
            return PadeExponentiator

        eigen = CheckedExponentiator if check_eigen else FastExponentiator

        if not allow_pade:
            return eigen
        else:
            return _EigenPade(eigen=eigen)

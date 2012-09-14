#!/usr/bin/env python
import numpy
import warnings

from cogent.recalculation.definition import PositiveParamDefn, RatioParamDefn, \
        CalculationDefn, MonotonicDefn, ProductDefn, ConstDefn, PartitionDefn, \
        NonParamDefn, CallDefn, SelectForDimension, \
        GammaDefn, WeightedPartitionDefn, CalcDefn
from cogent.maths.matrix_exponentiation import PadeExponentiator, \
        FastExponentiator, CheckedExponentiator, LinAlgError

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Custom subclasses of Defn (see cogent.recalulation) for use by substitution models.

class AlignmentAdaptDefn(CalculationDefn):
    name = 'leaf_likelihoods'
    def calc(self, model, alignment):
        return model.convertAlignment(alignment)       

class LengthDefn(PositiveParamDefn):
    name = 'length'
    valid_dimensions = ('edge',)
    independent_by_default = True
    upper = 10.0

class RateDefn(RatioParamDefn):
    name = 'rate'
    valid_dimensions = ('bin', 'locus')
    independent_by_default = True
    lower = 1e-3
    upper = 1e+3

class SubstitutionParameterDefn(RatioParamDefn):
    valid_dimensions = ('edge', 'bin', 'locus')
    independent_by_default = False

class ExpDefn(CalculationDefn):
    name = 'exp'
    
    def calc(self, expm):
        (allow_eigen, check_eigen, allow_pade) = {
            'eigen': (True, False, False),
            'checked': (True, True, False),
            'pade': (False, False, True),
            'either': (True, True, True),
            }[str(expm)]
        
        if not allow_eigen:
            return PadeExponentiator
        
        eigen = CheckedExponentiator if check_eigen else FastExponentiator
        
        if not allow_pade:
            return eigen
        else:
            def _both(Q, eigen=eigen):
                try:
                    return eigen(Q)
                except (ArithmeticError, LinAlgError), detail:
                    if not _both.given_expm_warning:
                        warnings.warn("using slow exponentiator because '%s'"
                                % str(detail))
                        _both.given_expm_warning = True
                    return PadeExponentiator(Q)
            _both.given_expm_warning = False
            return _both



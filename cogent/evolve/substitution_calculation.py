#!/usr/bin/env python
import numpy

Float = numpy.core.numerictypes.sctype2char(float)

import logging
LOG = logging.getLogger('cogent.evolve.substitution')

from cogent.recalculation.definition import PositiveParamDefn, RatioParamDefn, \
        CalculationDefn, MonotonicDefn, ProductDefn, ConstDefn, PartitionDefn, \
        NonParamDefn, CallDefn, SelectForDimension, SelectFromDimension, \
        GammaDefn, WeightedPartitionDefn, CalcDefn
from cogent.maths.matrix_exponentiation import PadeExponentiator, \
        chooseFastExponentiators, FastExponentiator, LinAlgError

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# These classes make up the algorithm for calculating substitution
# rate matricies.  It is broken up into small steps so that each
# step can be cached.

# Instances can set up working arrays for temp storage within one
# call to calc, but shouldn't assume any state beyond that.  The
# 'args' attribute must be a tuple of CalculationDefns or ParamDefns
# that match the arguments to 'calc'.  See the recalculation
# module for details.

class AlignmentAdaptDefn(CalculationDefn):
    name = 'leaf_likelihoods'
    def calc(self, model, alignment):
        return model.convertAlignment(alignment)

class RateMatrix(object):
    def __init__(self, rate_matrix, ident, symmetric):
        self.rate_matrix = rate_matrix
        self.ident = ident
        self.symmetric = symmetric
    
    def calcQ(self, word_probs, mprobs_matrix):
        Q = self.rate_matrix * mprobs_matrix
        row_totals = numpy.sum(Q, axis=1)
        # Set diagonal values so rows sum to 0
        Q -= self.ident * row_totals
        return Q
    
    def getFastEigenExponentiators(self, word_probs, mprobs_matrix):
        sampleQ = self.calcQ(word_probs, mprobs_matrix)
        return chooseFastExponentiators(sampleQ)
    
    def makeExponentiator(self, motif_probs):
        Q = self.calcQ(motif_probs, motif_probs)
        return FastExponentiator(Q)
    

class ScaledRateMatrix(RateMatrix):
    def calcQ(self, word_probs, mprobs_matrix):
        sum = numpy.sum
        Q = self.rate_matrix * mprobs_matrix
        row_totals = sum(Q, axis=1)
        scale = 1.0 / sum(word_probs * row_totals)
        Q -= self.ident * row_totals
        Q *= scale
        return Q
    
class ExpDefn(CalculationDefn):
    name = 'exp'
    
    def setup(self, model):
        self.model = model
        
    def calc(self, expm):
        (allow_eigen, check_eigen, allow_pade) = {
            'eigen': (True, False, False),
            'checked': (True, True, False),
            'pade': (False, False, True),
            'either': (True, True, True),
            }[str(expm)]
        
        if not allow_eigen:
            return PadeExponentiator
        
        eigen = self.model.suitableEigenExponentiators()[check_eigen]
        
        if not allow_pade:
            return eigen
        else:
            def _both(Q, eigen=eigen):
                try:
                    return eigen(Q)
                except (ArithmeticError, LinAlgError), detail:
                    if not _both.given_expm_warning:
                        LOG.warning("using slow exponentiator because '%s'"
                                % str(detail))
                        _both.given_expm_warning = True
                    return PadeExponentiator(Q)
            _both.given_expm_warning = False
            return _both


class QdDefn(CalculationDefn):
    """Produce an instantaneous rate matrix from the motif probabilities and
    the combined substitution parameter matricies, then diagonalise it"""
    name = 'Qd'
    def setup(self, calc_rate_matrix):
        self.calc_rate_matrix = calc_rate_matrix

    def calc(self, exp, word_probs, mprobs_matrix, *params):
        rate_matrix = self.calc_rate_matrix(*params)
        Q = rate_matrix.calcQ(word_probs, mprobs_matrix)
        return exp(Q)

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


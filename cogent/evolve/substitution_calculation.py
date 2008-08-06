#!/usr/bin/env python
import numpy

from cogent.recalculation.definition import PositiveParamDefn, RatioParamDefn, \
        CalculationDefn, MonotonicDefn, ProductDefn, ConstDefn, PartitionDefn, \
        NonParamDefn, CallDefn, SelectForDimension, GammaDefn, \
        WeightedPartitionDefn
from cogent.maths.matrix_exponentiation import ExponentiatorFactory, \
        chooseFastExponentiator, FastExponentiator

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.1"
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
    
    def calcQ(self, motif_probs):
        Q = self.rate_matrix * motif_probs
        row_totals = numpy.sum(Q, axis=1)
        # Set diagonal values so rows sum to 0
        Q -= self.ident * row_totals
        return Q
    
    def getFastExponentiatorClass(self, motif_probs):
        sampleQ = self.calcQ(motif_probs)
        return chooseFastExponentiator(sampleQ, self.symmetric)
    
    def makeExponentiator(self, motif_probs):
        Q = self.calcQ(motif_probs)
        return FastExponentiator(Q)
    

class ScaledRateMatrix(RateMatrix):
    def calcQ(self, motif_probs):
        sum = numpy.sum
        Q = self.rate_matrix * motif_probs
        row_totals = sum(Q, axis=1)
        scale = 1.0 / sum(motif_probs * row_totals, 0)
        Q -= self.ident * row_totals
        Q *= scale
        return Q
    

class QdDefn(CalculationDefn):
    """Produce an instantaneous rate matrix from the motif probabilities and
    the combined substitution parameter matricies, then diagonalise it"""
    name = 'Qd'
    
    def setup(self, calc_rate_matrix, check_eigenvalues, exponentiator_class):
        self.calc_rate_matrix = calc_rate_matrix
        self.check_eigenvalues = check_eigenvalues
        self.fast_exp = exponentiator_class
    
    def makeCalcFunction(self):
        exp = ExponentiatorFactory(
            check = self.check_eigenvalues,
            exponentiator_class = self.fast_exp)
        calc_rate_matrix = self.calc_rate_matrix
        
        def calc(motif_probs, *params):
            rate_matrix = self.calc_rate_matrix(*params)
            return exp(rate_matrix.calcQ(motif_probs))
        
        return calc
    

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


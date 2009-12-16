#!/usr/bin/env python
# 4 implementations of P = exp(Q*t)
# APIs along the lines of:
#   exponentiator = WhateverExponenentiator(Q or Q derivative(s))
#   P = exponentiator(t)
#
#           Class(Q)     instance(t)       Limitations
# Eigen      slow           fast           not too asymm
# SemiSym    slow           fast           mprobs > 0
# Pade       instant        slow
# Taylor     instant        very slow

from cogent.util.modules import importVersionedModule, ExpectedImportError
import numpy
Float = numpy.core.numerictypes.sctype2char(float)
from numpy.linalg import inv as _inv, eig as _eig,\
                        solve as solve_linear_equations, LinAlgError

import logging
LOG = logging.getLogger('cogent.maths.exponentiation')

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def eig(*a, **kw):
    """make eig to return the same eigvectors as numpy ones"""
    vals, vecs = _eig(*a, **kw)
    return vals, vecs.T
    
def inv(a):
    """make inv return the same contiguous matrix as numpy one"""
    return numpy.ascontiguousarray(_inv(a))

def inv(a):
    """make inv return the same contiguous matrix as numpy one"""
    return numpy.ascontiguousarray(_inv(a))

try:
    pyrex = importVersionedModule('_matrix_exponentiation', globals(),
            (1, 2), LOG, "pure Python/NumPy exponentiation")
except ExpectedImportError:
    pyrex = None
else:
    pyrex.setNumPy(numpy)
    if pyrex.version_info == (1, 2):
        def _pyrex_eigenvectors(q, orig=pyrex.eigenvectors):
            try:
                return orig(q)
            except RuntimeError, detail:
                raise ArithmeticError, detail
        pyrex.eigenvectors = _pyrex_eigenvectors

class _Exponentiator(object):
    def __init__(self, Q):
        self.Q = Q
    
    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, repr(self.Q))
    

class EigenExponentiator(_Exponentiator):
    """A matrix ready for fast exponentiation.  P=exp(Q*t)"""
    
    __slots__ = ['Q', 'ev', 'roots', 'evI', 'evT']
    
    def __init__(self, Q, roots, ev, evI=None):
        self.Q = Q
        if evI is None:
            evI = inv(ev)
        self.evI = evI
        self.evT = numpy.transpose(ev)
        self.ev = ev
        self.roots = roots
    
    def __call__(self, t):
        exp_roots = numpy.exp(t*self.roots)
        result = numpy.inner(self.evT * exp_roots, self.evI)
        if result.dtype.kind == "c":
            result = numpy.asarray(result.real)
        result = numpy.maximum(result, 0.0)
        return result
    

def SemiSymmetricExponentiator(motif_probs, Q, ex):
    """Like EigenExponentiator, but more numerically stable and
    30% faster when the rate matrix (Q/motif_probs) is symmetrical.
    Only usable when all motif probs > 0.  Unlike the others
    it needs to know the motif probabilities."""
    
    H = numpy.sqrt(motif_probs)
    H2 = numpy.divide.outer(H, H)
    #A = Q * H2
    #assert numpy.allclose(A, numpy.transpose(A)), A
    (roots, R) = ex.eigenvectors(Q*H2)
    ev = R / H2
    evI = numpy.transpose(R*H2)
    #self.evT = numpy.transpose(self.ev)
    return ex.exponentiator(Q, roots, ev, evI)


# These next two are slow exponentiators, they don't get any speed up
# from reusing Q with a new t, but for compatability with the diagonalising
# approach they look like they do.  They are derived from code in SciPy.

class TaylorExponentiator(_Exponentiator):
    def __init__(self, Q):
        self.Q = Q
        self.q = 21
    
    def __call__(self, t=1.0):
        """Compute the matrix exponential using a Taylor series of order q."""
        A = self.Q * t
        M = A.shape[0]
        eA = numpy.identity(M, Float)
        trm = eA
        for k in range(1, self.q):
            trm = numpy.dot(trm, A/float(k))
            eA += trm
        while not numpy.allclose(eA, eA-trm):
            k += 1
            trm = numpy.dot(trm, A/float(k))
            eA += trm
        if k >= self.q:
            LOG.warning("Taylor series lengthened from %s to %s" % (self.q, k+1))
            self.q = k + 1
        return eA
    

class PadeExponentiator(_Exponentiator):
    def __init__(self, Q):
        self.Q = Q
    
    def __call__(self, t=1.0):
        """Compute the matrix exponential using Pade approximation of order q.
        """
        A = self.Q * t
        M = A.shape[0]
        # Scale A so that norm is < 1/2
        norm = numpy.maximum.reduce(numpy.sum(numpy.absolute(A), axis=1))
        j = int(numpy.floor(numpy.log(max(norm, 0.5))/numpy.log(2.0))) + 1
        A = A / 2.0**j
        
        # How many iterations required
        e = 1.0
        q = 0
        qf = 1.0
        while e > 1e-12:
            q += 1
            q2 = 2.0 * q
            qf *= q**2 / (q2 * (q2-1) * q2 * (q2+1))
            e = 8 * (norm/(2**j))**(2*q) * qf
        
        # Pade Approximation for exp(A)
        X = A
        c = 1.0/2
        N = numpy.identity(M) + c*A
        D = numpy.identity(M) - c*A
        for k in range(2,q+1):
            c = c * (q-k+1) / (k*(2*q-k+1))
            X = numpy.dot(A,X)
            cX = c*X
            N = N + cX
            if not k % 2:
                D = D + cX;
            else:
                D = D - cX;
        F = solve_linear_equations(D,N)
        for k in range(1,j+1):
            F = numpy.dot(F,F)
        return F
    

import time
def _fastest(fs, *args):
    if len(fs) == 1:
        i = 0
    else:
        es = [[] for f in fs]
        samples = 0
        while samples < 10:
            for (f, e) in zip(fs, es):
                t0 = time.time()
                f(*args)
                t1 = time.time()
                e.append(t1-t0)
            samples += 1
        
        m = []
        for e in es:
            e.sort()
            m.append(e[samples/2])
        i = numpy.argmin(m, -1)
    return (fs[i], fs[i](*args))

def _chooseFastExponentiators(Q):
    if pyrex is not None:
        eigen_candidates = [eig, pyrex.eigenvectors]
        (eigenvectors, (roots,ev)) = _fastest(eigen_candidates, Q)
        inverse_candidates = [inv, pyrex.inverse]
        (inverse, evI) = _fastest(inverse_candidates, ev)
        exponentiator_candidates = [EigenExponentiator(Q, roots, ev, evI)]
        if roots.dtype.kind == 'f':
            assert ev.dtype.kind == 'f'
            pyx = pyrex.EigenExponentiator(Q, roots, ev, evI)
            exponentiator_candidates.append(pyx)
        (exponentiator, P) = _fastest(exponentiator_candidates, 1.1)
        FastestExponentiator = type(exponentiator)
    else:
        eigenvectors = eig
        inverse = inv
        FastestExponentiator = EigenExponentiator
    
    # FastestExponentiator may be the pyrex Exponentiator which
    # doesn't cope with complex inputs
    e_class = {'c': EigenExponentiator, 'f':FastestExponentiator}
    
    def Exp(Q, e_class=e_class):
        (roots, ev) = eigenvectors(Q)
        return e_class[roots.dtype.kind](Q, roots, ev, inverse(ev))
    
    def Exp2(Q, e_class=e_class):
        (roots, ev) = eigenvectors(Q)
        evI = inverse(ev)
        reQ = numpy.inner(ev.T * roots, evI).real
        if not numpy.allclose(Q, reQ):
            raise ArithmeticError, "eigen failed precision test"
        return e_class[roots.dtype.kind](Q, roots, ev, evI)
    
    # These function attributes are just for log / debug etc.
    Exp.eigenvectors = eigenvectors
    Exp.inverse = inverse
    Exp.exponentiator = FastestExponentiator
    
    return (Exp, Exp2)

def chooseFastExponentiators(Q):
    ex = _chooseFastExponentiators(Q)
    LOG.info('Strategy for Q Size %s: PyrexEig:%s PyrexInv:%s PyrexExp:%s'
        % (
        Q.shape[0],
        (ex[0].eigenvectors is not eig),
        (ex[0].inverse is not inv),
        (ex[0].exponentiator is not EigenExponentiator)))
    return ex

def FastExponentiator(Q):
    size = Q.shape[0]
    if pyrex is not None and size < 32:
        (roots, ev) = pyrex.eigenvectors(Q)
    else:
        (roots, ev) = eig(Q)
    ex = EigenExponentiator(Q, roots, ev)
    return ex

def RobustExponentiator(Q):
    return PadeExponentiator(Q)


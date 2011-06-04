"""P matrices for some DNA models can be calculated without going via the 
intermediate rate matrix Q.  A Cython implementation of this calculation can 
be used when Q is not required, for example during likelihood tree optimisation.
Equivalent pure python code is NOT provided because it is typically slower 
than the rate-matrix based alternative and provides no extra functionality.
"""

from cogent.evolve.substitution_model import Nucleotide, CalcDefn
from cogent.evolve.predicate import MotifChange
from cogent.maths.matrix_exponentiation import FastExponentiator
import numpy

from cogent.util.modules import importVersionedModule, ExpectedImportError
try:
    _solved_models = importVersionedModule('_solved_models', globals(), 
            (1, 0), "only matrix exponentiating DNA models")
except ExpectedImportError:
    _solved_models = None
 
class PredefinedNucleotide(Nucleotide):
    _default_expm_setting = None

    # Instead of providing calcExchangeabilityMatrix this subclass overrrides
    # makeContinuousPsubDefn to bypass the Q / Qd step.

    def makeContinuousPsubDefn(self, word_probs, mprobs_matrix, distance, rate_params):
        # Only one set of mprobs will be used
        assert word_probs is mprobs_matrix
        # Order of bases is assumed later, so check it really is Y,Y,R,R:
        alphabet = self.getAlphabet()
        assert set(list(alphabet)[:2]) == set(['T', 'C'])
        assert set(list(alphabet)[2:]) == set(['G', 'A'])
        # Should produce the same P as an ordinary Q based model would:
        self.checkPsubCalculationsMatch()
        return CalcDefn(self.calcPsubMatrix, name='psubs')(
                    word_probs, distance, *rate_params)

    def calcPsubMatrix(self, pi, time, kappa_y=1.0, kappa_r=None):
        """Is F81, HKY83 or TN93 when passed 0, 1 or 2 parameters"""
        if kappa_r is None:
            kappa_r = kappa_y
        result = numpy.empty([4,4], float)
        _solved_models.calc_TN93_P(self._do_scaling, pi, time, kappa_y, kappa_r, result)
        return result
        
    def checkPsubCalculationsMatch(self):
        pi = numpy.array([.1, .2, .3, .4])
        params = [4,6][:len(self.parameter_order)]
        Q = self.calcQ(pi, pi, *params)
        P1 = FastExponentiator(Q)(.5)
        P2 = self.calcPsubMatrix(pi, .5, *params)
        assert numpy.allclose(P1, P2)

def _solvedNucleotide(name, predicates, rate_matrix_required=True, **kw):
    if _solved_models is not None and not rate_matrix_required:
        klass = PredefinedNucleotide
    else:
        klass = Nucleotide
    return klass(name=name, predicates=predicates, model_gaps=False, **kw)

kappa_y = MotifChange('T', 'C').aliased('kappa_y')
kappa_r = MotifChange('A', 'G').aliased('kappa_r')
kappa = (kappa_y | kappa_r).aliased('kappa')

def TN93(**kw):
    """Tamura and Nei 1993 model"""
    return _solvedNucleotide('TN93', [kappa_y, kappa_r], recode_gaps=True, **kw)

def HKY85(**kw):
    """Hasegawa, Kishino and Yanamo 1985 model"""
    return _solvedNucleotide('HKY85', [kappa], recode_gaps=True, **kw)

def F81(**kw):
    """Felsenstein's 1981 model"""
    return _solvedNucleotide('F81', [], recode_gaps=True, **kw)



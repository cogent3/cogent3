#!/usr/bin/env python

import unittest, time, sys, os, numpy
from cogent.maths.optimisers import SimulatedAnnealing, Powell

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class NullFile(object):
    def write(self, x):
        pass
        
def quiet(f, *args, **kw):
    orig = sys.stdout
    try:
        sys.stdout = NullFile()
        result = f(*args, **kw)
    finally:
        sys.stdout = orig
    return result
    
def easilyOptimised(direction):
    def testoptparvector(vector):
        if 0.0 <= vector[0] <= 10.0:
            return -direction * (vector[0]-5.0)**2 + 3
        else:
            return -direction * numpy.inf
    bounds = ([0.0], [10.0])
    return (testoptparvector, [1.0], bounds)
    
class OptimiserTestCase(unittest.TestCase):
    def _test_optimisation(self, klass, kw, kw2, accuracy):
        direction = +1
        (f,x,b) = easilyOptimised(direction)
        o = klass(f,x,b, seed=1, direction=direction, **kw)
        t0 = time.time()
        (f, x) = quiet(o.run, **kw2)
        t1 = time.time()
        
        if 'max_evaluations' in kw:
            evals = o.getEvaluationCount()
            #print (evals, kw)
            #self.assertEqual(kw['max_evaluations'], evals)
        #print evals, 'in', t1-t0, '=', evals/(t1-t0), '/sec'
        
        self.assert_(abs(f-3) < .1, (klass.__name__, kw, x, f))
        self.assert_(abs(x[0]-5) < accuracy, (klass.__name__, kw, x, f))
        
    def test_simanneal(self):
        self._test_optimisation(SimulatedAnnealing, {}, {}, 0.0001)
        
    def test_simanneal_limited(self):
        self._test_optimisation(SimulatedAnnealing,
                {'max_evaluations': 50}, {'show_progress':True}, 0.5)

    def test_powell(self):
        self._test_optimisation(Powell, {}, {}, 0.0001)

    def test_powell_limited(self):
        self._test_optimisation(Powell,
                {'max_evaluations': 10}, {'show_progress':True}, 0.1)
                
    def test_checkpointing(self):
        (f,x,b) = easilyOptimised(1)
        filename = 'checkpoint.tmp.pickle'
        if os.path.exists(filename):
            os.remove(filename)
        for proc in [0,1]:
            o = SimulatedAnnealing(f,x,b, seed=1, direction=1)
            o.setCheckpointing(filename)
            quiet(o.run)
        o = SimulatedAnnealing(f,x,b, seed=1, direction=1, init_temp=3.21)
        o.setCheckpointing(filename)
        self.assertRaises(Exception, quiet, o.run)
        if os.path.exists(filename):
            os.remove(filename)
        
        
        
if __name__ == '__main__':
        unittest.main()

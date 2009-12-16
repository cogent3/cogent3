#!/usr/bin/env python
"""
Some tests for the likelihood function class.

tests to do:
    setting of parameters, by coord, by for-all, checking pars sets
    testing the likelihood for specified pars
    getting ancestral probs
    simulating sequence (not possible to verify values as random)
    
    checking that the object resets on tree change, model change, etc
"""
import warnings

warnings.filterwarnings("ignore", "Motif probs overspecified")
warnings.filterwarnings("ignore", "Model not reversible")

import os
from numpy import ones, dot

from cogent.evolve import substitution_model, predicate
from cogent import DNA, LoadSeqs, LoadTree
from cogent.util.unit_test import TestCase, main
from cogent.maths.matrix_exponentiation import PadeExponentiator as expm

Nucleotide = substitution_model.Nucleotide
MotifChange = predicate.MotifChange

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Matthew Wakefield", "Brett Easton"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

ALIGNMENT = LoadSeqs(
    moltype=DNA,
    filename = os.path.join(data_path,'brca1.fasta'))

OTU_NAMES = ["Human", "Mouse", "HowlerMon"]

########################################################
# some funcs for assembling Q-matrices for 'manual' calc


def isTransition(motif1, motif2):
    position = getposition(motif1, motif2)
    a, b = motif1[position], motif2[position]
    transitions = {('A', 'G') : 1, ('C', 'T'):1}
    pair = (min(a, b), max(a, b))
    
    return transitions.has_key(pair)

def numdiffs_position(motif1, motif2):
    assert len(motif1) == len(motif2),\
        "motif1[%s] & motif2[%s] have inconsistent length" %\
        (motif1, motif2)
    
    ndiffs, position = 0, -1
    for i in range(len(motif1)):
        if motif1[i] != motif2[i]:
            position = i
            ndiffs += 1
            
    return ndiffs == 1, position

def isinstantaneous(motif1, motif2):
    if motif1 != motif2 and (motif1 == '-' * len(motif1) or \
                             motif2 == '-' * len(motif1)):
        return True
    ndiffs, position = numdiffs_position(motif1, motif2)
    return ndiffs

def getposition(motif1, motif2):
    ndiffs, position = numdiffs_position(motif1, motif2)
    return position

##############################################################
# funcs for testing the monomer weighted substitution matrices
_root_probs = lambda x: dict([(n1+n2, p1*p2) \
            for n1,p1 in x.items() for n2,p2 in x.items()])

def make_p(length, coord, val):
    """returns a probability matrix with value set at coordinate in
    instantaneous rate matrix"""
    Q = ones((4,4), float)*0.25 # assumes equi-frequent mprobs at root
    for i in range(4):
        Q[i,i] = 0.0
    Q[coord] *= val
    row_sum = Q.sum(axis=1)
    scale = 1/(.25*row_sum).sum()
    for i in range(4):
        Q[i,i] -= row_sum[i]
    Q *= scale
    return expm(Q)(length)


class LikelihoodCalcs(TestCase):
    """tests ability to calculate log-likelihoods for several
    substitution models."""
    def setUp(self):
        self.alignment = ALIGNMENT.takeSeqs(OTU_NAMES)[0: 42]
        self.tree = LoadTree(tip_names=OTU_NAMES)
        self.par_values = {'kappa': 3.0}
        self.length = 1.0
    
    def _makeLikelihoodFunction(self, submod, alignment, **kw):
        calc = submod.makeLikelihoodFunction(self.tree, **kw)
        calc.setAlignment(alignment)
        return calc
    
    def test_no_seq_named_root(self):
        """root is a reserved name"""
        aln = self.alignment.takeSeqs(self.alignment.Names[:4])
        aln = aln.todict()
        one = aln.pop(aln.keys()[0])
        aln["root"] = one
        aln = LoadSeqs(data=aln)
        submod = Nucleotide()
        tree = LoadTree(treestring="%s" % str(tuple(aln.Names)))
        lf = submod.makeLikelihoodFunction(tree)
        try:
            lf.setAlignment(aln)
        except AssertionError:
            pass
        
        collection = aln.degap().NamedSeqs
        collection.pop("Human")
        tree = LoadTree(treestring="%s" % str(tuple(collection.keys())))
        lf = submod.makeLikelihoodFunction(tree, aligned=False)
        try:
            lf.setSequences(collection)
        except AssertionError:
            pass
        
    
    def test_binned_gamma(self):
        """just rate is gamma distributed"""
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='rate', distribution='gamma', mprob_model='tuple')
        lf = self._makeLikelihoodFunction(submod, self.alignment, bins=3)
        try:
            values = lf.getParamValueDict(['bin'])['omega_factor'].values()
        except KeyError:
            # there shouldn't be an omega factor
            pass
        values = lf.getParamValueDict(['bin'])['rate'].values()
        obs = round(sum(values) / len(values), 6)
        self.assertEqual(obs, 1.0)
        self.assertEqual(len(values), 3)
        shape = lf.getParamValue('rate_shape')
    
    def test_binned_gamma_ordered_param(self):
        """rate is gamma distributed omega follows"""
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='rate', partitioned_params='omega', 
            distribution='gamma', mprob_model='tuple')
        lf = self._makeLikelihoodFunction(submod, self.alignment,bins=3) 
        values = lf.getParamValueDict(['bin'])['omega_factor'].values()
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)
        shape = lf.getParamValue('rate_shape')
    
    def test_binned_partition(self):
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='rate', partitioned_params='omega', 
            distribution='free', mprob_model='tuple')
        lf = self._makeLikelihoodFunction(submod, self.alignment, bins=3)
        values = lf.getParamValueDict(['bin'])['omega_factor'].values()
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)
    
    def test_complex_binned_partition(self):
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='kappa', partitioned_params=['omega'], 
            mprob_model='tuple')
        lf = self._makeLikelihoodFunction(submod, self.alignment,
                    bins=['slow', 'fast'])
        lf.setParamRule('kappa', value=1.0, is_const=True)
        lf.setParamRule('kappa', edge="Human", init=1.0, is_const=False)
        values = lf.getParamValueDict(['bin'])['kappa_factor'].values()
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 2)
    
    def test_codon(self):
        """test a three taxa codon model."""
        submod = substitution_model.Codon(
            equal_motif_probs=True,
            do_scaling=False,
            motif_probs=None,
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            mprob_model='tuple')
        
        self.par_values.update({'omega':0.5})
        likelihood_function = self._makeLikelihoodFunction(
                submod, self.alignment)
                    
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertFloatEqual(evolve_lnL, -80.67069614541883)
    
    def test_nucleotide(self):
        """test a nucleotide model."""
        submod = Nucleotide(
            equal_motif_probs=True,
            do_scaling=False,
            motif_probs=None,
            predicates={'kappa': 'transition'})
        # now do using the evolve
        likelihood_function = self._makeLikelihoodFunction(
                submod, self.alignment)
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertFloatEqual(evolve_lnL, -157.49363874840455)
    
    def test_dinucleotide(self):
        """test a dinucleotide model."""
        submod = substitution_model.Dinucleotide(
                equal_motif_probs=True,
                do_scaling=False,
                motif_probs = None,
                predicates = {'kappa': 'transition'},
                mprob_model='tuple')
        likelihood_function = self._makeLikelihoodFunction(
                submod, self.alignment)
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertFloatEqual(evolve_lnL, -102.48145536663735)
    
    def test_protein(self):
        """test a protein model."""
        submod = substitution_model.Protein(
            do_scaling=False, equal_motif_probs=True)
        alignment = self.alignment.getTranslation()
        
        likelihood_function = self._makeLikelihoodFunction(submod, alignment)
        
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertFloatEqual(evolve_lnL, -89.830370754876185)
    

class LikelihoodFunctionTests(TestCase):
    """tests for a tree analysis class. Various tests to create a tree analysis class,
    set parameters, and test various functions.
    """
    def setUp(self):
        self.submodel = Nucleotide(
            do_scaling=True, model_gaps=False, equal_motif_probs=True,
            predicates = {'beta': 'transition'})
        
        self.data = LoadSeqs(
                filename = os.path.join(data_path, 'brca1_5.paml'),
                moltype = self.submodel.MolType)
        
        self.tree = LoadTree(
                filename = os.path.join(data_path, 'brca1_5.tree'))
    
    def _makeLikelihoodFunction(self, **kw):
        lf = self.submodel.makeLikelihoodFunction(self.tree, **kw)
        lf.setParamRule('beta', is_independent=True)
        lf.setAlignment(self.data)
        return lf
    
    def _setLengthsAndBetas(self, likelihood_function):
        for (species, length) in [
                ("DogFaced", 0.1),
                ("NineBande",  0.2),
                ("Human", 0.3),
                ("HowlerMon", 0.4),
                ("Mouse",  0.5)]:
            likelihood_function.setpar("length", length, edge=species)
        for (species1, species2, length) in [
                ("Human", "HowlerMon", 0.7),
                ("Human", "Mouse", 0.6)]:
            LCA = self.tree.getConnectingNode(species1, species2).Name
            likelihood_function.setpar("length", length, edge=LCA)
        
        likelihood_function.setpar("beta", 4.0)
    
    def test_result_str(self):
        # actualy more a test of self._setLengthsAndBetas()
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertEqual(str(likelihood_function), \
"""Likelihood Function Table\n\
======
  beta
------
4.0000
------
=============================
     edge    parent    length
-----------------------------
    Human    edge.0    0.3000
HowlerMon    edge.0    0.4000
   edge.0    edge.1    0.7000
    Mouse    edge.1    0.5000
   edge.1      root    0.6000
NineBande      root    0.2000
 DogFaced      root    0.1000
-----------------------------
===============
motif    mprobs
---------------
    T    0.2500
    C    0.2500
    A    0.2500
    G    0.2500
---------------""")
    
        likelihood_function = self._makeLikelihoodFunction(digits=2,space=2)
        self.assertEqual(str(likelihood_function), \
"""Likelihood Function Table\n\
===============================
     edge  parent  length  beta
-------------------------------
    Human  edge.0    1.00  1.00
HowlerMon  edge.0    1.00  1.00
   edge.0  edge.1    1.00  1.00
    Mouse  edge.1    1.00  1.00
   edge.1    root    1.00  1.00
NineBande    root    1.00  1.00
 DogFaced    root    1.00  1.00
-------------------------------
=============
motif  mprobs
-------------
    T    0.25
    C    0.25
    A    0.25
    G    0.25
-------------""")
    
    def test_calclikelihood(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertAlmostEquals(-250.686745262,
            likelihood_function.testfunction(),places=9)
    
    def test_ancestralsequences(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        result = likelihood_function.reconstructAncestralSeqs()['edge.0']
        a_column_with_mostly_Ts = -1
        motif_G = 2
        self.assertAlmostEquals(2.28460181711e-05,
                result[a_column_with_mostly_Ts][motif_G], places=8)
        lf = self.submodel.makeLikelihoodFunction(self.tree, bins=['low', 'high'])
        lf.setParamRule('beta', bin='low', value=0.1)
        lf.setParamRule('beta', bin='high', value=10.0)
        lf.setAlignment(self.data)
        result = lf.reconstructAncestralSeqs()
    
    def test_likely_ancestral(self):
        """excercising the most likely ancestral sequences"""
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        result = likelihood_function.likelyAncestralSeqs()
    
    def test_simulateAlignment(self):
        "Simulate DNA alignment"
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        simulated_alignment = likelihood_function.simulateAlignment(20, exclude_internal = False)
        self.assertEqual(len(simulated_alignment), 20)
        self.assertEqual(len(simulated_alignment.getSeqNames()), 8)
    
    def test_simulateHetergeneousAlignment(self):
        "Simulate substitution-heterogeneous DNA alignment"
        lf = self.submodel.makeLikelihoodFunction(self.tree, bins=['low', 'high'])
        lf.setParamRule('beta', bin='low', value=0.1)
        lf.setParamRule('beta', bin='high', value=10.0)
        simulated_alignment = lf.simulateAlignment(100)
    
    def test_simulatePatchyHetergeneousAlignment(self):
        "Simulate patchy substitution-heterogeneous DNA alignment"
        lf = self.submodel.makeLikelihoodFunction(self.tree, bins=['low', 'high'], sites_independent=False)
        lf.setParamRule('beta', bin='low', value=0.1)
        lf.setParamRule('beta', bin='high', value=10.0)
        simulated_alignment = lf.simulateAlignment(100)
    
    def test_simulateAlignment2(self):
        "Simulate alignment with dinucleotide model"
        al = LoadSeqs(data={'a':'ggaatt','c':'cctaat'})
        t = LoadTree(treestring="(a,c);")
        sm = substitution_model.Dinucleotide(mprob_model='tuple')
        pc = sm.makeParamController(t)
        lf = pc.makeCalculator(al)
        simalign = lf.simulateAlignment()
        self.assertEqual(len(simalign), 6)
    
    def test_simulateAlignment3(self):
        """Simulated alignment with gap-induced ambiguous positions
        preserved"""
        t = LoadTree(treestring='(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1)root;')
        al = LoadSeqs(data={
            'a':'g--cactat?',
            'b':'---c-ctcct',
            'c':'-a-c-ctat-',
            'd':'-a-c-ctat-'})
        sm = Nucleotide(recode_gaps=True)
        pc = sm.makeParamController(t)
        #pc.setConstantLengths()
        lf=pc.makeCalculator(al)
        #print lf.simulateAlignment(sequence_length=10)
        simulated = lf.simulateAlignment()
        self.assertEqual(len(simulated.getSeqNames()), 4)
        import re
        self.assertEqual(
            re.sub('[ATCG]', 'x', simulated.todict()['a']),
            'x??xxxxxx?')
        
    
    def test_simulateAlignment_root_sequence(self):
        """provide a root sequence for simulating an alignment"""
        def use_root_seq(root_sequence):
            al = LoadSeqs(data={'a':'ggaatt','c':'cctaat'})
            t = LoadTree(treestring="(a,c);")
            sm = substitution_model.Dinucleotide(mprob_model='tuple')
            pc = sm.makeParamController(t)
            lf = pc.makeCalculator(al)
            simalign = lf.simulateAlignment(exclude_internal=False,
                                            root_sequence=root_sequence)
            root = simalign.NamedSeqs['root']
            self.assertEqual(str(root), str(root_sequence))
        
        root_sequence = DNA.makeSequence('GTAATT')
        use_root_seq(root_sequence) # as a sequence instance
        use_root_seq('GTAATC') # as a string
    
    def test_pc_initial_parameters(self):
        """Default parameter values from original annotated tree"""
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        tree = likelihood_function.getAnnotatedTree()
        pc = self.submodel.makeParamController(tree)
        lf = pc.makeCalculator(self.data)
        self.assertEqual(lf.getParamValue("length", "Human"), 0.3)
        self.assertEqual(lf.getParamValue("beta", "Human"), 4.0)
    
    def test_set_par_all(self):
        likelihood_function = self._makeLikelihoodFunction()
        likelihood_function.setpar("length", 4.0)
        likelihood_function.setpar("beta", 6.0)
        self.assertEqual(str(likelihood_function), \
"""Likelihood Function Table
======
  beta
------
6.0000
------
=============================
     edge    parent    length
-----------------------------
    Human    edge.0    4.0000
HowlerMon    edge.0    4.0000
   edge.0    edge.1    4.0000
    Mouse    edge.1    4.0000
   edge.1      root    4.0000
NineBande      root    4.0000
 DogFaced      root    4.0000
-----------------------------
===============
motif    mprobs
---------------
    T    0.2500
    C    0.2500
    A    0.2500
    G    0.2500
---------------""")
        
        #self.submodel.setScaleRule("ts",['beta'])
        #self.submodel.setScaleRule("tv",['beta'], exclude_pars = True)
        self.assertEqual(str(likelihood_function),\
"""Likelihood Function Table
======
  beta
------
6.0000
------
=============================
     edge    parent    length
-----------------------------
    Human    edge.0    4.0000
HowlerMon    edge.0    4.0000
   edge.0    edge.1    4.0000
    Mouse    edge.1    4.0000
   edge.1      root    4.0000
NineBande      root    4.0000
 DogFaced      root    4.0000
-----------------------------
===============
motif    mprobs
---------------
    T    0.2500
    C    0.2500
    A    0.2500
    G    0.2500
---------------""")
    
    def test_getMotifProbs(self):
        likelihood_function = self._makeLikelihoodFunction()
        mprobs = likelihood_function.getMotifProbs()
        assert hasattr(mprobs, 'keys'), mprobs
        keys = mprobs.keys()
        keys.sort()
        obs = self.submodel.getMotifs()
        obs.sort()
        self.assertEqual(obs, keys)
    
    def test_getAnnotatedTree(self):
        likelihood_function = self._makeLikelihoodFunction()
        likelihood_function.setpar("length", 4.0, edge="Human")
        result = likelihood_function.getAnnotatedTree()
        self.assertEqual(result.getNodeMatchingName('Human').params['length'], 4.0)
        self.assertEqual(result.getNodeMatchingName('Human').Length, 4.0)
    
    def test_getstatsasdict(self):
        likelihood_function = self._makeLikelihoodFunction()
        likelihood_function.setName("TEST")
        self.assertEqual(str(likelihood_function),\
"""TEST
=======================================
     edge    parent    length      beta
---------------------------------------
    Human    edge.0    1.0000    1.0000
HowlerMon    edge.0    1.0000    1.0000
   edge.0    edge.1    1.0000    1.0000
    Mouse    edge.1    1.0000    1.0000
   edge.1      root    1.0000    1.0000
NineBande      root    1.0000    1.0000
 DogFaced      root    1.0000    1.0000
---------------------------------------
===============
motif    mprobs
---------------
    T    0.2500
    C    0.2500
    A    0.2500
    G    0.2500
---------------""")
        self.assertEqual(likelihood_function.getStatisticsAsDict(),
{'edge.parent': {'NineBande': 'root', 'edge.1': 'root', 'DogFaced': 'root',
         'Human': 'edge.0', 'edge.0': 'edge.1', 'Mouse': 'edge.1',
         'HowlerMon': 'edge.0'},
 'beta': {'NineBande': 1.0, 'edge.1': 1.0,'DogFaced': 1.0, 'Human': 1.0,
      'edge.0': 1.0, 'Mouse': 1.0, 'HowlerMon': 1.0},
 'length': {'NineBande': 1.0,'edge.1': 1.0, 'DogFaced': 1.0, 'Human': 1.0,
        'edge.0': 1.0, 'Mouse': 1.0,'HowlerMon': 1.0}})
    
    def test_constant_to_free(self):
        """excercise setting a constant param rule, then freeing it"""
        # checks by just trying to make the calculator
        lf = self.submodel.makeLikelihoodFunction(self.tree)
        lf.setAlignment(self.data)
        lf.setParamRule('beta', is_const=True, value=2.0, 
                        edges=['NineBande', 'DogFaced'], is_clade=True)
        lf.real_par_controller.makeCalculator()
        lf.setParamRule('beta', init=2.0, is_const=False,
                        edges=['NineBande', 'DogFaced'], is_clade=True)
        lf.real_par_controller.makeCalculator()
    

if __name__ == '__main__':
    main()

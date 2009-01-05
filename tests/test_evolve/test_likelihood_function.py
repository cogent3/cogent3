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
__version__ = "1.3.0.dev"
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
            ordered_param='rate', distribution='gamma')
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
            ordered_param='rate', partitioned_params='omega', distribution='gamma')
        lf = self._makeLikelihoodFunction(submod, self.alignment,bins=3) 
        values = lf.getParamValueDict(['bin'])['omega_factor'].values()
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)
        shape = lf.getParamValue('rate_shape')
    
    def test_binned_partition(self):
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='rate', partitioned_params='omega', distribution='free')
        lf = self._makeLikelihoodFunction(submod, self.alignment, bins=3)
        values = lf.getParamValueDict(['bin'])['omega_factor'].values()
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)
    
    def test_complex_binned_partition(self):
        submod = substitution_model.Codon(
            predicates={'kappa': 'transition', 'omega': 'replacement'},
            ordered_param='kappa', partitioned_params=['omega'])
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
            do_scaling=False,
            motif_probs=None,
            predicates={'kappa': 'transition', 'omega': 'replacement'})
        
        self.par_values.update({'omega':0.5})
        likelihood_function = self._makeLikelihoodFunction(submod, self.alignment)
                    
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertEqual("%.6f" % -57.8379659216, "%.6f" % evolve_lnL)
    
    def test_nucleotide(self):
        """test a nucleotide model."""
        submod = Nucleotide(
            do_scaling=False,
            motif_probs=None,
            predicates={'kappa': 'transition'})
        # now do using the evolve
        likelihood_function = self._makeLikelihoodFunction(submod, self.alignment)
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertEqual("%.6f" % -155.775725365, "%.6f" % evolve_lnL)
    
    def test_dinucleotide(self):
        """test a dinucleotide model."""
        submod = substitution_model.Dinucleotide(
                do_scaling=False,
                motif_probs = None,
                predicates = {'kappa': 'transition'})
        likelihood_function = self._makeLikelihoodFunction(submod, self.alignment)
        for par, val in self.par_values.items():
            likelihood_function.setpar(par, val)
            
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertEqual("%.6f" % -85.2399172216, "%.6f" % evolve_lnL)
    
    def test_protein(self):
        """test a protein model."""
        submod = substitution_model.Protein(
            do_scaling=False, motif_probs=None)
        alignment = self.alignment.getTranslation()
        
        likelihood_function = self._makeLikelihoodFunction(submod, alignment)
        
        likelihood_function.setpar("length", self.length)
        evolve_lnL = likelihood_function.testfunction()
        self.assertEqual("%.6f" % -76.301896714, "%.6f" % evolve_lnL)
    

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
    
    def _makeLikelihoodFunction(self):
        lf = self.submodel.makeLikelihoodFunction(self.tree)
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
    
    def test_calclikelihood(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertAlmostEquals(-250.686745262,
            likelihood_function.testfunction(),places=9)
    
    def test_ancestralsequences(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        result = likelihood_function.reconstructAncestralSequences()['edge.0']
        a_column_with_mostly_Ts = -1
        motif_G = 2
        self.assertAlmostEquals(2.28460181711e-05,
                result[a_column_with_mostly_Ts][motif_G], places=8)
        lf = self.submodel.makeLikelihoodFunction(self.tree, bins=['low', 'high'])
        lf.setParamRule('beta', bin='low', value=0.1)
        lf.setParamRule('beta', bin='high', value=10.0)
        lf.setAlignment(self.data)
        result = lf.reconstructAncestralSequences()
    
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
        sm = substitution_model.Dinucleotide()
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
            sm = substitution_model.Dinucleotide()
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
    
class NewQ(TestCase):
    aln = LoadSeqs(data={
    'seq1': 'TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
    'seq2': 'TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT'},
    moltype=DNA)
    tree = LoadTree(tip_names=['seq1', 'seq2'])
    
    asymm_nuc_probs_GC = dict(A=0.1,T=0.1,C=0.4,G=0.4)
    asymm_root_probs = _root_probs(asymm_nuc_probs_GC)
    symm_nuc_probs = dict(A=0.25,T=0.25,C=0.25,G=0.25)
    symm_root_probs = _root_probs(symm_nuc_probs)
    
    def test_newQ_is_nuc_process(self):
        """newQ is an extension of an independent nucleotide process"""
        nuc = Nucleotide(motif_probs = self.asymm_nuc_probs_GC)
        new_di = Nucleotide(motif_length=2, mprob_model='monomer',
            motif_probs = self.asymm_nuc_probs_GC)
        
        nuc_lf = nuc.makeLikelihoodFunction(self.tree)
        new_di_lf = new_di.makeLikelihoodFunction(self.tree)
        # newQ branch length is exactly motif_length*nuc branch length
        nuc_lf.setParamRule('length', is_independent=False, init=0.2)
        new_di_lf.setParamRule('length', is_independent=False, init=0.4)
        
        nuc_lf.setAlignment(self.aln)
        new_di_lf.setAlignment(self.aln)
        self.assertFloatEqual(nuc_lf.getLogLikelihood(),
                                new_di_lf.getLogLikelihood())
    
    def test_newQ_is_not_oldQ(self):
        """newQ produces a likelihood different to oldQ when monomers are
        not equi-frequent"""
        new_di = Nucleotide(motif_length=2, mprob_model='monomer',
                                motif_probs = self.asymm_nuc_probs_GC)
        old_di = Nucleotide(motif_length=2, mprob_model=None,
                                motif_probs = self.asymm_root_probs)
        
        new_di_lf = new_di.makeLikelihoodFunction(self.tree)
        old_di_lf = old_di.makeLikelihoodFunction(self.tree)
        
        new_di_lf.setParamRule('length', is_independent=False, init=0.4)
        old_di_lf.setParamRule('length', is_independent=False, init=0.4)
        
        new_di_lf.setAlignment(self.aln)
        old_di_lf.setAlignment(self.aln)
        self.failIfAlmostEqual(new_di_lf.getLogLikelihood(),
                               old_di_lf.getLogLikelihood(), places=2)
    
    def test_newQ_eq_oldQ(self):
        """when monomers are equi-frequent, new and old Q are the same"""
        new_di = Nucleotide(motif_length=2, mprob_model='monomer',
                                motif_probs = self.symm_nuc_probs)
        old_di = Nucleotide(motif_length=2, mprob_model=None,
                                motif_probs = self.symm_root_probs)
        
        new_di_lf = new_di.makeLikelihoodFunction(self.tree)
        old_di_lf = old_di.makeLikelihoodFunction(self.tree)
        
        new_di_lf.setParamRule('length', is_independent=False, init=0.4)
        old_di_lf.setParamRule('length', is_independent=False, init=0.4)
        
        new_di_lf.setAlignment(self.aln)
        old_di_lf.setAlignment(self.aln)
        self.assertFloatEqual(old_di_lf.getLogLikelihood(),
                                new_di_lf.getLogLikelihood())
    
    def test_nuc_models_unaffected(self):
        """nuc models should not be affected by the mprob_model arg"""
        gtr = [MotifChange(x,y) for x,y in 'AC AG AT CG CT'.split()]
        mg_gtr = Nucleotide(predicates=gtr, mprob_model='monomer')
        gy_gtr = Nucleotide(predicates=gtr, mprob_model=None)
        mg_gtr = mg_gtr.makeLikelihoodFunction(self.tree)
        gy_gtr = gy_gtr.makeLikelihoodFunction(self.tree)
        for param, val in zip('A/C A/G A/T C/G C/T'.split(), [.1,2.,.5,.5,3.]):
            for lf in mg_gtr, gy_gtr:
                lf.setParamRule(param, value=val)
        mg_gtr.setAlignment(self.aln)
        gy_gtr.setAlignment(self.aln)
        self.assertFloatEqual(mg_gtr.getLogLikelihood(), gy_gtr.getLogLikelihood())
    
    def test_position_specific_mprobs(self):
        """correctly compute likelihood when positions have distinct
        probabilities"""
        aln_len = len(self.aln)
        posn1 = []
        posn2 = []
        for name, seq in self.aln.todict().items():
            p1 = [seq[i] for i in range(0,aln_len,2)]
            p2 = [seq[i] for i in range(1,aln_len,2)]
            posn1.append([name, ''.join(p1)])
            posn2.append([name, ''.join(p2)])
        
        # the position specific alignments
        posn1 = LoadSeqs(data=posn1)
        posn2 = LoadSeqs(data=posn2)
        
        # NOTE: I had to set do_scaling=False to get this to work
        # NOTE: the argument word_length will be deprecated in a future release
        sm = Nucleotide(word_length=2, do_scaling=False) # a newQ dinucleotide model
        lf = sm.makeLikelihoodFunction(self.tree)
        lf.setAlignment(posn1)
        posn1_lnL = lf.getLogLikelihood()
        lf.setAlignment(posn2)
        posn2_lnL = lf.getLogLikelihood()
        expect_lnL = posn1_lnL+posn2_lnL
        # the joint model
        lf.setAlignment(self.aln)
        aln_lnL = lf.getLogLikelihood()
        
        # setting the full alignment, which has different motif probs, should
        # produce a different lnL
        self.failIfAlmostEqual(expect_lnL, aln_lnL)
        
        # the position specific model
        
        # set the arguments for taking position specific mprobs
        sm = Nucleotide(word_length=2, mprob_model='monomers', do_scaling=False)
        lf = sm.makeLikelihoodFunction(self.tree)
        lf.setAlignment(self.aln)
        posn12_lnL = lf.getLogLikelihood()
        self.assertFloatEqual(expect_lnL, posn12_lnL)
    
    def test_getting_node_mprobs(self):
        """return correct motif probability vector for tree nodes"""
        tree = LoadTree(treestring='(a:.2,b:.2,(c:.1,d:.1):.1)')
        aln = LoadSeqs(data={
        'a': 'TGTG',
        'b': 'TGTG',
        'c': 'TGTG',
        'd': 'TGTG',
        })
        
        motifs = ['T', 'C', 'A', 'G']
        aX = MotifChange(motifs[0], motifs[3], forward_only=True).aliased('aX')
        bX = MotifChange(motifs[3], motifs[0], forward_only=True).aliased('bX')
        edX = MotifChange(motifs[1], motifs[2], forward_only=True).aliased('edX')
        cX = MotifChange(motifs[2], motifs[1], forward_only=True).aliased('cX')
        sm = Nucleotide(predicates=[aX, bX, edX, cX], equal_motif_probs=True)
        
        lf = sm.makeLikelihoodFunction(tree)
        lf.setParamRule('aX', edge='a', value=8.0)
        lf.setParamRule('bX', edge='b', value=8.0)
        lf.setParamRule('edX', edge='edge.0', value=2.0)
        lf.setParamRule('cX', edge='c', value=0.5)
        lf.setParamRule('edX', edge='d', value=4.0)
        lf.setAlignment(aln)
        
        # we construct the hand calc variants
        mprobs = ones(4, float) * .25
        a = make_p(.2, (0,3), 8)
        a = dot(mprobs, a)
        
        b = make_p(.2, (3, 0), 8)
        b = dot(mprobs, b)
        
        e = make_p(.1, (1, 2), 2)
        e = dot(mprobs, e)
        
        c = make_p(.1, (2, 1), 0.5)
        c = dot(e, c)
        
        d = make_p(.1, (1, 2), 4)
        d = dot(e, d)
        
        prob_vectors = lf.getMotifProbsByNode()
        self.assertFloatEqual(prob_vectors['a'].array, a)
        self.assertFloatEqual(prob_vectors['b'].array, b)
        self.assertFloatEqual(prob_vectors['c'].array, c)
        self.assertFloatEqual(prob_vectors['d'].array, d)
        self.assertFloatEqual(prob_vectors['edge.0'].array, e)
    
    def test_get_statistics(self):
        """get statistics should correctly apply arguments"""
        for mprob_model in None, 'monomer', 'monomers':
            di = Nucleotide(motif_length=2, mprob_model=mprob_model)
            lf = di.makeLikelihoodFunction(self.tree)
            for wm, wt in [(True, True), (True, False), (False, True),
                           (False, False)]:
                stats = lf.getStatistics(with_motif_probs=wm, with_titles=wt)
        
    


if __name__ == '__main__':
    main()

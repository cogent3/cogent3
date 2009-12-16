#!/usr/bin/env python

import os

from cogent import LoadSeqs, CodonAlphabet, DNA, LoadTable
from cogent.core import genetic_code
from cogent.evolve import substitution_model, substitution_calculation
from cogent.util.unit_test import TestCase, main

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

class NucleotideModelTestMethods(TestCase):
    def setUp(self):
        self.submodel = substitution_model.Nucleotide(
                do_scaling=True, model_gaps=False)
            
    def test_isTransition(self):
        """testing isTransition"""
        isTransition = self.submodel.getPredefinedPredicate('transition')
        assert isTransition('A', 'G')
        assert isTransition('C', 'T')
        assert not isTransition('A', 'T')
        assert not isTransition('C', 'G')

    def test_isTransversion(self):
        """testing isTransversion"""
        isTransversion = self.submodel.getPredefinedPredicate('transversion')
        assert not isTransversion('A', 'G')
        assert not isTransversion('C', 'T')
        assert isTransversion('A', 'T')
        assert isTransversion('C', 'G')
        
    def test_isIndel(self):
        """testing indel comparison nucleotide model"""
        model = substitution_model.Nucleotide(
                do_scaling=True, model_gaps=True)
        isIndel = model.getPredefinedPredicate('indel')
        assert isIndel('A', '-')
        assert isIndel('-', 'G')
        #assert not self.submodel.isIndel('-', '-')
        assert not isIndel('a', 't')
        
    def test_PredicateChecks(self):
        # overparameterisation
        self.assertRaises(ValueError, substitution_model.Nucleotide,
                model_gaps=False, predicates=['transition', 'transversion'])
        
class MultiLetterMotifSubstModelTests(TestCase):
    def setUp(self):
        self.submodel = substitution_model.Dinucleotide(do_scaling=True, 
                model_gaps=True, mprob_model='tuple')
        
    def test_asciiArt(self):
        model = substitution_model.Dinucleotide(mprob_model='tuple', 
            predicates=['k:transition'])    
        model.asciiArt()
        model = substitution_model.Dinucleotide(mprob_model='tuple')
        model.asciiArt()
        
    def test_isIndel(self):
        """testing indel comparison for dinucleotide model"""
        # these are non-instantaneous
        isIndel = self.submodel.getPredefinedPredicate('indel')
        assert not isIndel('AA', '--')
        assert not isIndel('--', 'CT')

        #assert not self.submodel.isIndel('--', '--')
        assert not isIndel('AT', 'AA')
        
        assert isIndel('AA', 'A-')
        assert isIndel('-A', 'CA')
        assert isIndel('TA', '-A')

        # isIndel can now assume it won't get any non-instantaneous pairs
        # assert self.submodel.isIndel('-a', 'a-') == 0
    

class TupleModelMotifProbFuncs(TestCase):
    dinucs = ('TT', 'CT', 'AT', 'GT',
              'TC', 'CC', 'AC', 'GC',
              'TA', 'CA', 'AA', 'GA',
              'TG', 'CG', 'AG', 'GG')
    nuc_probs = [('T', 0.1), ('C', 0.2), ('A', 0.3), ('G', 0.4)]
    dinuc_probs=[(m2+m1,p1*p2) for m1,p1 in nuc_probs for m2,p2 in nuc_probs]
    mat_indices = dict(
        C=set([(0,1),(0,4),(1,5),(2,1),(2,6),(3,1),(3,7),(4,5),(6,5),
               (7,5),(8,4),(8,9),(9,5),(10,6),(10,9),(11,7),(11,9),(12,4),
               (12,13),(13,5),(14,6),(14,13),(15,7),(15,13)]),
        A=set([(0,2),(0,8),(1,2),(1,9),(2,10),(3,2),(3,11),(4,6),(4,8),
               (5,6),(5,9),(6,10),(7,6),(7,11),(8,10),(9,10),(11,10),
               (12,8),(12,14),(13,9),(13,14),(14,10),(15,11),(15,14)]),
        G=set([(0,3),(0,12),(1,3),(1,13),(2,3),(2,14),(3,15),(4,7),(4,12),
               (5,7),(5,13),(6,7),(6,14),(7,15),(8,11),(8,12),(9,11),
               (9,13),(10,11),(10,14),(11,15),(12,15),(13,15),(14,15)]),
        T=set([(1,0),(2,0),(3,0),(4,0),(5,1),(5,4),(6,2),(6,4),(7,3),
               (7,4),(8,0),(9,1),(9,8),(10,2),(10,8),(11,3),(11,8),(12,0),
               (13,1),(13,12),(14,2),(14,12),(15,3),(15,12)])
           )

class ThreeLetterMotifSubstModelTests(TestCase):
    def setUp(self):
        self.submodel = substitution_model.Nucleotide(motif_length=3,
            mprob_model='tuple')
        
    def test_isIndel(self):
        """testing indel comparison for trinucleotide model"""
        isIndel = self.submodel.getPredefinedPredicate('indel')
        assert isIndel('AAA', 'AA-')
        assert isIndel('-CA', 'CCA')
        assert isIndel('TAC', 'T-C')

        # isIndel can now assume it won't get any non-instantaneous pairs
        assert not isIndel('AAA', '---')
        assert not isIndel('---', 'CTT')
        assert not isIndel('AAA', '--A')
        assert not isIndel('C--', 'CTT')

class CodonSubstModelTests(TestCase):
    def setUp(self):
        self.standardcode = substitution_model.Codon(model_gaps=True, gc=1,
            mprob_model='tuple')
        self.mitocode = substitution_model.Codon(model_gaps=False, gc=2,
            mprob_model='tuple')
        
    def test_isTransition(self):
        """testing codon isTransition"""
        isTransition = self.standardcode.getPredefinedPredicate('transition')
        # first position
        assert isTransition('TGC', 'CGC')
        assert isTransition('GGC', 'AGC')
        # second position   
        assert isTransition('CTT', 'CCT')
        assert isTransition('CAT', 'CGT')
        # thirs position    
        assert isTransition('CTT', 'CTC')
        assert isTransition('CTA', 'CTG')
        # mito code         
        assert isTransition('CTT', 'CTC')
        assert isTransition('CTA', 'CTG')
                            
        assert not isTransition('GAG', 'GTG')
        assert not isTransition('CCC', 'CGC')
        
    def test_isReplacement(self):
        """test isReplacement for the two major genetic codes"""
        isReplacement = self.standardcode.getPredefinedPredicate('replacement')
        # for the standard code, a replacement
        assert isReplacement('CTG', 'ATG')
        assert not isReplacement('AGT','TCC')
        assert not isReplacement('CTG', '---')
        assert not isReplacement('---', 'CTA')
        # for the vert mitocho code, instantaneous replacement
        isReplacement = self.mitocode.getPredefinedPredicate('replacement')
        assert isReplacement('AAA', 'AAC')
        
    def test_isSilent(self):
        """testing isSilent for the two major genetic codes"""
        isSilent = self.standardcode.getPredefinedPredicate('silent')
        assert isSilent('CTA', 'CTG')
        assert not isSilent('AGT','AAG')
        assert not isSilent('CTG', '---')
        assert not isSilent('---', 'CTG')
        # for vert mito code
        isSilent = self.mitocode.getPredefinedPredicate('silent')
        assert isSilent('TCC', 'TCA')
         
    def test_isIndel(self):
        """test isIndel for codon model"""
        isIndel = self.standardcode.getPredefinedPredicate('indel')
        assert isIndel('CTA', '---')
        assert not isIndel('---', '---')
        assert isIndel('---', 'TTC')
                
    def test_str_(self):
        """str() and repr() of a substitution model"""
        s = str(self.standardcode)
        r = repr(self.standardcode)

    
class ModelDataInteractionTestMethods(TestCase):
        
    def test_excludeinggaps(self):
        """testing excluding gaps from model"""
        model = substitution_model.Nucleotide(model_gaps=False)
        assert len(model.getAlphabet()) == 4
                
    def test_includinggaps(self):
        """testing excluding gaps from model"""
        model = substitution_model.Nucleotide(model_gaps=True)
        assert len(model.getAlphabet()) == 5
                
    def test_getMotifs(self):
        """testing return of motifs"""
        model_motifs = substitution_model.Nucleotide().getMotifs()
    
    def test_getParamList(self):
        """testing getting the parameter list"""
        model = substitution_model.Nucleotide()
        self.assertEqual(model.getParamList(), [])

        model = substitution_model.Nucleotide(predicates=['beta:transition'])
        self.assertEqual(model.getParamList(), ['beta'])
        
    
    # need to ensure entering motif probs that sum to 1, that motif sets are the same

if __name__ == '__main__':
        main()


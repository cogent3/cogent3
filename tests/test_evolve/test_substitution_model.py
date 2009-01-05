#!/usr/bin/env python

import unittest
import os

from cogent import LoadSeqs, CodonAlphabet, DNA
from cogent.core import genetic_code
from cogent.evolve import substitution_model

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

class NucleotideModelTestMethods(unittest.TestCase):
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
        
class MultiLetterMotifSubstModelTests(unittest.TestCase):
    def setUp(self):
        self.submodel = substitution_model.Dinucleotide(do_scaling=True, 
                model_gaps=True)
        
    def test_asciiArt(self):
        model = substitution_model.Dinucleotide(predicates=['k:transition'])
        model.asciiArt()
        model = substitution_model.Dinucleotide()
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
    
            
class ThreeLetterMotifSubstModelTests(unittest.TestCase):
    def setUp(self):
        self.submodel = substitution_model.Nucleotide(motif_length=3)
        
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

class CodonSubstModelTests(unittest.TestCase):
    def setUp(self):
        self.standardcode = substitution_model.Codon(model_gaps=True, gc=1)
        self.mitocode = substitution_model.Codon(model_gaps=False, gc=2)
        
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

    
class ModelDataInteractionTestMethods(unittest.TestCase):
        
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
        unittest.main()


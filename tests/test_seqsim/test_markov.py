#!/usr/bin/env/python
"""test_markov.py: tests of the MarkovGenerator class.

"""
from cogent.seqsim.markov import MarkovGenerator
from StringIO import StringIO
from operator import mul
from sys import path
from cogent.util.unit_test import TestCase, main
from numpy import array

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Jesse Zaneveld", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class MarkovGeneratorTests(TestCase):
    """Tests of the MarkovGenerator class."""
    def setUp(self):
        """Define a few well-known frequencies."""
        self.single = MarkovGenerator(['UUUUUUUUUU'], order=0)
        self.equal = MarkovGenerator(['UUUUUCCCCC'], order=0)
        self.unequal_0 = MarkovGenerator(['UCCCCCCCCC'], order=0)
        self.unequal_1 = MarkovGenerator(['UCCCCCCCCC'], order=-1)
        self.pairs = MarkovGenerator(['UCAUCAUCAUCAUCA'], order=1)
        self.randquads=MarkovGenerator(['AACCUAUCUACUACUAUCUUCAUAUUCC']\
            ,order=3, calc_entropy=True, delete_bad_suffixes=False)
        self.empty = MarkovGenerator('', order=0)
        self.linebreaks= MarkovGenerator(StringIO('abb\nbcc\nd\n'))
        self.dinucs=MarkovGenerator(['ATACATAC'],order=1)
        self.orderfive=MarkovGenerator(['AAAAAGAAAAATAAAAAGAAAAAT'],order=5)
        
    def test_init(self):
        """MarkovGenerator init should give right frequency distributions."""
        self.assertEqual(self.empty.Frequencies, {})
        self.assertEqual(self.single.Frequencies, {'':{'U':1.0}})
        self.assertEqual(self.equal.Frequencies, {'':{'U':0.5,'C':0.5}})
        self.assertEqual(self.unequal_0.Frequencies, {'':{'U':0.1,'C':0.9}})
        self.assertEqual(self.unequal_1.Frequencies, {'':{'U':0.5, 'C':0.5}})
        self.assertEqual(self.pairs.Frequencies, \
            {'U':{'C':1},'C':{'A':1},'A':{'U':1}})
        #check that recalculating the frequencies doesn't break anything
        self.pairs.calcFrequencies()
        self.assertEqual(self.pairs.Frequencies, \
            {'U':{'C':1},'C':{'A':1},'A':{'U':1}})
        exp={'AAC':{'C':1},'ACC':{'U':1},'CCU':{'A':1},'CUA':{'U':0.5,'C':0.5},\
             'UAU':{'U':1/3.0,'C':2/3.0},'AUC':{'U':1},'UCU':{'U':0.5,'A':0.5},\
             'UAC':{'U':1},'ACU':{'A':1},'CUU':{'C':1},'UUC':{'C':0.5,'A':0.5},\
             'UCA':{'U':1},'CAU':{'A':1},'AUA':{'U':1},'AUU':{'C':1},
             }
        obs = self.randquads.Frequencies
        self.assertFloatEqual(obs, exp)
        #check that resetting linebreaks has the desired effect
        self.assertEqual(self.linebreaks.Frequencies, \
            {'a':{'b':1},'b':{'b':0.5,'c':0.5},'c':{'c':1}})
        self.linebreaks.Linebreaks = True
        self.linebreaks.Text.seek(0)
        self.linebreaks.calcFrequencies()
        #NOTE: current algorithm won't extend over line breaks. If you want
        #to force use of line breaks, read into a single string.
        self.assertEqual(self.linebreaks.Frequencies, \
            {'a':{'b':1},'b':{'b':0.5,'c':0.5},'c':{'c':1}})
    
    def test_next(self):
        """MarkovGenerator.next should generate text with expected properties"""
        #haven't figured how to do this for longer correlation lengths yet
        pass

    def test_entropy(self):
        """MarkovGenerator._entropy() should correctly calculate average H"""
        self.assertFloatEqual(self.randquads.Entropy, \
        3.0/25 * 0.91829583405448956 + 8.0/25)
    
    def test_evaluateProbability(self):
        """Should calculate  proper P value for seq"""
        self.dinucs.Prior=1
        q=self.dinucs.evaluateProbability('AT')
        self.assertFloatEqual(q,.50)
        z=self.dinucs.evaluateProbability('ATAT')
        self.assertFloatEqual(z,.25)
        p=self.dinucs.evaluateProbability('ATATAT')
        self.assertFloatEqual(p,.125)
        j=self.dinucs.evaluateProbability('ATACAT')
        self.assertFloatEqual(j,.125)
        h=self.orderfive.evaluateProbability('AAAAAT')
        self.assertFloatEqual(h,.50)
    
    def test_replaceDegenerateBases(self):
        """strips degenerate bases...."""
        text = 'AATCGCRRCCYAATC'
        m=MarkovGenerator([text],order=2)
        self.assertEqual(m.Text, [text])
        m.replaceDegenerateBases()
        self.assertEqual(m.Text[0][0:6],'aatcgc')
        p=m.Text[0][6]
        q= p in ['a','t','c','g']
        self.assertEqual(q,True)

    def test_wordToUniqueKey(self):
        """wordToUniqueKey should generate proper integers"""
        m=MarkovGenerator(['aataacaataac'],order=2)
        word='gca'
        uniqueKey=m.wordToUniqueKey(word)
        #a=0 c=1 t=2 g=3
        #should be (4^0*3)+(4^1*1)+(4^2*0)=3+4+0=7
        self.assertEqual(uniqueKey,7)
   
    def test_evaluateArrayProbability(self):
        """evaluateArrayProbability should calc prob from array indices"""
        m=MarkovGenerator(['aaaaaaaatt'],order=0)
        #8 a's, 2 t's
        m.calcFrequencies()
        prob=m.evaluateArrayProbability(array([0,2]))
        self.assertFloatEqual(prob,0.16) #0.8*0.2
        prob=m.evaluateArrayProbability(array([0,1]))
        self.assertFloatEqual(prob,0) #0.8*0



if __name__ == '__main__':
    main()

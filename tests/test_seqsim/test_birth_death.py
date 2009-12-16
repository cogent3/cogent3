#!/usr/bin/env python
#file test_birth_death.py

"""Unit tests of birth_death.py: implementation of the birth-death model.
"""
from cogent.seqsim.birth_death import ExtinctionError, TooManyTaxaError, \
    BirthDeathModel, DoubleBirthDeathModel
from cogent.util.unit_test import TestCase, main, FakeRandom

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Mike Robeson"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class BirthDeathModelTests(TestCase):
    """Tests of the BirthDeathModel class, which makes birth-death trees."""
    
    def test_init_deafults(self):
        """BirthDeathModel should init correctly w/ default params"""
        m = BirthDeathModel(0.1, 0.2, 0.3)
        self.assertEqual(m.BirthProb, 0.1)
        self.assertEqual(m.DeathProb, 0.2)
        self.assertEqual(m.TimePerStep, 0.3)
        self.assertEqual(m.MaxStep, 1000)
        self.assertEqual(m.MaxTaxa, None)
        self.assertEqual(m.CurrStep, 0)
        self.assertEqual(m.Tree.__class__, m.NodeClass)
        self.assertEqual(m.CurrTaxa, [m.Tree])
        self.assertEqual(m.ChangedBirthProb,None)
        self.assertEqual(m.ChangedDeathProb,None)
        self.assertEqual(m.ChangedBirthStep,None)
        self.assertEqual(m.ChangedDeathStep,None)
        self.assertEqual(m.CurrBirthProb, 0.1)
        self.assertEqual(m.CurrDeathProb, 0.2)

    def test_init_bad(self):
        """BirthDeathModel should raise exceptions on init with bad data"""
        #BirthProb and DeathProb must be probabilities between 0 and 1
        self.assertRaises(ValueError, BirthDeathModel, -1, 0.2, 0.3)
        self.assertRaises(ValueError, BirthDeathModel, 2, 0.2, 0.3)
        self.assertRaises(ValueError, BirthDeathModel, 0.1, -1, 0.3)
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 2, 0.3)
        #TimePerStep can't be negative or 0
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, -1)
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, 0)

    def test_init_extras(self):
        """BirthDeathModel should init OK with extra params"""
        m = BirthDeathModel(BirthProb=0.1, DeathProb=0.2, TimePerStep=0.3, \
            ChangedBirthProb=0.4, ChangedDeathProb=0.3, ChangedBirthStep=3,\
            ChangedDeathStep=4, MaxStep=5, MaxTaxa=10)
        self.assertEqual(m.BirthProb, 0.1)
        self.assertEqual(m.DeathProb, 0.2)
        self.assertEqual(m.TimePerStep, 0.3)
        self.assertEqual(m.ChangedBirthProb, 0.4)
        self.assertEqual(m.ChangedDeathProb, 0.3)
        self.assertEqual(m.ChangedBirthStep, 3)
        self.assertEqual(m.ChangedDeathStep, 4)
        self.assertEqual(m.MaxStep, 5)
        self.assertEqual(m.MaxTaxa, 10)
        self.assertEqual(m.CurrStep, 0)
        self.assertEqual(m.Tree.__class__, m.NodeClass)
        self.assertEqual(m.CurrTaxa, [m.Tree])

    def test_step(self):
        """BirthDeathModel step should match hand-calculated results"""
        m = BirthDeathModel(BirthProb=0.1, DeathProb=0.2, TimePerStep=1)
        born_and_died = FakeRandom([0],True)
        born_only = FakeRandom([1,0],True)
        died_only = FakeRandom([0,1],True)
        neither = FakeRandom([1],True)
        kill_alternate = FakeRandom([0,1,1,1], True)
        born_alternate = FakeRandom([1,1,1,0], True)
        #check that with neither birth nor death, we just continue
        m.step(neither)
        self.assertEqual(len(m.Tree.Children), 0)
        #check that with born_only we get a duplication
        m.step(born_only)
        self.assertEqual(len(m.Tree.Children), 2)
        assert m.Tree not in m.CurrTaxa
        for i in m.CurrTaxa:
            assert i.Parent is m.Tree
            self.assertEqual(i.Length, 1)
        #check that with a second round of born_only we duplicate again
        m.step(born_only)
        self.assertEqual(len(m.Tree.Children), 2)
        self.assertEqual(len(list(m.Tree.traverse())), 4)
        for i in m.Tree.traverse():
            self.assertEqual(i.Length, 1)
        for i in m.Tree.Children:
            self.assertEqual(i.Length, 1)
        #check that branch lengths add correctly
        for i in range(4):
            m.step(neither)
        self.assertEqual(len(m.CurrTaxa), 4)
        self.assertEqual(len(m.Tree.Children), 2)
        self.assertEqual(len(list(m.Tree.traverse())), 4)
        for i in m.Tree.traverse():
            self.assertEqual(i.Length, 5)
        for i in m.Tree.Children:
            self.assertEqual(i.Length, 1)
        #check that we can kill offspring correctly
        m.step(kill_alternate)
        self.assertEqual(len(m.CurrTaxa), 2)
        #make sure we killed the right children
        m.Tree.assignIds()
        for i in m.Tree.Children:
            #note that killing a child doesn't remove it, just stops it changing
            self.assertEqual(len(i.Children), 2)
            self.assertEqual(i.Children[0].Length, 5)
            self.assertEqual(i.Children[1].Length, 6)
        self.assertEqual([i.Length for i in m.Tree.traverse()], \
            [5,6,5,6])
        #make sure that born_and_died does the same thing as neither
        m.step(born_and_died)
        self.assertEqual([i.Length for i in m.Tree.traverse()], \
            [5,7,5,7])
        m.step(neither)
        self.assertEqual([i.Length for i in m.Tree.traverse()], \
            [5,8,5,8])
        #check that only CurrTaxa are brought forward
        self.assertEqual([i.Length for i in m.CurrTaxa], [8,8])
        #check that we can duplicate a particular taxon
        m.step(born_alternate)
        self.assertEqual([i.Length for i in m.CurrTaxa], [9,1,1])
        self.assertEqual(m.CurrTaxa[1].Parent.Length, 8)
        #check that we can kill 'em all
        m.step(died_only)
        self.assertEqual(len(m.CurrTaxa), 0)

    def test_prob_step_check(self):
        """prob_check and step_check should return error when out of bounds.
        Prob values should be between zero and one
        Step values should be greater than zero
        """
        #ChangedBirthProb = -0.1 , raises ValueError
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, 0.3,\
        ChangedBirthProb=-0.1,ChangedBirthStep=3,ChangedDeathProb=0.3,\
        ChangedDeathStep=4, MaxStep=5)
        #ChangedBirthStep = 0 , raises ValueError
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, 0.3,\
        ChangedBirthProb=0.6,ChangedBirthStep=0,ChangedDeathProb=0.3,\
        ChangedDeathStep=4, MaxStep=5)
        #ChangedDeathProb = 2 , raises ValueError
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, 0.3,\
        ChangedBirthProb=0.6,ChangedBirthStep=3,ChangedDeathProb=2,\
        ChangedDeathStep=4, MaxStep=5)
        #ChangedDeathStep = -1 , raises ValueError
        self.assertRaises(ValueError, BirthDeathModel, 0.1, 0.2, 0.3,\
        ChangedBirthProb=0.6,ChangedBirthStep=3,ChangedDeathProb=0.3,\
        ChangedDeathStep=-1, MaxStep=5)
 

    def test_timeOk(self):
        """BirthDeathModel TimeOk should return True if time not exceeded"""
        b = BirthDeathModel(0.1, 0.2, 0.3, MaxStep=5)
        assert b.timeOk()
        b.CurrStep = 4
        assert b.timeOk()
        b.CurrStep = 5
        assert not b.timeOk()
        b.CurrStep = 1000
        assert not b.timeOk()
        b.MaxStep = None
        assert b.timeOk()
        b.MaxStep = 1001
        assert b.timeOk()
        b.step()
        assert not b.timeOk()

    def test_taxaOk(self):
        """BirthDeathModel TaxaOk should return True if taxa not exceeded"""
        b = BirthDeathModel(0.1, 0.2, 0.3, MaxTaxa=5)
        born_alternate = FakeRandom([1,1,1,0], True)
        born_only = FakeRandom([1,0],True)
        kill_only = FakeRandom([0,1,0,1], True)
        #start off with single taxon
        assert b.taxaOk()
        #taxa are OK if there are a few
        b.step(born_only)   #now 2 taxa
        assert b.taxaOk()   
        b.step(born_only)   #now 4 taxa
        assert b.taxaOk()
        b.step(born_only)   #now 8 taxa
        assert not b.taxaOk()
        b.MaxTaxa = 8
        assert not b.taxaOk()
        b.MaxTaxa = 9
        assert b.taxaOk()
        b.MaxTaxa = 17
        assert b.taxaOk()
        b.step(born_only)
        assert b.taxaOk()
        b.step(born_only)
        assert not b.taxaOk()
        #ok if no maximum
        b.MaxTaxa = None
        assert b.taxaOk()
        #not ok if there are no taxa left
        b.step(kill_only)
        assert not b.taxaOk()
        #still not OK if not MaxTaxa
        b.MaxTaxa = None
        assert not b.taxaOk()
        
    def test_call_exact(self):
        """BirthDeathModel call should produce right # taxa when exact"""
        m = BirthDeathModel(0.01, 0.005, 0.1, MaxTaxa=10)
        for i in range(10):
            try:
                result = m(filter=True, exact=True)
                self.assertEqual(len(list(result.traverse())), 10)
            except (TooManyTaxaError, ExtinctionError), e:
                pass
    
    def test_call(self):
        """BirthDeathModel call should produce hand-calculated trees"""
        m = BirthDeathModel(0.01, 0.005, 0.1, MaxTaxa=10)
        r = FakeRandom(\
        [1,0,\
         1,1, 1,1,\
         1,0, 0,0,\
         0,0, 0,0, 1,0,\
         0,0, 0,0, 0,1, 0,0, \
         1,0, 0,0, 0,0,\
         1,0, 0,0, 0,0, 1,0, \
         1,0, 1,0, 0,1, 1,1, 1,0, 1,0, \
         1,1, 1,1, 1,1, 1,1, 1,0, 1,1, 1,1, 1,1, 1,1], True)
        m = BirthDeathModel(0.1, 0.5, 1, MaxTaxa=10)
        result = m(filter=False, random_f=r)
        self.assertEqual([i.Length for i in result.traverse()], \
            [2,2,2,2,2,1,1,1,2,2,2,2])
        #try it with pruning
        m = BirthDeathModel(0.1, 0.5, 1, MaxTaxa=10)
        result = m(filter=True, random_f=r)
        self.assertEqual([i.Length for i in result.traverse()], \
            [2,2,2,2,1,1,2,2,2,2])
       #try it with fewer taxa
        m = BirthDeathModel(0.1, 0.5, 1, MaxTaxa=4)
        result = m(filter=True, random_f=r)
        self.assertEqual([i.Length for i in result.traverse()], \
            [2,2,1,1])

    def test_changed_values_step(self):
        """Tests if values changed at specified steps in step().
        Note, in m.step() CurrStep is logically tested one step later. 
        """
        m = BirthDeathModel( 0.1, 0.2, 0.3,ChangedBirthProb=0.6,\
            ChangedBirthStep=3,ChangedDeathProb=0.3,ChangedDeathStep=4,\
            MaxStep=5)
        # all values should be as initialized
        m.step()
        assert m.CurrStep == 1
        assert m.BirthProb == 0.1
        assert m.DeathProb == 0.2
        assert m.CurrBirthProb == 0.1
        assert m.CurrDeathProb == 0.2
        # continue 2 steps
        m.step()
        m.step()
        # when logically evaluated CurrBirthProb should change
        # from 0.1 to 0.6
        m.step()
        assert m.CurrStep == 4
        assert m.BirthProb == 0.1
        assert m.DeathProb == 0.2
        assert m.CurrBirthProb == 0.6
        assert m.CurrDeathProb == 0.2
        # All values other than CurrStep should be as above
        # except that CurrDeathProb should change from 0.2 to 0.3
        m.step()
        assert m.CurrStep == 5
        assert m.BirthProb == 0.1
        assert m.DeathProb == 0.2
        assert m.CurrBirthProb == 0.6
        assert m.CurrDeathProb == 0.3


class DoubleBirthDeathTests(TestCase):
    """Tests of the double birth-death model."""
    def test_double_birth_death(self):
        """double_birth_death should run without errors"""
        pass 


if __name__ == "__main__":
    main()

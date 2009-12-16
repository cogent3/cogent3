#!/usr/bin/env python

"""Tests Filter, Organizer and filterfunctions, classes for filtering
"""

from cogent.util.organizer import Filter, Organizer, GroupList, regroup
from cogent.util.transform import find_any, find_no, find_all,\
    keep_if_more, exclude_if_more, keep_if_more_other, exclude_if_more_other
from cogent.util.unit_test import TestCase,main

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class Sequence(object):
    """Simple sequence class for tests."""
    def __init__(self, s, info):
        self.s = s
        self.__dict__.update(info)
        if not hasattr(self, 'Gene'):
            self.Gene = None
    def __contains__(self, s):
        return s in self.s
    def __repr__(self):
        return `self.s`
    def __iter__(self):
        return iter(self.s)
    def __nonzero__(self):
        return bool(self.s)
    def lower(self):
        return self.s.lower()
    def __cmp__(self, other):
        return cmp(self.s,other.s)

class FilterTests(TestCase):
    """Tests of Filter class"""

    def test_init(self):
        """Filter should init as expected"""
        empty_filter = Filter('',{})
        named_empty_filter = Filter('Archaea',{})

        self.assertEqual(empty_filter,{})
        self.assertEqual(empty_filter.Name,'')
        self.assertEqual(named_empty_filter,{})
        self.assertEqual(named_empty_filter.Name,'Archaea')
        
        f = find_all('abcd')
        g = keep_if_more_other('ab',7)
        fil = Filter('Archaea',{'Arch':[f,g]})
        assert fil['Arch'][0] is f
        assert fil['Arch'][1] is g

    def test_call_empty(self):
        """Empty Filter should return True when called on anything"""
        f = Filter('',{})
        data = ['aa','bb','cc']
        self.assertEqual(f(data),True)
    
    def test_call_full(self):
        """Filter should return True if the object satisfies all criteria"""
        seq1 = Sequence('ACGU',{'Gene':'LSU'})
        seq2 = Sequence('ACGUACGU',{'Gene':'SSU'})
        seq3 = Sequence('ACGUN',{'Gene':'LSU'})
        seq4 = Sequence('ACG',{'Gene':'LSU'})
        seq5 = Sequence('ACGU',{})
        seq6 = Sequence('',{})
        
        f = Filter('valid',{None:[find_all('AGCU'),find_no('N')],\
                            'Gene':[find_any(['LSU'])]})
        self.assertEqual(f(seq1),True)
        self.assertEqual(f(seq2),False)
        self.assertEqual(f(seq3),False)
        self.assertEqual(f(seq4),False)
        self.assertEqual(f(seq5),False)
        self.assertEqual(f(seq6),False)
        
class GroupListTests(TestCase):
    """Tests of GroupList class"""

    def test_init_empty(self):
        """Empty GroupList should init OK"""
        g = GroupList([])
        self.assertEqual(len(g),0)
        self.assertEqual(g.Groups,[])

    def test_init_full(self):
        """GroupList should init OK with data and groups"""
        data = ['a','b','c']
        groups = [1,2,3]
        g = GroupList(data,groups)
        self.assertEqual(g,data)
        self.assertEqual(g.Groups,groups)
        self.assertEqual(len(g),3)

        
class OrganizerTests(TestCase):
    """Tests of Classifier class"""

    def setUp(self):
        """Define some standard Organizers for testing"""
        self.Empty = Organizer([])
        
        self.a = Filter('a',{None:[find_any('a')]})
        self.b = Filter('b',{None:[find_any('b')]})
        self.Ab_org = Organizer([self.a,self.b])
        
        lsu = Filter('LSU',{None:[exclude_if_more('N',5)],\
                            'Gene':[find_any(['LSU'])]})
        ssu = Filter('SSU',{None:[exclude_if_more('N',5)],\
                            'Gene':[find_any(['SSU'])]})
        self.Gene_org = Organizer([lsu,ssu])
        
        self.Ab_seq = ['aa','bb','abab','cc','']
        self.seq1 = Sequence('ACGU',{'Gene':'LSU'})
        self.seq2 = Sequence('ACGUACGU',{'Gene':'SSU'})
        self.seq3 = Sequence('ACGUNNNNNN',{'Gene':'LSU'})
        self.seq4 = Sequence('ACGUNNNNNN',{'Gene':'SSU'})
        self.seq5 = Sequence('ACGU',{})
        self.seq6 = Sequence('',{})
        self.seq7 = Sequence('ACGU',{'Gene':'unit'})
        
        self.Gene_seq = [self.seq1,self.seq2,self.seq3,self.seq4,\
                self.seq5,self.seq6,self.seq7]

        f = Filter('valid',{None:[find_all('AGCU'),find_no('N')],\
                            'Gene':[find_any(['LSU'])]})
        self.Mult_func_org = Organizer([f])
    
    def test_init_empty(self):
        """Empty Organizer should init correctly"""
        org = self.Empty
        self.assertEqual(len(org),0)

    def test_init_full(self):
        """Organizer should init correctly with multiple functions"""
        org = Organizer([self.a,self.b])
        self.assertEqual(org[0],self.a)
        self.assertEqual(org[1],self.b)
        self.assertEqual(len(org),2)

    def test_empty_org_empty_list(self):
        """Empty Organizer should return [] when applied to []"""
        org = self.Empty
        l = []
        self.assertEqual(org(l),[])

    def test_empty_org_full_list(self):
        """Empty organizer, applied to full list, should return the original"""
        org = self.Empty
        l = self.Ab_seq
        obs = org(l)
        self.assertEqual(obs,[l])
        self.assertEqual(obs[0].Groups,[None])

    def test_full_org_empty_list(self):
        """Organizer should return [] when applied to []"""
        org = self.Ab_org
        l = []
        obs = org(l)
        self.assertEqual(obs,[])

    def test_full_org_full_list(self):
        """Organizer should return correct organization"""
        org = self.Ab_org
        l = self.Ab_seq
        obs = org(l)
        obs.sort()
        exp = [['aa','abab'],['bb'],['cc','']]
        self.assertEqual(obs,exp)
        self.assertEqual(obs[0].Groups,['a'])
        self.assertEqual(obs[1].Groups,['b'])
        self.assertEqual(obs[2].Groups,[None])

    def test_double_org_empty_list(self):
        """Organizer should return [] when applied to []"""
        org = self.Gene_org
        l = []
        obs = org(l)
        self.assertEqual(obs,[])

    def test_double_org_full_list(self):
        """Organizer should handle multiple filters correctly"""
        org = self.Gene_org
        l = self.Gene_seq
        obs = org(l)
        obs.sort()
        exp = [[self.seq1],[self.seq2],[self.seq3,self.seq4,\
                    self.seq5,self.seq6,self.seq7]]
        self.assertEqual(obs,exp)
        self.assertEqual(obs[0].Groups,['LSU'])
        self.assertEqual(obs[1].Groups,['SSU'])
        self.assertEqual(obs[2].Groups,[None])

    def test_multiple_func(self):
        """Organizer should handle filter with multiple functions correctly"""
        org = self.Mult_func_org
        l = self.Gene_seq
        obs = org(l)
        obs.sort()
        exp = [[self.seq1],[self.seq2,self.seq3,self.seq4,self.seq5,\
                    self.seq6,self.seq7]]
        self.assertEqual(obs,exp)
        self.assertEqual(obs[0].Groups,['valid'])
        self.assertEqual(obs[1].Groups,[None])


class organizerTests(TestCase):
    """Tests for module-level functions"""

    def test_regroup(self):
        """regroup: should groups with identical hierarchy-info together"""
        g1 = GroupList([1,2,3],['a'])
        g2 = GroupList([4,5,6],['b'])
        g3 = GroupList([7,7,7],['a','b'])
        g4 = GroupList([8,8,8],['a'])
        all = [g1, g2, g3, g4]
        self.assertEqualItems(regroup(all), [[1,2,3,8,8,8],[7,7,7],[4,5,6]])


        
if __name__ == "__main__":
    main()   



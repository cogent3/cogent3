#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.core.sequence import RnaSequence, frac_same, ModelSequence, Sequence
from cogent.maths.stats.util import Freqs, Numbers
from cogent.core.moltype import RNA, DNA, PROTEIN, BYTES
from cogent.struct.rna2d import ViennaStructure

from cogent.core.alignment import SequenceCollection, \
    make_gap_filter, coerce_to_string, \
    seqs_from_array, seqs_from_model_seqs, seqs_from_generic, seqs_from_fasta, \
    seqs_from_dict, seqs_from_aln, seqs_from_kv_pairs, seqs_from_empty, \
    aln_from_array, aln_from_model_seqs, aln_from_collection,\
    aln_from_generic, aln_from_fasta, aln_from_dense_aln, aln_from_empty, \
    DenseAlignment, Alignment, DataError

from cogent.core.moltype import AB, DNA
from cogent.parse.fasta import MinimalFastaParser
from numpy import array, arange, transpose
from tempfile import mktemp
from os import remove
import re

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Catherine Lozuopone", "Gavin Huttley",
                    "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class alignment_tests(TestCase):
    """Tests of top-level functions."""

    def test_seqs_from_array(self):
        """seqs_from_array should return chars, and successive indices."""
        a = array([[0,1,2],[2,1,0]])    #three 2-char seqs
        obs_a, obs_labels = seqs_from_array(a)
        #note transposition
        self.assertEqual(obs_a, [array([0,2]), array([1,1]), array([2,0])])
        self.assertEqual(obs_labels, None)

    def test_seqs_from_model_seqs(self):
        """seqs_from_model_seqs should return model seqs + names."""
        s1 = ModelSequence('ABC', Name='a')
        s2 = ModelSequence('DEF', Name='b')
        obs_a, obs_labels = seqs_from_model_seqs([s1, s2])
        self.assertEqual(obs_a, [s1,s2]) #seq -> numbers
        self.assertEqual(obs_labels, ['a','b'])

    def test_seqs_from_generic(self):
        """seqs_from_generic should initialize seqs from list of lists, etc."""
        s1 = 'ABC'
        s2 = 'DEF'
        obs_a, obs_labels = seqs_from_generic([s1, s2])
        self.assertEqual(obs_a, ['ABC','DEF'])
        self.assertEqual(obs_labels, [None, None])
    
    def test_seqs_from_fasta(self):
        """seqs_from_fasta should initialize seqs from fasta-format string"""
        s = '>aa\nAB\nC\n>bb\nDE\nF\n'
        obs_a, obs_labels = seqs_from_fasta(s)
        self.assertEqual(obs_a, ['ABC','DEF'])
        self.assertEqual(obs_labels, ['aa','bb'])

    def test_seqs_from_aln(self):
        """seqs_from_aln should initialize from existing alignment"""
        c = SequenceCollection(['abc','def'])
        obs_a, obs_labels = seqs_from_aln(c)
        self.assertEqual(obs_a, ['abc','def'])
        self.assertEqual(obs_labels, ['seq_0','seq_1'])

    def test_seqs_from_kv_pairs(self):
        """seqs_from_kv_pairs should initialize from key-value pairs"""
        c = [['a', 'abc'], ['b', 'def']]
        obs_a, obs_labels = seqs_from_kv_pairs(c)
        self.assertEqual(obs_a, ['abc','def'])
        self.assertEqual(obs_labels, ['a','b'])

    def test_seqs_from_empty(self):
        """seqs_from_empty should always raise ValueError"""
        self.assertRaises(ValueError, seqs_from_empty, 'xyz')
    
    def test_aln_from_array(self):
        """aln_from_array should return same array, and successive indices."""
        a = array([[0,1,2],[3,4,5]])    #three 2-char seqs
        obs_a, obs_labels = aln_from_array(a)
        self.assertEqual(obs_a, transpose(a))
        self.assertEqual(obs_labels, None)

    def test_aln_from_model_seqs(self):
        """aln_from_model_seqs should initialize aln from sequence objects."""
        s1 = ModelSequence('ACC', Name='a', Alphabet=RNA.Alphabet)
        s2 = ModelSequence('GGU', Name='b', Alphabet=RNA.Alphabet)
        obs_a, obs_labels = aln_from_model_seqs([s1, s2], \
            Alphabet=BYTES.Alphabet)
        self.assertEqual(obs_a, array([[2,1,1],[3,3,0]], 'b'))
        #seq -> numbers
        self.assertEqual(obs_labels, ['a','b'])

    def test_aln_from_generic(self):
        """aln_from_generic should initialize aln from list of lists, etc."""
        s1 = 'AAA'
        s2 = 'GGG'
        obs_a, obs_labels = aln_from_generic([s1, s2], 'b', \
            Alphabet=RNA.Alphabet) #specify array type
        self.assertEqual(obs_a, array([[2,2,2],[3,3,3]], 'b')) #str -> chars
        self.assertEqual(obs_labels, [None, None])
    
    def test_aln_from_fasta(self):
        """aln_from_fasta should initialize aln from fasta-format string"""
        s = '>aa\nAB\nC\n>bb\nDE\nF\n'
        obs_a, obs_labels = aln_from_fasta(s.splitlines())
        self.assertEqual(obs_a, array(['ABC','DEF'], 'c').view('B')) #seq -> numbers
        self.assertEqual(obs_labels, ['aa','bb'])

    def test_aln_from_dense_aln(self):
        """aln_from_dense_aln should initialize from existing alignment"""
        a = DenseAlignment(array([[0,1,2],[3,4,5]]), conversion_f=aln_from_array)
        obs_a, obs_labels = aln_from_dense_aln(a)
        self.assertEqual(obs_a, a.SeqData)
        self.assertEqual(obs_labels, a.Names)

    def test_aln_from_collection(self):
        """aln_from_collection should initialize from existing alignment"""
        a = SequenceCollection(['AAA','GGG'])
        obs_a, obs_labels = aln_from_collection(a, Alphabet=RNA.Alphabet)
        self.assertEqual(a.toFasta(), '>seq_0\nAAA\n>seq_1\nGGG')
        self.assertEqual(obs_a, array([[2,2,2],[3,3,3]]))
    
    def test_aln_from_empty(self):
        """aln_from_empty should always raise ValueError"""
        self.assertRaises(ValueError, aln_from_empty, 'xyz')


class SequenceCollectionBaseTests(object):
    """Base class for testing the SequenceCollection object.
    
    Unlike Alignments, SequenceCollections can have sequences that are not equal
    length. This module contains all the code that _doesn't_ depend on being
    able to look at "ragged" SequenceCollections. It is intended that all
    classes that inherit from SequenceCollection should have test classes that
    inherit from this class, but that the SequenceCollection tests themselves
    will additionally contain code to deal with SequenceCollections of unequal
    length.

    set self.Class in subclasses to generate the rught constructor.
    """
    Class = SequenceCollection

    def setUp(self):
        """Define some standard SequenceCollection objects."""
        self.one_seq = self.Class({'a':'AAAAA'})
        self.ragged_padded = self.Class({'a':'AAAAAA','b':'AAA---', \
            'c':'AAAA--'})
        self.identical = self.Class({'a':'AAAA','b':'AAAA'})
        self.gaps = self.Class({'a':'AAAAAAA','b':'A--A-AA', \
            'c':'AA-----'})
        self.gaps_rna = self.Class({'a':RnaSequence('AAAAAAA'), \
            'b':RnaSequence('A--A-AA'), \
            'c':RnaSequence('AA-----')})
        self.unordered = self.Class({'a':'AAAAA','b':'BBBBB'})
        self.ordered1 = self.Class({'a':'AAAAA','b':'BBBBB'}, \
            Names=['a','b'])
        self.ordered2 = self.Class({'a':'AAAAA','b':'BBBBB'}, \
            Names=['b','a'])
        self.mixed = self.Class({'a':'ABCDE', 'b':'LMNOP'})
        self.end_gaps = self.Class({'a':'--A-BC-', 'b':'-CB-A--', \
            'c':'--D-EF-'}, Names=['a','b','c'])
        self.many = self.Class({
            'a': RnaSequence('UCAGUCAGUU'),
            'b': RnaSequence('UCCGUCAAUU'),
            'c': RnaSequence('ACCAUCAGUC'),
            'd': RnaSequence('UCAAUCGGUU'),
            'e': RnaSequence('UUGGUUGGGU'),
            'f': RnaSequence('CCGGGCGGCC'),
            'g': RnaSequence('UCAACCGGAA'),
            })
        #Additional SequenceCollections for tests added 6/4/04 by Jeremy Widmann
        self.sequences = self.Class(map(RnaSequence, ['UCAG', 'UCAG', 'UCAG']))
        self.structures = self.Class(map(ViennaStructure, 
                                ['(())..', '......', '(....)']), MolType=BYTES)
        self.labeled = self.Class(['ABC', 'DEF'], ['1st', '2nd'])
        #Additional SequenceCollection for tests added 1/30/06 by Cathy Lozupone
        self.omitSeqsTemplate_aln = self.Class({
            's1':RnaSequence('UC-----CU---C'),
            's2':RnaSequence('UC------U---C'),
            's3':RnaSequence('UUCCUUCUU-UUC'),
            's4':RnaSequence('UU-UUUU-UUUUC'),
            's5':RnaSequence('-------------')
            })

        self.a = DenseAlignment(['AAA','AAA'])
        self.b = Alignment(['AAA','AAA'])
        self.c = SequenceCollection(['AAA','AAA'])

    def test_guess_input_type(self):
        """SequenceCollection  _guess_input_type should figure out data type correctly"""
        git = self.a._guess_input_type
        self.assertEqual(git(self.a), 'dense_aln')
        self.assertEqual(git(self.b), 'aln')
        self.assertEqual(git(self.c), 'collection')
        self.assertEqual(git('>ab\nabc'), 'fasta')
        self.assertEqual(git(['>ab','abc']), 'fasta')
        self.assertEqual(git(['abc','def']), 'generic')
        self.assertEqual(git([[1,2],[4,5]]), 'kv_pairs') #precedence over generic
        self.assertEqual(git([[1,2,3],[4,5,6]]), 'generic')
        self.assertEqual(git([ModelSequence('abc')]), 'model_seqs')
        self.assertEqual(git(array([[1,2,3],[4,5,6]])), 'array')
        self.assertEqual(git({'a':'aca'}), 'dict')
        self.assertEqual(git([]), 'empty')

    def test_init_aln(self):
        """ SequenceCollection should init from existing alignments"""
        exp = self.Class(['AAA','AAA'])
        x = self.Class(self.a)
        y = self.Class(self.b)
        z = self.Class(self.c)
        self.assertEqual(x, exp)
        self.assertEqual(z, exp)
        self.assertEqual(y, exp)
    test_init_aln.__doc__ = Class.__name__ + test_init_aln.__doc__
 
    def test_init_dict(self):
        """SequenceCollection init from dict should work as expected"""
        d = {'a':'AAAAA', 'b':'BBBBB'}
        a = self.Class(d)
        self.assertEqual(a, d)
        self.assertEqual(a.NamedSeqs.items(), d.items())

    def test_init_name_mapped(self):
        """SequenceCollection init should allow name mapping function"""
        d = {'a':'AAAAA', 'b':'BBBBB'}
        f = lambda x: x.upper()
        a = self.Class(d, label_to_name=f)
        self.assertNotEqual(a, d)
        self.assertNotEqual(a.NamedSeqs.items(), d.items())
        d_upper = {'A':'AAAAA','B':'BBBBB'}
        self.assertEqual(a, d_upper)
        self.assertEqual(a.NamedSeqs.items(), d_upper.items())

    def test_init_seq(self):
        """SequenceCollection init from list of sequences should use indices as keys"""
        seqs = ['AAAAA', 'BBBBB', 'CCCCC']
        a = self.Class(seqs)
        self.assertEqual(len(a.NamedSeqs), 3)
        self.assertEqual(a.NamedSeqs['seq_0'], 'AAAAA')
        self.assertEqual(a.NamedSeqs['seq_1'], 'BBBBB')
        self.assertEqual(a.NamedSeqs['seq_2'], 'CCCCC')
        self.assertEqual(a.Names, ['seq_0','seq_1','seq_2'])
        self.assertEqual(list(a.Seqs), ['AAAAA','BBBBB','CCCCC'])

    def test_init_annotated_seq(self):
        """SequenceCollection init from seqs w/ Info should preserve data"""
        a = Sequence('AAA', Name='a', Info={'x':3})
        b = Sequence('CCC', Name='b', Info={'x':4})
        c = Sequence('GGG', Name='c', Info={'x':5})
        seqs = [c,b,a]
        a = self.Class(seqs)
        self.assertEqual(list(a.Names), ['c','b','a'])
        self.assertEqual(map(str, a.Seqs), ['GGG','CCC','AAA'])
        if self.Class is not DenseAlignment:
            #DenseAlignment is allowed to strip Info objects
            self.assertEqual([i.Info.x for i in a.Seqs], [5,4,3])
        #check it still works if constructed from same class
        b = self.Class(a)
        self.assertEqual(list(b.Names), ['c','b','a'])
        self.assertEqual(map(str, b.Seqs), ['GGG','CCC','AAA'])
        if self.Class is not DenseAlignment:
            #DenseAlignment is allowed to strip Info objects
            self.assertEqual([i.Info.x for i in b.Seqs], [5,4,3])


    def test_init_pairs(self):
        """SequenceCollection init from list of (key,val) pairs should work correctly"""
        seqs = [['x', 'XXX'], ['b','BBB'], ['c','CCC']]
        a = self.Class(seqs)
        self.assertEqual(len(a.NamedSeqs), 3)
        self.assertEqual(a.NamedSeqs['x'], 'XXX')
        self.assertEqual(a.NamedSeqs['b'], 'BBB')
        self.assertEqual(a.NamedSeqs['c'], 'CCC')
        self.assertEqual(a.Names, ['x','b','c'])
        self.assertEqual(list(a.Seqs), ['XXX','BBB','CCC'])

    def test_init_duplicate_keys(self):
        """SequenceCollection init from (key, val) pairs should fail on dup. keys"""
        seqs = [['x', 'XXX'], ['b','BBB'],['x','CCC'], ['d','DDD'], ['a','AAA']]
        self.assertRaises(ValueError, self.Class, seqs)
        aln = self.Class(seqs, remove_duplicate_names=True)
        self.assertEqual(str(self.Class(seqs, remove_duplicate_names=True)),
            '>x\nXXX\n>b\nBBB\n>d\nDDD\n>a\nAAA\n')

    def test_init_ordered(self):
        """SequenceCollection should iterate over seqs correctly even if ordered"""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        self.assertEqual(first.Names, ['a','b'])
        self.assertEqual(sec.Names, ['b', 'a'])
        self.assertEqual(un.Names, un.NamedSeqs.keys())

        first_list = list(first.Seqs)
        sec_list = list(sec.Seqs)
        un_list = list(un.Seqs)

        self.assertEqual(first_list, ['AAAAA','BBBBB'])
        self.assertEqual(sec_list, ['BBBBB','AAAAA'])
    
        #check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)

    def test_init_ambig(self):
        """SequenceCollection should tolerate ambiguous chars"""
        aln = self.Class(['AAA','CCC'],MolType=DNA)
        aln = self.Class(['ANS','CWC'],MolType=DNA)
        aln = self.Class(['A-A','CC-'],MolType=DNA)
        aln = self.Class(['A?A','CC-'],MolType=DNA)

    def test_aln_from_fasta_parser(self):
        """aln_from_fasta_parser should init from iterator"""
        s = '>aa\nAC\n>bb\nAA\n>c\nGG\n'.splitlines()
        p = MinimalFastaParser(s)
        aln = self.Class(p, MolType=DNA)
        self.assertEqual(aln.NamedSeqs['aa'], 'AC')
        self.assertEqual(aln.toFasta(), '>aa\nAC\n>bb\nAA\n>c\nGG')
        s2_ORIG = '>x\nCA\n>b\nAA\n>>xx\nGG'
        s2 = '>aa\nAC\n>bb\nAA\n>c\nGG\n'
        d = DenseAlignment(MinimalFastaParser(s2.splitlines()))
        self.assertEqual(d.toFasta(), aln.toFasta())

    def test_aln_from_fasta(self):
        """SequenceCollection should init from fasta-format string"""
        s = '>aa\nAC\n>bb\nAA\n>c\nGG\n'
        aln = self.Class(s)
        self.assertEqual(aln.toFasta(), s.strip())

    def test_SeqLen_get(self):
        """SequenceCollection SeqLen should return length of longest seq"""
        self.assertEqual(self.one_seq.SeqLen, 5)
        self.assertEqual(self.identical.SeqLen, 4)
        self.assertEqual(self.gaps.SeqLen, 7)

    def test_Seqs(self):
        """SequenceCollection Seqs property should return seqs in correct order."""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        first_list = list(first.Seqs)
        sec_list = list(sec.Seqs)
        un_list = list(un.Seqs)

        self.assertEqual(first_list, ['AAAAA','BBBBB'])
        self.assertEqual(sec_list, ['BBBBB','AAAAA'])
    
        #check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)

        
    def test_iterSeqs(self):
        """SequenceCollection iterSeqs() method should support reordering of seqs"""
        self.ragged_padded = self.Class(self.ragged_padded.NamedSeqs, \
            Names=['a','b','c'])
        seqs = list(self.ragged_padded.iterSeqs())
        self.assertEqual(seqs, ['AAAAAA', 'AAA---', 'AAAA--'])
        seqs = list(self.ragged_padded.iterSeqs(seq_order=['b','a','a']))
        self.assertEqual(seqs, ['AAA---', 'AAAAAA', 'AAAAAA'])
        self.assertSameObj(seqs[1], seqs[2])
        self.assertSameObj(seqs[0], self.ragged_padded.NamedSeqs['b'])
        
    def test_Items(self):
        """SequenceCollection Items should iterate over items in specified order."""
        #should work if one row
        self.assertEqual(list(self.one_seq.Items), ['A']*5)
        #should take order into account
        self.assertEqual(list(self.ordered1.Items), ['A']*5 + ['B']*5)
        self.assertEqual(list(self.ordered2.Items), ['B']*5 + ['A']*5)
        
    def test_iterItems(self):
        """SequenceCollection iterItems() should iterate over items in correct order"""
        #should work if one row
        self.assertEqual(list(self.one_seq.iterItems()), ['A']*5)
        #should take order into account
        self.assertEqual(list(self.ordered1.iterItems()), ['A']*5 + ['B']*5)
        self.assertEqual(list(self.ordered2.iterItems()), ['B']*5 + ['A']*5)
        #should allow row and/or col specification
        r = self.ragged_padded
        self.assertEqual(list(r.iterItems(seq_order=['c','b'], \
            pos_order=[5,1,3])), list('-AA-A-'))
        #should not interfere with superclass iteritems()
        i = list(r.NamedSeqs.iteritems()) 
        i.sort()
        self.assertEqual(i, [('a','AAAAAA'),('b','AAA---'),('c','AAAA--')])
    
    def test_takeSeqs(self):
        """SequenceCollection takeSeqs should return new SequenceCollection with selected seqs."""
        a = self.ragged_padded.takeSeqs('bc')
        self.assertTrue(isinstance(a, SequenceCollection))
        self.assertEqual(a, {'b':'AAA---','c':'AAAA--'})
        #should be able to negate
        a = self.ragged_padded.takeSeqs('bc', negate=True)
        self.assertEqual(a, {'a':'AAAAAA'})

    def test_getSeqIndices(self):
        """SequenceCollection getSeqIndices should return names of seqs where f(row) is True"""
        srp = self.ragged_padded
        is_long = lambda x: len(x) > 10
        is_med = lambda x: len(str(x).replace('-','')) > 3 #strips gaps
        is_any = lambda x: len(x) > 0
        self.assertEqual(srp.getSeqIndices(is_long), [])
        srp.Names = 'cba'
        self.assertEqual(srp.getSeqIndices(is_med), ['c','a'])
        srp.Names = 'bac'
        self.assertEqual(srp.getSeqIndices(is_med), ['a','c'])
        self.assertEqual(srp.getSeqIndices(is_any),['b','a','c'])
        #should be able to negate
        self.assertEqual(srp.getSeqIndices(is_med, negate=True), ['b'])
        self.assertEqual(srp.getSeqIndices(is_any, negate=True), [])

    def test_takeSeqsIf(self):
        """SequenceCollection takeSeqsIf should return seqs where f(row) is True"""
        is_long = lambda x: len(x) > 10
        is_med = lambda x: len(str(x).replace('-','')) > 3
        is_any = lambda x: len(x) > 0
        srp = self.ragged_padded
        self.assertEqual(srp.takeSeqsIf(is_long), {})
        srp.Names = 'cba'
        self.assertEqual(srp.takeSeqsIf(is_med), \
            {'c':'AAAA--','a':'AAAAAA'})
        srp.Names = srp.NamedSeqs.keys()
        self.assertEqual(srp.takeSeqsIf(is_med), \
            {'c':'AAAA--','a':'AAAAAA'})
        self.assertEqual(srp.takeSeqsIf(is_any), srp)
        self.assertTrue(isinstance(srp.takeSeqsIf(is_med), SequenceCollection))
        #should be able to negate
        self.assertEqual(srp.takeSeqsIf(is_med, negate=True), \
            {'b':'AAA---'})

    def test_getItems(self):
        """SequenceCollection getItems should return list of items from k,v pairs"""
        self.assertEqual(self.mixed.getItems([('a',3),('b',4),('a',0)]), \
            ['D','P','A'])
        self.assertRaises(KeyError, self.mixed.getItems, [('x','y')])
        self.assertRaises(IndexError, self.mixed.getItems, [('a',1000)])
        #should be able to negate -- note that results will have seqs in
        #arbitrary order
        self.assertEqualItems(self.mixed.getItems([('a',3),('b',4),('a',0)], \
            negate=True), ['B','C','E','L','M','N','O'])

    def test_getItemIndices(self):
        """SequenceCollection getItemIndices should return coordinates of matching items"""
        is_vowel = lambda x: x in 'AEIOU'
        #reverse name order to test that it's not alphabetical
        self.mixed = self.Class(self.mixed.NamedSeqs, Names=['b','a'])    
        self.assertEqual(self.mixed.getItemIndices(is_vowel), \
            [('b',3),('a',0),('a',4)])
        is_lower = lambda x: x.islower()
        self.assertEqual(self.ragged_padded.getItemIndices(is_lower), [])
        #should be able to negate
        self.assertEqualItems(self.mixed.getItemIndices(is_vowel, negate=True),\
            [('a',1),('a',2),('a',3),('b',0),('b',1),('b',2),('b',4)])

    def test_getItemsIf(self):
        """SequenceCollection getItemsIf should return matching items"""
        is_vowel = lambda x: x in 'AEIOU'
        #reverse name order to test that it's not alphabetical
        self.mixed = self.Class(self.mixed.NamedSeqs, Names=['b','a'])    
        self.assertEqual(self.mixed.getItemsIf(is_vowel), ['O','A','E'])
        self.assertEqual(self.one_seq.getItemsIf(is_vowel), list('AAAAA'))
        #should be able to negate
        self.assertEqualItems(self.mixed.getItemsIf(is_vowel, negate=True), \
            list('BCDLMNP'))

    def test_getSimilar(self):
        """SequenceCollection getSimilar should get all sequences close to target seq"""
        aln = self.many
        x = RnaSequence('GGGGGGGGGG')
        y = RnaSequence('----------')
        #test min and max similarity ranges
        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0.4,\
            max_similarity=0.7)
        for seq in 'cefg':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 4)
        
        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0.95, \
            max_similarity=1)
        for seq in 'a':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 1)

        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0.75, \
            max_similarity=0.85)
        for seq in 'bd':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 2)

        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0, \
            max_similarity=0.2)
        self.assertEqual(result, {})

        #test some sequence transformations
        transform = lambda s: s[1:4]
        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0.5, \
            transform=transform)
        for seq in 'abdfg':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 5)

        transform = lambda s: s[-3:]
        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=0.5, \
            transform=transform)
        for seq in 'abcde':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 5)

        #test a different distance metric
        metric = lambda x, y: str(x).count('G') + str(y).count('G')
        result = aln.getSimilar(aln.NamedSeqs['a'], min_similarity=5, \
            max_similarity=10, metric=metric)
        for seq in 'ef':
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 2)

        #test the combination of a transform and a distance metric
        aln = self.Class(dict(enumerate(map(RnaSequence, \
            ['GA-GU','A-GAC','GG-GG']))), MolType=RNA)
        transform = lambda s: RnaSequence(str(s).replace('G','A'\
            ).replace('U','C'))
        metric = RnaSequence.fracSameNonGaps
        null_transform = lambda s: RnaSequence(str(s))
        #first, do it without the transformation
        try:
            result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.5, \
                metric=metric)
        except TypeError:   #need to coerce to RNA seq w/ null_transform
            result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.5, \
                metric=metric, transform=null_transform)
        for seq in [0,2]:
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 2)
        #repeat with higher similarity
        try:
            result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.8, \
                metric=metric)
        except TypeError:   #need to coerce to RNA
            result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.8, \
                metric=metric, transform=null_transform)
        for seq in [0]:
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 1)
        #then, verify that the transform changes the results         
        result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.5, \
            metric=metric, transform=transform)
        for seq in [0,1,2]:
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 3)
        
        result = aln.getSimilar(aln.NamedSeqs[0], min_similarity=0.8, \
            metric=metric, transform=transform)
        for seq in [0,1]:
            self.assertContains(result.NamedSeqs, seq)
            self.assertEquals(result.NamedSeqs[seq], aln.NamedSeqs[seq])
        self.assertEqual(len(result.NamedSeqs), 2)
         
    def test_distanceMatrix(self):
        """SequenceCollection distanceMatrix should produce correct scores"""
        self.assertEqual(self.one_seq.distanceMatrix(frac_same), {'a':{'a':1}})
        self.assertEqual(self.gaps.distanceMatrix(frac_same), 
            {   'a':{'a':7/7.0,'b':4/7.0,'c':2/7.0},
                'b':{'a':4/7.0,'b':7/7.0,'c':3/7.0},
                'c':{'a':2/7.0,'b':3/7.0,'c':7/7.0},
            })


    def test_isRagged(self):
        """SequenceCollection isRagged should return true if ragged alignment"""
        assert(not self.identical.isRagged())
        assert(not self.gaps.isRagged())


    def test_toPhylip(self):
        """SequenceCollection should return PHYLIP string format correctly"""
        align_norm = self.Class( ['ACDEFGHIKLMNPQRSTUVWY-',
                                     'ACDEFGHIKLMNPQRSUUVWF-',
                                     'ACDEFGHIKLMNPERSKUVWC-',
                                     'ACNEFGHIKLMNPQRS-UVWP-',                                     
                                     ])

        phylip_str, id_map =  align_norm.toPhylip()

        self.assertEqual(phylip_str, """4 22\nseq0000001 ACDEFGHIKLMNPQRSTUVWY-\nseq0000002 ACDEFGHIKLMNPQRSUUVWF-\nseq0000003 ACDEFGHIKLMNPERSKUVWC-\nseq0000004 ACNEFGHIKLMNPQRS-UVWP-""")
        self.assertEqual(id_map, {'seq0000004':'seq_3', 'seq0000001':'seq_0', \
            'seq0000003': 'seq_2', 'seq0000002': 'seq_1'})

    def test_toFasta(self):
        """SequenceCollection should return correct FASTA string"""
        aln = self.Class(['AAA','CCC'])
        self.assertEqual(aln.toFasta(), '>seq_0\nAAA\n>seq_1\nCCC')

        #NOTE THE FOLLOWING SURPRISING BEHAVIOR BECAUSE OF THE TWO-ITEM
        #SEQUENCE RULE:
        aln = self.Class(['AA','CC'])
        self.assertEqual(aln.toFasta(), '>A\nA\n>C\nC')

    def test_toNexus(self):
        """SequenceCollection should return correct Nexus string format"""
        align_norm = self.Class( ['ACDEFGHIKLMNPQRSTUVWY-',
                                  'ACDEFGHIKLMNPQRSUUVWF-',
                                  'ACDEFGHIKLMNPERSKUVWC-',
                                  'ACNEFGHIKLMNPQRS-UVWP-'])
        expect = '#NEXUS\n\nbegin data;\n    dimensions ntax=4 nchar=22;\n'+\
        '    format datatype=protein interleave=yes missing=? gap=-;\n'+\
        '    matrix\n    seq_1    ACDEFGHIKLMNPQRSUUVWF-\n    seq_0'+\
        '    ACDEFGHIKLMNPQRSTUVWY-\n    seq_3    ACNEFGHIKLMNPQRS-UVWP-\n   '+\
        ' seq_2    ACDEFGHIKLMNPERSKUVWC-\n\n    ;\nend;'
        self.assertEqual(align_norm.toNexus('protein'), expect)
    
    def test_getIntMap(self):
        """SequenceCollection.getIntMap should return correct mapping."""
        aln = self.Class({'seq1':'ACGU','seq2':'CGUA','seq3':'CCGU'})
        int_keys = {'seq_0':'seq1','seq_1':'seq2','seq_2':'seq3'}
        int_map = {'seq_0':'ACGU','seq_1':'CGUA','seq_2':'CCGU'}
        im,ik = aln.getIntMap()
        self.assertEqual(ik,int_keys)
        self.assertEqual(im,int_map)
        #test change prefix from default 'seq_'
        prefix='seqn_'
        int_keys = {'seqn_0':'seq1','seqn_1':'seq2','seqn_2':'seq3'}
        int_map = {'seqn_0':'ACGU','seqn_1':'CGUA','seqn_2':'CCGU'}
        im,ik = aln.getIntMap(prefix=prefix)
        self.assertEqual(ik,int_keys)
        self.assertEqual(im,int_map)

    def test_getNumSeqs(self):
        """SequenceCollection.getNumSeqs should count seqs."""
        aln = self.Class({'seq1':'ACGU','seq2':'CGUA','seq3':'CCGU'})
        self.assertEqual(aln.getNumSeqs(), 3)

    def test_copyAnnotations(self):
        """SequenceCollection copyAnnotations should copy from seq objects"""
        aln = self.Class({'seq1':'ACGU','seq2':'CGUA','seq3':'CCGU'})
        seq_1 = Sequence('ACGU', Name='seq1')
        seq_1.addFeature('xyz','abc', [(1,2)])
        seq_5 = Sequence('ACGUAAAAAA', Name='seq5')
        seq_5.addFeature('xyzzz','abc', [(1,2)])
        annot = {'seq1': seq_1, 'seq5':seq_5}
        aln.copyAnnotations(annot)
        aln_seq_1 = aln.NamedSeqs['seq1']
        if not hasattr(aln_seq_1, 'annotations'):
            aln_seq_1 = aln_seq_1.data
        aln_seq_2 = aln.NamedSeqs['seq2']
        if not hasattr(aln_seq_2, 'annotations'):
            aln_seq_2 = aln_seq_2.data
        self.assertEqual(len(aln_seq_1.annotations), 1)
        self.assertEqual(aln_seq_1.annotations[0].Name,'abc')
        self.assertEqual(len(aln_seq_2.annotations), 0)

    def test_annotateFromGff(self):
        """SequenceCollection.annotateFromGff should read gff features"""
        aln = self.Class({'seq1':'ACGU','seq2':'CGUA','seq3':'CCGU'})
        gff = [
            ['seq1', 'prog1', 'snp', '1', '2', '1.0', '+', '1','"abc"'],
            ['seq5', 'prog2', 'snp', '2', '3', '1.0', '+', '1','"yyy"'],
            ]
        gff = map('\t'.join, gff)
        aln.annotateFromGff(gff)
        aln_seq_1 = aln.NamedSeqs['seq1']
        if not hasattr(aln_seq_1, 'annotations'):
            aln_seq_1 = aln_seq_1.data
        aln_seq_2 = aln.NamedSeqs['seq2']
        if not hasattr(aln_seq_2, 'annotations'):
            aln_seq_2 = aln_seq_2.data
        self.assertEqual(len(aln_seq_1.annotations), 1)
        self.assertEqual(aln_seq_1.annotations[0].Name,'abc')
        self.assertEqual(len(aln_seq_2.annotations), 0)
          
    def test_replaceSeqs(self):
        """replaceSeqs should replace 1-letter w/ 3-letter seqs"""
        a = Alignment({'seq1':'ACGU','seq2':'C-UA','seq3':'C---'})
        seqs = {'seq1':'AAACCCGGGUUU','seq2':'CCCUUUAAA','seq3':'CCC'}
        result = a.replaceSeqs(seqs)
        self.assertEqual(result.toFasta(), \
        ">seq1\nAAACCCGGGUUU\n>seq2\nCCC---UUUAAA\n>seq3\nCCC---------")

    def test_getGappedSeq(self):
        """SequenceCollection.getGappedSeq should return seq, with gaps"""
        aln = self.Class({'seq1': '--TTT?', 'seq2': 'GATC??'})
        self.assertEqual(str(aln.getGappedSeq('seq1')), '--TTT?')
 
    def test_add(self):
        """__add__ should concatenate sequence data, by name"""
        align1= self.Class({'a': 'AAAA', 'b': 'TTTT', 'c': 'CCCC'})
        align2 = self.Class({'a': 'GGGG', 'b': '----', 'c': 'NNNN'})
        align = align1 + align2
        concatdict = align.todict()
        self.assertEqual(concatdict, {'a': 'AAAAGGGG', 'b': 'TTTT----', 'c': 'CCCCNNNN'})

    def test_addSeqs(self):
        """addSeqs should return an alignment with the new sequences appended"""
        a = [('s4', 'ACDEFGHIKLMNPQRSTUVWY-'), ('s3', 'ACDEFGHIKLMNPQRSUUVWF-')]
        b = [('s1', 'ACDEFGHIKLMNPERSKUVWC-'), ('s2', 'ACNEFGHIKLMNPQRS-UVWP-')]
        aln1 = self.Class(a)
        aln2 = self.Class(b)
        self.assertEqual(aln1.addSeqs(aln2).toFasta(), self.Class(a+b).toFasta())
        if isinstance(aln1, Alignment) or isinstance(aln1, DenseAlignment):
            self.assertRaises((DataError, ValueError), aln1.addSeqs, aln2+aln2)
        else:
            exp = set([seq for name, seq in a])
            exp.update([seq+seq for name, seq in b])
            got = set()
            for seq in aln1.addSeqs(aln2+aln2).Seqs:
                got.update([str(seq).strip()])
            self.assertEqual(got, exp)

    def test_writeToFile(self):
        """SequenceCollection.writeToFile should write in correct format"""
        aln = self.Class([('a','AAAA'),( 'b','TTTT'),('c','CCCC')])
        fn = mktemp(suffix='.fasta')
        aln.writeToFile(fn)
        result = open(fn, 'U').read()
        self.assertEqual(result, '>a\nAAAA\n>b\nTTTT\n>c\nCCCC\n')
        remove(fn)

    def test_len(self):
        """len(SequenceCollection) returns length of longest sequence"""
        aln = self.Class([('a','AAAA'),( 'b','TTTT'),('c','CCCC')])
        self.assertEqual(len(aln), 4) 

    def test_getTranslation(self):
        """SequenceCollection.getTranslation translates each seq"""
        for seqs in [
                {'seq1': 'GATTTT', 'seq2': 'GATC??'},
                {'seq1': 'GAT---', 'seq2': '?GATCT'}]:
            alignment = self.Class(data=seqs, MolType=DNA)
            self.assertEqual(len(alignment.getTranslation()), 2)
            # check for a failure when no moltype specified
            alignment = self.Class(data=seqs)
            try:
                peps = alignment.getTranslation()
            except AttributeError:
                pass

    def test_getSeq(self):
        """SequenceCollection.getSeq should return specified seq"""
        aln = self.Class({'seq1': 'GATTTT', 'seq2': 'GATC??'})
        self.assertEqual(aln.getSeq('seq1'), 'GATTTT')
        self.assertRaises(KeyError, aln.getSeq, 'seqx')
       

    def test_todict(self):
        """SequenceCollection.todict should return dict of strings (not obj)"""
        aln = self.Class({'seq1': 'GATTTT', 'seq2': 'GATC??'})
        self.assertEqual(aln.todict(), {'seq1':'GATTTT','seq2':'GATC??'})
        for i in aln.todict().values():
            assert isinstance(i, str)

    def test_getPerSequenceAmbiguousPositions(self):
        """SequenceCollection.getPerSequenceAmbiguousPositions should return pos"""
        aln = self.Class({'s1':'ATGRY?','s2':'T-AG??'}, MolType=DNA)
        self.assertEqual(aln.getPerSequenceAmbiguousPositions(), \
            {'s2': {4: '?', 5: '?'}, 's1': {3: 'R', 4: 'Y', 5: '?'}})

    def test_degap(self):
        """SequenceCollection.degap should strip gaps from each seq"""
        aln = self.Class({'s1':'ATGRY?','s2':'T-AG??'}, MolType=DNA)
        self.assertEqual(aln.degap(), {'s1':'ATGRY','s2':'TAG'})

    def test_withModifiedTermini(self):
        """SequenceCollection.withModifiedTermini should code trailing gaps as ?"""
        aln = self.Class({'s1':'AATGR--','s2':'-T-AG?-'}, MolType=DNA)
        self.assertEqual(aln.withModifiedTermini(), \
            {'s1':'AATGR??','s2':'?T-AG??'})



    def test_omitSeqsTemplate(self):
        """SequenceCollection.omitSeqsTemplate returns new aln with well-aln to temp"""
        aln = self.omitSeqsTemplate_aln
        result = aln.omitSeqsTemplate('s3', 0.9, 5)
        self.assertEqual(result, {'s3': 'UUCCUUCUU-UUC', \
                's4': 'UU-UUUU-UUUUC'})
        result2 = aln.omitSeqsTemplate('s4', 0.9, 4)
        self.assertEqual(result2, {'s3': 'UUCCUUCUU-UUC', \
                's4': 'UU-UUUU-UUUUC'})
        result3 = aln.omitSeqsTemplate('s1', 0.9, 4)
        self.assertEqual(result3, {'s2': 'UC------U---C', \
                's1': 'UC-----CU---C', 's5': '-------------'})
        result4 = aln.omitSeqsTemplate('s3', 0.5, 13)
        self.assertEqual(result4, {'s3': 'UUCCUUCUU-UUC', \
                's4': 'UU-UUUU-UUUUC'})
        
    def test_make_gap_filter(self):
        """make_gap_filter returns f(seq) -> True if aligned ok w/ query"""
        s1 = RnaSequence('UC-----CU---C')
        s3 = RnaSequence('UUCCUUCUU-UUC')
        s4 = RnaSequence('UU-UUUU-UUUUC')
        #check that the behavior is ok for gap runs
        f1 = make_gap_filter(s1, 0.9, 5)
        f3 = make_gap_filter(s3, 0.9, 5)
        #Should return False since s1 has gap run >= 5 with respect to s3
        self.assertEqual(f3(s1), False)
        #Should return False since s3 has an insertion run >= 5 to s1
        self.assertEqual(f1(s3), False)
        #Should retun True since s4 does not have a long enough gap or ins run
        self.assertEqual(f3(s4), True)
        f3 = make_gap_filter(s3, 0.9, 6)
        self.assertEqual(f3(s1), True)
        
        #Check that behavior is ok for gap_fractions
        f1 = make_gap_filter(s1, 0.5, 6)
        f3 = make_gap_filter(s3, 0.5, 6)
        #Should return False since 0.53% of positions are diff for gaps
        self.assertEqual(f3(s1), False)
        self.assertEqual(f1(s3), False)
        self.assertEqual(f3(s4), True)

    def test_omitGapSeqs(self):
        """SequenceCollection omitGapSeqs should return alignment w/o seqs with gaps"""
        #check default params
        self.assertEqual(self.gaps.omitGapSeqs(), self.gaps.omitGapSeqs(0))
        #check for boundary effects
        self.assertEqual(self.gaps.omitGapSeqs(-1), {})
        self.assertEqual(self.gaps.omitGapSeqs(0), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps.omitGapSeqs(0.1), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps.omitGapSeqs(3.0/7 - 0.01), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps.omitGapSeqs(3.0/7), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapSeqs(3.0/7 + 0.01), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapSeqs(5.0/7 - 0.01), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapSeqs(5.0/7 + 0.01), self.gaps)
        self.assertEqual(self.gaps.omitGapSeqs(0.99), self.gaps)
        #check new object creation
        self.assertNotSameObj(self.gaps.omitGapSeqs(0.99), self.gaps)
        self.assertTrue(isinstance(self.gaps.omitGapSeqs(3.0/7), 
                                   SequenceCollection))
        #repeat tests for object that supplies its own gaps
        self.assertEqual(self.gaps_rna.omitGapSeqs(-1), {})
        self.assertEqual(self.gaps_rna.omitGapSeqs(0), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(0.1), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(3.0/7 - 0.01), \
            {'a':'AAAAAAA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(3.0/7), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(3.0/7 + 0.01), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(5.0/7 - 0.01), \
            {'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps_rna.omitGapSeqs(5.0/7 + 0.01), self.gaps_rna)
        self.assertEqual(self.gaps_rna.omitGapSeqs(0.99), self.gaps_rna)
        self.assertNotSameObj(self.gaps_rna.omitGapSeqs(0.99), self.gaps_rna)
        self.assertTrue(isinstance(self.gaps_rna.omitGapSeqs(3.0/7), 
                                   SequenceCollection))

    def test_omitGapRuns(self):
        """SequenceCollection omitGapRuns should return alignment w/o runs of gaps"""
        #negative value will still let through ungapped sequences
        self.assertEqual(self.gaps.omitGapRuns(-5), {'a':'AAAAAAA'})
        #test edge effects
        self.assertEqual(self.gaps.omitGapRuns(0), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps.omitGapRuns(1), {'a':'AAAAAAA'})
        self.assertEqual(self.gaps.omitGapRuns(2),{'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapRuns(3),{'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapRuns(4),{'a':'AAAAAAA','b':'A--A-AA'})
        self.assertEqual(self.gaps.omitGapRuns(5), self.gaps)
        self.assertEqual(self.gaps.omitGapRuns(6), self.gaps)
        self.assertEqual(self.gaps.omitGapRuns(1000), self.gaps)
        #test new object creation
        self.assertNotSameObj(self.gaps.omitGapRuns(6), self.gaps)
        self.assertTrue(isinstance(self.gaps.omitGapRuns(6), 
                                   SequenceCollection))
    
    def test_consistent_gap_degen_handling(self):
        """gap degen character should be treated consistently"""
        # the degen character '?' can be a gap, so when we strip gaps it should
        # be gone too
        raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
        raw_ungapped = re.sub("[-?]", "", raw_seq)
        raw_no_ambigs = re.sub("[N?]+", "", raw_seq)
        dna = DNA.makeSequence(raw_seq)
        
        aln = self.Class(data=[("a", dna),("b", dna)])
        expect = self.Class(data=[("a", raw_ungapped),("b", raw_ungapped)]).toFasta()
        self.assertEqual(aln.degap().toFasta(), expect)
        seqs = self.Class(data=[("a", dna),("b", dna)])
        self.assertEqual(seqs.degap().toFasta(), expect)
    
    def test_padSeqs(self):
        """SequenceCollection padSeqs should work on alignment."""
        #pad to max length
        padded1 = self.ragged_padded.padSeqs()
        seqs1 = list(padded1.iterSeqs(seq_order=['a','b','c']))
        self.assertEqual(map(str,seqs1),['AAAAAA', 'AAA---', 'AAAA--'])
        
        #pad to alternate length
        padded1 = self.ragged_padded.padSeqs(pad_length=10)
        seqs1 = list(padded1.iterSeqs(seq_order=['a','b','c']))
        self.assertEqual(map(str,seqs1),['AAAAAA----', 'AAA-------',\
            'AAAA------'])
        
        #assertRaises error when pad_length is less than max seq length
        self.assertRaises(ValueError, self.ragged_padded.padSeqs, 5)

class SequenceCollectionTests(SequenceCollectionBaseTests, TestCase):
    """Tests of the SequenceCollection object. Includes ragged collection tests.

    Should not test alignment-specific features.
    """

    def setUp(self):
        """Adds self.ragged for ragged collection tests."""
        self.ragged = SequenceCollection({'a':'AAAAAA', 'b':'AAA', 'c':'AAAA'})
        super(SequenceCollectionTests, self).setUp()

    def test_SeqLen_get_ragged(self):
        """SequenceCollection SeqLen get should work for ragged seqs"""
        self.assertEqual(self.ragged.SeqLen, 6)

    def test_isRagged_ragged(self):
        """SequenceCollection isRagged should return True if ragged"""
        self.assertTrue(self.ragged.isRagged())

    def test_Seqs_ragged(self):
        """SequenceCollection Seqs should work on ragged alignment"""
        self.ragged.Names = 'bac'
        self.assertEqual(list(self.ragged.Seqs), ['AAA', 'AAAAAA', 'AAAA'])

    def test_iterSeqs_ragged(self):
        """SequenceCollection iterSeqs() method should support reordering of seqs"""
        self.ragged.Names = ['a','b','c']
        seqs = list(self.ragged.iterSeqs())
        self.assertEqual(seqs, ['AAAAAA', 'AAA', 'AAAA'])
        seqs = list(self.ragged.iterSeqs(seq_order=['b','a','a']))
        self.assertEqual(seqs, ['AAA', 'AAAAAA', 'AAAAAA'])
        self.assertSameObj(seqs[1], seqs[2])
        self.assertSameObj(seqs[0], self.ragged.NamedSeqs['b'])

    def test_toPHYLIP_ragged(self):
        """SequenceCollection should refuse to convert ragged seqs to phylip"""
        align_rag = self.Class( ['ACDEFGHIKLMNPQRSTUVWY-',
                                     'ACDEFGHIKLMNPQRSUUVWF-',
                                     'ACDEFGHIKLMNPERSKUVWC-',
                                     'ACNEFGHIKLMNUVWP-',                                     
                                     ])


        self.assertRaises(ValueError,  align_rag.toPhylip)
    
    def test_padSeqs_ragged(self):
        """SequenceCollection padSeqs should work on ragged alignment."""
        #pad to max length
        padded1 = self.ragged.padSeqs()
        seqs1 = list(padded1.iterSeqs(seq_order=['a','b','c']))
        self.assertEqual(map(str,seqs1),['AAAAAA', 'AAA---', 'AAAA--'])
        
        #pad to alternate length
        padded1 = self.ragged.padSeqs(pad_length=10)
        seqs1 = list(padded1.iterSeqs(seq_order=['a','b','c']))
        self.assertEqual(map(str,seqs1),['AAAAAA----', 'AAA-------',\
            'AAAA------'])
        
        #assertRaises error when pad_length is less than max seq length
        self.assertRaises(ValueError, self.ragged.padSeqs, 5)

class AlignmentBaseTests(SequenceCollectionBaseTests):
    """Tests of basic Alignment functionality. All Alignments should pass these.

    Note that this is not a TestCase: need to subclass to test each specific
    type of Alignment. Override self.Constructor with your alignment class
    as a constructor.
    """
    def test_Positions(self):
        """SequenceCollection Positions property should iterate over positions, using self.Names"""
        r = self.Class({'a':'AAAAAA','b':'AAA---','c':'AAAA--'})
        r.Names = ['a','b','c']
        self.assertEqual(list(r.Positions), map(list, \
            ['AAA','AAA','AAA', 'A-A', 'A--', 'A--']))
        
    def test_iterPositions(self):
        #"""SequenceCollection iterPositions() method should support reordering of #cols"""
        r = self.Class(self.ragged_padded.NamedSeqs, Names=['c','b'])
        self.assertEqual(list(r.iterPositions(pos_order=[5,1,3])),\
            map(list,['--','AA','A-']))
        #reorder names
        r = self.Class(self.ragged_padded.NamedSeqs, Names=['a','b','c'])
        cols = list(r.iterPositions())
        self.assertEqual(cols, map(list, ['AAA','AAA','AAA','A-A','A--','A--']))
 
    def test_takePositions(self):
        """SequenceCollection takePositions should return new alignment w/ specified pos"""
        self.assertEqual(self.gaps.takePositions([5,4,0], \
            seq_constructor=coerce_to_string), \
            {'a':'AAA','b':'A-A','c':'--A'})
        self.assertTrue(isinstance(self.gaps.takePositions([0]), 
                                   SequenceCollection))
        #should be able to negate
        self.assertEqual(self.gaps.takePositions([5,4,0], negate=True, \
            seq_constructor=coerce_to_string),
            {'a':'AAAA','b':'--AA','c':'A---'})
    
    def test_getPositionIndices(self):
        """SequenceCollection getPositionIndices should return names of cols where f(col)"""
        gap_1st = lambda x: x[0] == '-'
        gap_2nd = lambda x: x[1] == '-'
        gap_3rd = lambda x: x[2] == '-'
        is_list =  lambda x: isinstance(x, list)
        self.gaps = self.Class(self.gaps.NamedSeqs, Names=['a','b','c'])

        self.assertEqual(self.gaps.getPositionIndices(gap_1st), [])
        self.assertEqual(self.gaps.getPositionIndices(gap_2nd), [1,2,4])
        self.assertEqual(self.gaps.getPositionIndices(gap_3rd), [2,3,4,5,6])
        self.assertEqual(self.gaps.getPositionIndices(is_list), [0,1,2,3,4,5,6])
        #should be able to negate
        self.assertEqual(self.gaps.getPositionIndices(gap_2nd, negate=True), \
            [0,3,5,6])
        self.assertEqual(self.gaps.getPositionIndices(gap_1st, negate=True), \
            [0,1,2,3,4,5,6])
        self.assertEqual(self.gaps.getPositionIndices(is_list, negate=True), [])

    def test_takePositionsIf(self):
        """SequenceCollection takePositionsIf should return cols where f(col) is True"""
        gap_1st = lambda x: x[0] == '-'
        gap_2nd = lambda x: x[1] == '-'
        gap_3rd = lambda x: x[2] == '-'
        is_list =  lambda x: isinstance(x, list)
        self.gaps.Names = 'abc'

        self.assertEqual(self.gaps.takePositionsIf(gap_1st,seq_constructor=coerce_to_string),\
            {'a':'', 'b':'', 'c':''})
        self.assertEqual(self.gaps.takePositionsIf(gap_2nd,seq_constructor=coerce_to_string),\
            {'a':'AAA', 'b':'---', 'c':'A--'})
        self.assertEqual(self.gaps.takePositionsIf(gap_3rd,seq_constructor=coerce_to_string),\
            {'a':'AAAAA', 'b':'-A-AA', 'c':'-----'})
        self.assertEqual(self.gaps.takePositionsIf(is_list,seq_constructor=coerce_to_string),\
            self.gaps)

        self.assertTrue(isinstance(self.gaps.takePositionsIf(gap_1st), 
                                   SequenceCollection))
        #should be able to negate
        self.assertEqual(self.gaps.takePositionsIf(gap_1st, seq_constructor=coerce_to_string,\
            negate=True), self.gaps)
        self.assertEqual(self.gaps.takePositionsIf(gap_2nd, seq_constructor=coerce_to_string,\
            negate=True), {'a':'AAAA','b':'AAAA','c':'A---'})
        self.assertEqual(self.gaps.takePositionsIf(gap_3rd, seq_constructor=coerce_to_string,\
            negate=True), {'a':'AA','b':'A-','c':'AA'})
    
    def test_omitGapPositions(self):
        """SequenceCollection omitGapPositions should return alignment w/o positions of gaps"""
        aln = self.end_gaps
        
        #first, check behavior when we're just acting on the cols (and not
        #trying to delete the naughty seqs).
        
        #default should strip out cols that are 100% gaps
        self.assertEqual(aln.omitGapPositions(seq_constructor=coerce_to_string), \
            {'a':'-ABC', 'b':'CBA-', 'c':'-DEF'})
        #if allowed_gap_frac is 1, shouldn't delete anything
        self.assertEqual(aln.omitGapPositions(1, seq_constructor=coerce_to_string), \
            {'a':'--A-BC-', 'b':'-CB-A--', 'c':'--D-EF-'})
        #if allowed_gap_frac is 0, should strip out any cols containing gaps
        self.assertEqual(aln.omitGapPositions(0, seq_constructor=coerce_to_string), \
            {'a':'AB', 'b':'BA', 'c':'DE'})
        #intermediate numbers should work as expected
        self.assertEqual(aln.omitGapPositions(0.4, seq_constructor=coerce_to_string), \
            {'a':'ABC', 'b':'BA-', 'c':'DEF'})
        self.assertEqual(aln.omitGapPositions(0.7, seq_constructor=coerce_to_string), \
            {'a':'-ABC', 'b':'CBA-', 'c':'-DEF'})

        #second, need to check behavior when the naughty seqs should be
        #deleted as well.

        #default should strip out cols that are 100% gaps
        self.assertEqual(aln.omitGapPositions(seq_constructor=coerce_to_string, \
            del_seqs=True), {'a':'-ABC', 'b':'CBA-', 'c':'-DEF'})
        #if allowed_gap_frac is 1, shouldn't delete anything
        self.assertEqual(aln.omitGapPositions(1, seq_constructor=coerce_to_string, \
            del_seqs=True), {'a':'--A-BC-', 'b':'-CB-A--', 'c':'--D-EF-'})
        #if allowed_gap_frac is 0, should strip out any cols containing gaps
        self.assertEqual(aln.omitGapPositions(0, seq_constructor=coerce_to_string, \
            del_seqs=True), {}) #everything has at least one naughty non-gap
        #intermediate numbers should work as expected
        self.assertEqual(aln.omitGapPositions(0.4, seq_constructor=coerce_to_string,
            del_seqs=True), {'a':'ABC', 'c':'DEF'}) #b has a naughty non-gap
        #check that does not delete b if allowed_frac_bad_calls higher than 0.14
        self.assertEqual(aln.omitGapPositions(0.4, seq_constructor=coerce_to_string,
            del_seqs=True, allowed_frac_bad_cols=0.2), \
                    {'a':'ABC', 'b':'BA-','c':'DEF'})
        self.assertEqual(aln.omitGapPositions(0.4, seq_constructor=coerce_to_string,
            del_seqs=True), {'a':'ABC', 'c':'DEF'}) #b has a naughty non-gap
        
        self.assertEqual(aln.omitGapPositions(0.7, seq_constructor=coerce_to_string,
            del_seqs=True), {'a':'-ABC', 'b':'CBA-', 'c':'-DEF'}) #all ok

        #when we increase the number of sequences to 6, more differences
        #start to appear.
        new_aln_data = aln.NamedSeqs.copy()
        new_aln_data['d'] = '-------'
        new_aln_data['e'] = 'XYZXYZX'
        new_aln_data['f'] = 'AB-CDEF'
        aln = self.Class(new_aln_data)

        #if no gaps are allowed, everything is deleted...
        result = aln.omitGapPositions(seq_constructor=coerce_to_string)
        self.assertEqual(aln.omitGapPositions(0, del_seqs=False), \
                {'a':'', 'b':'', 'c':'', 'd':'', 'e':'', 'f':''})
        #...though not a sequence that's all gaps, since it has no positions
        #that are not gaps. This 'feature' should possibly be considered a bug.
        self.assertEqual(aln.omitGapPositions(0, del_seqs=True), {'d':''})
        #if we're deleting only full positions of gaps, del_seqs does nothing.
        self.assertEqual(aln.omitGapPositions(del_seqs=True, \
            seq_constructor=coerce_to_string), aln)
        #at 50%, should delete a bunch of minority sequences
        self.assertEqual(aln.omitGapPositions(0.5, del_seqs=True, \
            seq_constructor=coerce_to_string), \
            {'a':'-ABC','b':'CBA-','c':'-DEF','d':'----'})
        #shouldn't depend on order of seqs
        aln.Names = 'fadbec'
        self.assertEqual(aln.omitGapPositions(0.5, del_seqs=True, \
            seq_constructor=coerce_to_string), \
            {'a':'-ABC','b':'CBA-','c':'-DEF','d':'----'})
        

    def test_IUPACConsensus_RNA(self):
        """SequenceCollection IUPACConsensus should use RNA IUPAC symbols correctly"""
        alignmentUpper = self.Class( ['UCAGN-UCAGN-UCAGN-UCAGAGCAUN-',
                                     'UUCCAAGGNN--UUCCAAGGNNAGCAG--',
                                     'UUCCAAGGNN--UUCCAAGGNNAGCUA--',
                                     'UUUUCCCCAAAAGGGGNNNN--AGCUA--',
                                     'UUUUCCCCAAAAGGGGNNNN--AGCUA--',
                                     ], MolType=RNA)
        
        #following IUPAC consensus calculated by hand
        #Test all uppper
        self.assertEqual(alignmentUpper.IUPACConsensus(),
                         'UYHBN?BSNN??KBVSN?NN??AGCWD?-')
    
    def test_IUPACConsensus_DNA(self):
        """SequenceCollection IUPACConsensus should use DNA IUPAC symbols correctly"""
        alignmentUpper = self.Class( ['TCAGN-TCAGN-TCAGN-TCAGAGCATN-',
                                     'TTCCAAGGNN--TTCCAAGGNNAGCAG--',
                                     'TTCCAAGGNN--TTCCAAGGNNAGCTA--',
                                     'TTTTCCCCAAAAGGGGNNNN--AGCTA--',
                                     'TTTTCCCCAAAAGGGGNNNN--AGCTA--',
                                     ])
        #following IUPAC consensus calculated by hand
        #Test all uppper
        self.assertEqual(alignmentUpper.IUPACConsensus(DNA),
                         'TYHBN?BSNN??KBVSN?NN??AGCWD?-')

    def test_IUPACConsensus_Protein(self):
        """SequenceCollection IUPACConsensus should use protein IUPAC symbols correctly"""
        alignmentUpper = self.Class( ['ACDEFGHIKLMNPQRSTUVWY-',
                                     'ACDEFGHIKLMNPQRSUUVWF-',
                                     'ACDEFGHIKLMNPERSKUVWC-',
                                     'ACNEFGHIKLMNPQRS-UVWP-',                                     
                                     ])
        #following IUPAC consensus calculated by hand
        #Test all uppper
        self.assertEqual(alignmentUpper.IUPACConsensus(PROTEIN),
                         'ACBEFGHIKLMNPZRS?UVWX-')
    def test_isRagged(self):
        """SequenceCollection isRagged should return true if ragged alignment"""
        assert(not self.identical.isRagged())
        assert(not self.gaps.isRagged())
    
    def test_columnProbs(self):
        """SequenceCollection.columnProbs should find Pr(symbol) in each column"""
        #make an alignment with 4 seqs (easy to calculate probabilities)
        align = self.Class(["AAA", "ACA", "GGG", "GUC"])
        cp = align.columnProbs()
        #check that the column probs match the counts we expect
        self.assertEqual(cp, map(Freqs, [   
            {'A':0.5, 'G':0.5},
            {'A':0.25, 'C':0.25, 'G':0.25, 'U':0.25},
            {'A':0.5, 'G':0.25, 'C':0.25},
            ]))

    def test_majorityConsensus(self):
        """SequenceCollection.majorityConsensus should return commonest symbol per column"""
        #Check the exact strings expected from string transform
        self.assertEqual(self.sequences.majorityConsensus(str), 'UCAG')
        self.assertEqual(self.structures.majorityConsensus(str), '(.....')

    
    def test_uncertainties(self):
        """SequenceCollection.uncertainties should match hand-calculated values"""
        aln = self.Class(['ABC', 'AXC'])
        obs = aln.uncertainties()
        self.assertFloatEqual(obs, [0, 1, 0])
        #check what happens with only one input sequence
        aln = self.Class(['ABC'])
        obs = aln.uncertainties()
        self.assertFloatEqual(obs, [0, 0, 0])
        #check that we can screen out bad items OK
        aln = self.Class(['ABC', 'DEF', 'GHI', 'JKL', '333'], MolType=BYTES)
        obs = aln.uncertainties('ABCDEFGHIJKLMNOP')
        self.assertFloatEqual(obs, [2.0] * 3)

    def test_columnFreqs(self):
        """Alignment.columnFreqs should count symbols in each column"""
        #calculate by hand what the first and last positions should look like in
        #each case
        firstvalues = [ 
                        [self.sequences, Freqs('UUU')],
                        [self.structures, Freqs('(.(')],
                    ]
        
        lastvalues = [ 
                        [self.sequences, Freqs('GGG')],
                        [self.structures, Freqs('..)')],
                    ]
        #check that the first positions are what we expected
        for obj, result in firstvalues:
            freqs = obj.columnFreqs()
            self.assertEqual(str(freqs[0]), str(result))
        #check that the last positions are what we expected
        for obj, result in lastvalues:
            freqs = obj.columnFreqs()
            self.assertEqual(str(freqs[-1]), str(result))

    def test_scoreMatrix(self):
        """Alignment scoreMatrix should produce position specific score matrix."""
        scoreMatrix = {
            0:{'A':1.0,'C':1.0,'U':5.0},
            1:{'C':6.0,'U':1.0},
            2:{'A':3.0,'C':2.0,'G':2.0},
            3:{'A':3.0,'G':4.0},
            4:{'C':1.0,'G':1.0,'U':5.0},
            5:{'C':6.0,'U':1.0},
            6:{'A':3.0,'G':4.0},
            7:{'A':1.0,'G':6.0},
            8:{'A':1.0,'C':1.0,'G':1.0,'U':4.0},
            9:{'A':1.0,'C':2.0,'U':4.0},
            }
        self.assertEqual(self.many.scoreMatrix(), scoreMatrix)


    def test_sample(self):
        """Alignment.sample should permute alignment by default"""
        alignment = self.Class({'seq1': 'ABCDEFGHIJKLMNOP',
                                    'seq2': 'ABCDEFGHIJKLMNOP'})
        # effectively permute columns, preserving length
        shuffled = alignment.sample()
        # ensure length correct
        sample = alignment.sample(10)
        self.assertEqual(len(sample), 10)
        # test columns alignment preserved
        seqs = sample.todict().values()
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs once as sampling without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 1)

    def test_sample_with_replacement(self):
        #test with replacement -- just verify that it rnus
        alignment = self.Class({'seq1': 'gatc', 'seq2': 'gatc'})
        sample = alignment.sample(1000, with_replacement=True)
        self.assertEqual(len(sample), 1000)
        # ensure that sampling with replacement works on single col alignment
        alignment1 = self.Class({'seq1': 'A',
                                    'seq2': 'A'})
        result = alignment1.sample(with_replacement=True)
        self.assertEqual(len(result), 1)

    def test_sample_tuples(self):
        ##### test with motif size != 1 #####
        alignment = self.Class({'seq1': 'AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP',
                                    'seq2': 'AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP'})
        shuffled = alignment.sample(motif_length=2)
        # ensure length correct
        sample = alignment.sample(10,motif_length=2)
        self.assertEqual(len(sample), 20)
        # test columns alignment preserved
        seqs = sample.todict().values()
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs twice as sampling dinucs without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 2) 


class DenseAlignmentTests(AlignmentBaseTests, TestCase):
    Class = DenseAlignment

    def test_get_freqs(self):
        """DenseAlignment getSeqFreqs: should work on positions and sequences 
        """
        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        da = DenseAlignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        seq_exp = array([[1,1,1,1],[0,3,1,0],[1,0,2,1]])
        pos_exp = array([[1,1,1,0],[0,2,0,1],[0,0,3,0],[1,1,0,1]])
        self.assertEqual(da._get_freqs(index=1), pos_exp)
        self.assertEqual(da._get_freqs(index=0), seq_exp)
        
    def test_getSeqFreqs(self):
        """DenseAlignment getSeqFreqs: should work with DnaSequences and strings
        """
        exp = array([[1,1,1,1],[0,3,1,0],[1,0,2,1]])

        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        da = DenseAlignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = da.getSeqFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

        s1 = 'TCAG'
        s2 = 'CCAC'
        s3 = 'AGAT'
        da = DenseAlignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = da.getSeqFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

    def test_getPosFreqs_sequence(self):
        """DenseAlignment getPosFreqs: should work with DnaSequences and strings
        """
        exp = array([[1,1,1,0],[0,2,0,1],[0,0,3,0],[1,1,0,1]])

        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        da = DenseAlignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = da.getPosFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

        s1 = 'TCAG'
        s2 = 'CCAC'
        s3 = 'AGAT'
        da = DenseAlignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = da.getPosFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))


class AlignmentTests(AlignmentBaseTests, TestCase):
    Class = Alignment

    def test_get_freqs(self):
        """Alignment _get_freqs: should work on positions and sequences 
        """
        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        aln = Alignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        seq_exp = array([[1,1,1,1],[0,3,1,0],[1,0,2,1]])
        pos_exp = array([[1,1,1,0],[0,2,0,1],[0,0,3,0],[1,1,0,1]])
        self.assertEqual(aln._get_freqs(index=1), pos_exp)
        self.assertEqual(aln._get_freqs(index=0), seq_exp)

    def test_getSeqFreqs(self):
        """Alignment getSeqFreqs: should work with DnaSequences and strings
        """
        exp = array([[1,1,1,1],[0,3,1,0],[1,0,2,1]])

        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        aln = Alignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = aln.getSeqFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

        s1 = 'TCAG'
        s2 = 'CCAC'
        s3 = 'AGAT'
        aln = Alignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = aln.getSeqFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

    def test_getPosFreqs(self):
        """Alignment getPosFreqs: should work with DnaSequences and strings
        """
        exp = array([[1,1,1,0],[0,2,0,1],[0,0,3,0],[1,1,0,1]])

        s1 = DNA.Sequence('TCAG', Name='s1')
        s2 = DNA.Sequence('CCAC', Name='s2')
        s3 = DNA.Sequence('AGAT', Name='s3')
        aln = Alignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = aln.getPosFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))

        s1 = 'TCAG'
        s2 = 'CCAC'
        s3 = 'AGAT'
        aln = Alignment([s1,s2,s3], MolType=DNA, Alphabet=DNA.Alphabet)
        obs = aln.getPosFreqs()
        self.assertEqual(obs.Data, exp)
        self.assertEqual(obs.Alphabet, DNA.Alphabet)
        self.assertEqual(obs.CharOrder, list("TCAG"))


    def make_and_filter(self, raw, expected, motif_length):
        # a simple filter func
        func = lambda x: re.findall("[-N?]", " ".join(x)) == []
        aln = self.Class(raw)
        result = aln.filtered(func,motif_length=motif_length,log_warnings=False)
        self.assertEqual(result.todict(), expected)
    
    def test_filtered(self):
        """filtered should return new alignment with positions consistent with
        provided callback function"""
        # a simple filter option
        raw = {'a':'ACGACGACG',
               'b':'CCC---CCC',
               'c':'AAAA--AAA'}
        self.make_and_filter(raw, {'a':'ACGACG','b':'CCCCCC','c':'AAAAAA'}, 1)
        # check with motif_length = 2
        self.make_and_filter(raw, {'a':'ACAC','b':'CCCC','c':'AAAA'}, 2)
        # check with motif_length = 3
        self.make_and_filter(raw, {'a':'ACGACG','b':'CCCCCC','c':'AAAAAA'}, 3)
        
    def test_slidingWindows(self):
        """slidingWindows should return slices of alignments."""
        alignment = self.Class({'seq1': 'ACGTACGT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'})
        result = []
        for bit in alignment.slidingWindows(5,2):
            result+=[bit]
        self.assertEqual(result[0].todict(), {'seq3': 'ACGTA', 'seq2': 'ACGTA', 'seq1': 'ACGTA'})
        self.assertEqual(result[1].todict(), {'seq3': 'GTACG', 'seq2': 'GTACG', 'seq1': 'GTACG'})

        result = []
        for bit in alignment.slidingWindows(5,1):
            result+=[bit]
        self.assertEqual(result[0].todict(), {'seq3': 'ACGTA', 'seq2': 'ACGTA', 'seq1': 'ACGTA'})
        self.assertEqual(result[1].todict(), {'seq3': 'CGTAC', 'seq2': 'CGTAC', 'seq1': 'CGTAC'})
        self.assertEqual(result[2].todict(), {'seq3': 'GTACG', 'seq2': 'GTACG', 'seq1': 'GTACG'})
        self.assertEqual(result[3].todict(), {'seq3': 'TACGT', 'seq2': 'TACGT', 'seq1': 'TACGT'})
 
    def test_withGapsFrom(self):
        """withGapsFrom should overwrite with gaps."""
        gapless   = self.Class({'seq1': 'TCG', 'seq2': 'TCG'})
        pregapped = self.Class({'seq1': '-CG', 'seq2': 'TCG'})
        template  = self.Class({'seq1': 'A-?', 'seq2': 'ACG'})
        r1 = gapless.withGapsFrom(template).todict()
        r2 = pregapped.withGapsFrom(template).todict()
        self.assertEqual(r1, {'seq1': 'T-G', 'seq2': 'TCG'}) 
        self.assertEqual(r2, {'seq1': '--G', 'seq2': 'TCG'}) 
        
class DenseAlignmentSpecificTests(TestCase):
    """Tests of the DenseAlignment object and its methods"""

    def setUp(self):
        """Define some standard alignments."""
        self.a = DenseAlignment(array([[0,1,2],[3,4,5]]), \
            conversion_f=aln_from_array)
        self.a2 = DenseAlignment(['ABC','DEF'], Names=['x','y'])
        class ABModelSequence(ModelSequence):
            Alphabet = AB.Alphabet
        self.ABModelSequence = ABModelSequence
        self.a = DenseAlignment(map(ABModelSequence, ['abaa','abbb']), \
            Alphabet=AB.Alphabet)
        self.b = Alignment(['ABC','DEF'])
        self.c = SequenceCollection(['ABC','DEF'])
    
    def test_init(self):
        """DenseAlignment init should work from a sequence"""
        a = DenseAlignment(array([[0,1,2],[3,4,5]]), conversion_f=aln_from_array)
        self.assertEqual(a.SeqData, array([[0,3],[1,4],[2,5]], 'B'))
        self.assertEqual(a.ArrayPositions, array([[0,1,2],[3,4,5]], 'B'))
        self.assertEqual(a.Names, ['seq_0','seq_1','seq_2'])
    
    def test_guess_input_type(self):
        """DenseAlignment _guess_input_type should figure out data type correctly"""
        git = self.a._guess_input_type
        self.assertEqual(git(self.a), 'dense_aln')
        self.assertEqual(git(self.b), 'aln')
        self.assertEqual(git(self.c), 'collection')
        self.assertEqual(git('>ab\nabc'), 'fasta')
        self.assertEqual(git(['>ab','abc']), 'fasta')
        self.assertEqual(git(['abc','def']), 'generic')
        self.assertEqual(git([[1,2],[4,5]]), 'kv_pairs') #precedence over generic
        self.assertEqual(git([[1,2,3],[4,5,6]]), 'generic')
        self.assertEqual(git([ModelSequence('abc')]), 'model_seqs')
        self.assertEqual(git(array([[1,2,3],[4,5,6]])), 'array')
        self.assertEqual(git({'a':'aca'}), 'dict')
        self.assertEqual(git([]), 'empty')

    def test_init_seqs(self):
        """DenseAlignment init should work from ModelSequence objects."""
        s = map(ModelSequence, ['abc','def'])
        a = DenseAlignment(s)
        self.assertEqual(a.SeqData, array(['abc','def'], 'c').view('B'))

    def test_init_generic(self):
        """DenseAlignment init should work from generic objects."""
        s = ['abc','def']
        a = DenseAlignment(s)
        self.assertEqual(a.SeqData, array(['abc','def'], 'c').view('B'))

    def test_init_aln(self):
        """DenseAlignment init should work from another alignment."""
        s = ['abc','def']
        a = DenseAlignment(s)
        b = DenseAlignment(a)
        self.assertNotSameObj(a.SeqData, b.SeqData)
        self.assertEqual(b.SeqData, array(['abc','def'], 'c').view('B'))

    def test_init_dict(self):
        """DenseAlignment init should work from dict."""
        s = {'abc':'aaaccc','xyz':'gcgcgc'}
        a = DenseAlignment(s)
        self.assertEqual(a.SeqData, array(['aaaccc','gcgcgc'], 'c').view('B'))
        self.assertEqual(tuple(a.Names), ('abc','xyz'))

    def test_init_empty(self):
        """DenseAlignment init should fail if empty."""
        self.assertRaises(TypeError, DenseAlignment)
        self.assertRaises(ValueError, DenseAlignment, 3)
    
    def test_get_alphabet_and_moltype(self):
        """DenseAlignment should figure out correct alphabet and moltype"""
        s1 = 'A'
        s2 = RNA.Sequence('AA')
        
        d = DenseAlignment(s1)
        self.assertSameObj(d.MolType, BYTES)
        self.assertSameObj(d.Alphabet, BYTES.Alphabet)
        
        d = DenseAlignment(s1, MolType=RNA)
        self.assertSameObj(d.MolType, RNA)
        self.assertSameObj(d.Alphabet, RNA.Alphabets.DegenGapped)
        
        d = DenseAlignment(s1, Alphabet=RNA.Alphabet)
        self.assertSameObj(d.MolType, RNA)
        self.assertSameObj(d.Alphabet, RNA.Alphabet)
        
        d = DenseAlignment(s2)
        self.assertSameObj(d.MolType, RNA)
        self.assertSameObj(d.Alphabet, RNA.Alphabets.DegenGapped)
        
        d = DenseAlignment(s2, MolType=DNA)
        self.assertSameObj(d.MolType, DNA)
        self.assertSameObj(d.Alphabet, DNA.Alphabets.DegenGapped)
        #checks for containers
        d = DenseAlignment([s2])
        self.assertSameObj(d.MolType, RNA)
        d = DenseAlignment({'x':s2})
        self.assertSameObj(d.MolType, RNA)
        d = DenseAlignment(set([s2]))
        self.assertSameObj(d.MolType, RNA)
    
    def test_iter(self):
        """DenseAlignment iter should iterate over positions"""
        result = list(iter(self.a2))
        for i, j in zip(result, [list(i) for i in ['AD', 'BE', 'CF']]):
            self.assertEqual(i,j)

    def test_getitem(self):
        """DenseAlignment getitem should default to positions as chars"""
        a2 = self.a2
        self.assertEqual(a2[1], ['B','E'])
        self.assertEqual(a2[1:], [['B','E'],['C','F']])

    def test_getSubAlignment(self):
        """DenseAlignment getSubAlignment should get requested part of alignment."""
        a = DenseAlignment('>x ABCE >y FGHI >z JKLM'.split())
        #passing in positions should keep all seqs, but just selected positions
        b = DenseAlignment('>x BC >y GH >z KL'.split())
        a_1 = a.getSubAlignment(pos=[1,2])
        self.assertEqual(a_1.Names, b.Names)
       
        self.assertEqual(a_1.Seqs, b.Seqs)
        #...and with invert_pos, should keep all except the positions passed in
        a_2 = a.getSubAlignment(pos=[0,3], invert_pos=True)
        self.assertEqual(a_2.Seqs, b.Seqs)
        self.assertEqual(a_2.Names, b.Names)
        #passing in seqs should keep all positions, but just selected seqs
        c = DenseAlignment('>x ABCE >z JKLM'.split())
        a_3 = a.getSubAlignment(seqs=[0,2])
        self.assertEqual(a_3.Seqs, c.Seqs)
        #check that labels were updates as well...
        self.assertEqual(a_3.Names, c.Names)
        #...and should work with invert_seqs to exclude just selected seqs
        a_4 = a.getSubAlignment(seqs=[1], invert_seqs=True)
        self.assertEqual(a_4.Seqs, c.Seqs)
        self.assertEqual(a_4.Names, c.Names)
        #should be able to do both seqs and positions simultaneously
        d = DenseAlignment('>x BC >z KL'.split())
        a_5 = a.getSubAlignment(seqs=[0,2], pos=[1,2])
        self.assertEqual(a_5.Seqs, d.Seqs)
        self.assertEqual(a_5.Names, d.Names)

    def test_str(self):
        """DenseAlignment str should return FASTA representation of aln"""
        self.assertEqual(str(self.a2), '>x\nABC\n>y\nDEF\n')
        #should work if labels diff length
        self.a2.Names[-1] = 'yyy'
        self.assertEqual(str(self.a2), '>x\nABC\n>yyy\nDEF\n')
    
    def test_get_freqs(self):
        """DenseAlignment _get_freqs should get row or col freqs"""
        ABModelSequence = self.ABModelSequence
        a = self.a
        self.assertEqual(a._get_freqs(0), array([[3,1],[1,3]]))
        self.assertEqual(a._get_freqs(1), array([[2,0],[0,2],[1,1],[1,1]]))
        
    def test_getSeqFreqs(self):
        """DenseAlignment getSeqFreqs should get profile of freqs in each seq"""
        ABModelSequence = self.ABModelSequence
        a = self.a
        f = a.getSeqFreqs()
        self.assertEqual(f.Data, array([[3,1],[1,3]]))

    def test_getPosFreqs(self):
        """DenseAlignment getPosFreqs should get profile of freqs at each pos"""
        ABModelSequence = self.ABModelSequence
        a = self.a
        f = a.getPosFreqs()
        self.assertEqual(f.Data, array([[2,0],[0,2],[1,1],[1,1]]))

    def test_getSeqEntropy(self):
        """DenseAlignment getSeqEntropy should get entropy of each seq"""
        ABModelSequence = self.ABModelSequence
        a = DenseAlignment(map(ABModelSequence, ['abab','bbbb','abbb']), \
            Alphabet=AB.Alphabet)
        f = a.getSeqEntropy()
        e = 0.81127812445913283 #sum(p log_2 p) for p = 0.25, 0.75
        self.assertFloatEqual(f, array([1,0,e]))

    def test_getPosEntropy(self):
        """DenseAlignment getPosEntropy should get entropy of each pos"""
        ABModelSequence = self.ABModelSequence
        a = self.a
        f = a.getPosEntropy()
        e = array([0,0,1,1])
        self.assertEqual(f, e)

class IntegrationTests(TestCase):
    """Test for integration between regular and model seqs and alns"""
    def setUp(self):
        """Intialize some standard sequences"""
        self.r1 = RNA.Sequence('AAA', Name='x')
        self.r2 = RNA.Sequence('CCC', Name='y')
        self.m1 = RNA.ModelSeq('AAA', Name='xx')
        self.m2 = RNA.ModelSeq('CCC', Name='yy')
    
    def test_model_to_model(self):
        """Model seq should work with dense alignment"""
        a = DenseAlignment([self.m1, self.m2])
        self.assertEqual(str(a), '>xx\nAAA\n>yy\nCCC\n')
        a = DenseAlignment([self.m1, self.m2], MolType=DNA)
        self.assertEqual(str(a), '>xx\nAAA\n>yy\nCCC\n')
        self.assertEqual(self.m1.Name, 'xx')
    
    def test_regular_to_model(self):
        """Regular seq should work with dense alignment"""
        a = DenseAlignment([self.r1, self.r2])
        self.assertEqual(str(a), '>x\nAAA\n>y\nCCC\n')
        a = DenseAlignment([self.r1, self.r2], MolType=DNA)
        self.assertEqual(str(a), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    
    def test_model_to_regular(self):
        """Model seq should work with regular alignment"""
        a = Alignment([self.m1, self.m2])
        self.assertEqual(str(a), '>xx\nAAA\n>yy\nCCC\n')
        a = Alignment([self.m1, self.m2], MolType=DNA)
        self.assertEqual(str(a), '>xx\nAAA\n>yy\nCCC\n')
        self.assertEqual(self.m1.Name, 'xx')
    
    def test_regular_to_regular(self):
        """Regular seq should work with regular alignment"""
        a = Alignment([self.r1, self.r2])
        self.assertEqual(str(a), '>x\nAAA\n>y\nCCC\n')
        a = Alignment([self.r1, self.r2], MolType=DNA)
        self.assertEqual(str(a), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    
    def test_model_aln_to_regular_aln(self):
        """Dense aln should convert to regular aln"""
        a = DenseAlignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        d = Alignment(a, MolType=DNA)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    
    def test_regular_aln_to_model_aln(self):
        """Regular aln should convert to model aln"""
        a = Alignment([self.r1, self.r2])
        d = DenseAlignment(a)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        d = DenseAlignment(a, MolType=DNA)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    
    def test_regular_aln_to_regular_aln(self):
        """Regular aln should convert to regular aln"""
        a = Alignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        d = Alignment(a, MolType=DNA)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    
    def test_model_aln_to_model_aln(self):
        """Model aln should convert to model aln"""
        a = Alignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        d = Alignment(a, MolType=DNA)
        self.assertEqual(str(d), '>x\nAAA\n>y\nCCC\n')
        self.assertEqual(self.r1.Name, 'x')
    

#run tests if invoked from command line
if __name__ == '__main__':
    main()

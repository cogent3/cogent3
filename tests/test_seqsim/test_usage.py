#!/usr/bin/env python
"""Unit tests for usage and substitution matrices.
"""
from cogent.util.unit_test import TestCase, main
from cogent.core.moltype import RNA
from cogent.core.usage import RnaBases, DnaBases, RnaPairs, DnaPairs
from cogent.core.alphabet import Alphabet
from cogent.core.sequence import ModelRnaSequence as RnaSequence, \
    ModelRnaCodonSequence
from cogent.seqsim.usage import Usage, DnaUsage, RnaUsage, PairMatrix, Counts,\
    Probs, Rates, goldman_q_dna_pair, goldman_q_rna_pair,\
    goldman_q_dna_triple, goldman_q_rna_triple
from numpy import average, asarray, sqrt, identity, diagonal, trace, \
                  array, sum
from cogent.maths.matrix_logarithm import logm
from cogent.maths.matrix_exponentiation import FastExponentiator as expm

#need to find test directory to get access to the tests of the Freqs interface
try:
    from os import getcwd
    from sys import path
    from os.path import sep,join
    test_path = getcwd().split(sep)
    index = test_path.index('tests')
    fields = test_path[:index+1] + ["test_maths"]
    test_path = sep + join(*fields)
    path.append(test_path)
    from test_stats.test_util import StaticFreqsTestsI

    my_alpha = Alphabet('abcde')
    class myUsage(Usage):
        Alphabet = my_alpha

    class UsageAsFreqsTests(StaticFreqsTestsI, TestCase):
        """Note that the remaining Usage methods are tested here."""
        ClassToTest=myUsage

except ValueError:  #couldn't find directory
    pass

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

NUM_TESTS = 10  #for randomization trials

class UsageTests(TestCase):
    """Tests of the Usage object."""
    def setUp(self):
        """Defines some standard test items."""
        self.ab = Alphabet('ab')

        class abUsage(Usage):
            Alphabet = self.ab
            
        self.abUsage = abUsage
        
    def test_init(self):
        """Usage init should succeed only in subclasses"""
        self.assertRaises(TypeError, Usage, [1,2,3,4])
        self.assertEqual(self.abUsage().items(), [('a',0),('b',0)])
        self.assertEqual(self.abUsage([5,6]).items(), [('a',5.0),('b',6.0)])
        #should also construct from seq, if not same length as freqs
        self.assertEqual(self.abUsage([0,0,1,1,1,0,1,1]).items(), \
            [('a',3),('b',5)])

    def test_getitem(self):
        """Usage getitem should get item via alphabet"""
        u = self.abUsage([3,4])
        self.assertEqual(u['a'], 3)
        self.assertEqual(u['b'], 4)

    def test_setitem(self):
        """Usage setitem should set item via alphabet"""
        u = self.abUsage([3,4])
        self.assertEqual(u['a'], 3)
        u['a'] = 10
        self.assertEqual(u['a'], 10)
        u['b'] += 5
        self.assertEqual(u['a'], 10)
        self.assertEqual(u['b'], 9)
   
    def test_str(self):
        """Usage str should print like equivalent list"""
        u = self.abUsage()
        self.assertEqual(str(u), "[('a', 0.0), ('b', 0.0)]")

        u = self.abUsage([1,2.0])
        self.assertEqual(str(u), \
            "[('a', 1.0), ('b', 2.0)]")

    def test_iter(self):
        """Usage iter should iterate over keys"""
        u = self.abUsage([1,2])
        x = tuple(u)
        self.assertEqual(x, ('a', 'b'))
        #should be able to convert to dict via iter
        d = dict(u)
        self.assertEqual(dict(u), {'a':1,'b':2})

    def test_cmp(self):
        """Usage cmp should work as expected"""
        a = self.abUsage([3,4])
        b = self.abUsage([3,2])
        c = self.abUsage([3,4])

        self.assertEqual(a, a)
        self.assertNotEqual(a,b)
        self.assertEqual(a,c)

        self.assertEqual(a==a, True)
        self.assertEqual(a!=a, False)
        self.assertEqual(a==b, False)
        self.assertEqual(a!=b, True)
        self.assertEqual(a==c, True)
        self.assertEqual(a!=c, False)
        self.assertEqual(a==3, False)
        self.assertEqual(a!=3, True)

    def test_add(self):
        """Usage add should add two sets of counts together"""
        u, v = self.abUsage([1,2]), self.abUsage([6,4])
        x = self.abUsage([7,6])
        y = self.abUsage([7,6])
        self.assertEqual(x, y)
        self.assertEqual(x, u+v)
        self.assertEqual(u + v, self.abUsage([7,6]))

    def test_sub(self):
        """Usage sub should subtract one set of counts from the other"""
        u, v = self.abUsage([1,2]), self.abUsage([6,4])
        self.assertEqual(v-u, self.abUsage([5,2]))

    def test_mul(self):
        """Usage mul should multiply usage by a scalar"""
        u = self.abUsage([0,4])
        self.assertEqual(u*3, self.abUsage([0,12]))

    def test_div(self):
        """Usage div should divide usage by scalar (unsafely)"""
        u = self.abUsage([0,4])
        self.assertEqual(u/2, self.abUsage([0,2]))
        self.assertEqual(u/8, self.abUsage([0.0,0.5]))
        #note: don't need to divide by floating point to get fractions
        self.assertEqual(u/8.0, self.abUsage([0.0,0.5]))
    
    def test_scale_sum(self):
        """Usage scale_sum should scale usage to specified sum"""
        u = self.abUsage([1,3])
        self.assertEqual(u.scale_sum(12), self.abUsage([3.0, 9.0]))
        self.assertEqual(u.scale_sum(1), self.abUsage([0.25,0.75]))
        #default is sum to 1
        self.assertEqual(u.scale_sum(), self.abUsage([0.25,0.75]))
    
    def test_scale_max(self):
        """Usage scale_max should scale usage to specified max"""
        u = self.abUsage([1,3])
        self.assertEqual(u.scale_max(12), self.abUsage([4.0, 12.0]))
        self.assertEqual(u.scale_max(1), self.abUsage([1/3.0,1.0]))
        #default is max to 1
        self.assertEqual(u.scale_max(), self.abUsage([1/3.0,1.0]))
    
    def test_probs(self):
        """Usage probs should scale usage to sum to 1"""
        u = self.abUsage([1,3])
        self.assertEqual(u.probs(), self.abUsage([0.25,0.75]))

    def test_randomIndices(self):
        """Usage randomIndices should return correct sequence."""
        d = DnaUsage([0.25, 0.5, 0.1, 0.15])
        s = d.randomIndices(7, [0, 0.49, 1, 0.74, 0.76, 0.86, 0.2])
        self.assertEqual(s, array([0,1,3,1,2,3,0]))

        s = d.randomIndices(10000)
        u, c, a, g = [asarray(s==i, 'int32') for i in [0,1,2,3]]
        assert 2300 < sum(u) < 2700
        assert 4800 < sum(c) < 5200
        assert 800 < sum(a) < 1200
        assert 1300 < sum(g) < 1700
    
    def test_fromSeqData(self):
        """Usage fromSeqData should construct from a sequence object w/ data"""
        class o(object): pass
        s = o()
        s._data = array([0,0,0,1])
        self.assertEqual(self.abUsage.fromSeqData(s), self.abUsage([3,1]))
    
    def test_fromArray(self):
        """Usage fromArray should construct from array holding seq of indices"""
        s = array([0,0,0,1])
        self.assertEqual(self.abUsage.fromArray(s), self.abUsage([3,1]))

    def test_get(self):
        """Usage get should behave like dict"""
        u = self.abUsage([3,4])
        self.assertEqual(u.get('a', 5), 3)
        self.assertEqual(u.get('b', 5), 4)
        self.assertEqual(u.get('x', 5), 5)

    def test_values(self):
        """Usage values should return list of values in alphabet order"""
        u = self.abUsage([3,4])
        self.assertEqual(u.values(), [3,4])

    def test_keys(self):
        """Usage keys should return list of symbols in alphabet order"""
        u = self.abUsage([3,4])
        self.assertEqual(u.keys(), ['a','b'])

    def test_items(self):
        """Usage items should return list of key-value pairs"""
        u = self.abUsage([3,4])
        self.assertEqual(u.items(), [('a',3),('b',4)])

    def test_entropy(self):
        """Usage items should calculate their Shannon entropy"""
        #two equal choices implies one bit of entropy
        u = RnaUsage([1,1,0,0])
        self.assertEqual(u.entropy(), 1)
        u = RnaUsage([10,10,0,0])
        self.assertEqual(u.entropy(), 1)
        #four equal choices implies two bits
        u = RnaUsage([3,3,3,3])
        self.assertEqual(u.entropy(), 2)
        #only one choice -> no entropy
        u = RnaUsage([3,0,0,0])
        self.assertEqual(u.entropy(), 0)
        #empty usage also has no entropy
        u = RnaUsage([0,0,0,0])
        self.assertEqual(u.entropy(), 0)
        #calculated this one by hand
        u = RnaUsage([.5,.3,.1,.1])
        self.assertFloatEqual(u.entropy(),1.6854752972273346) 

class PairMatrixTests(TestCase):
    """Tests of the PairMatrix base class."""
    def setUp(self):
        """Define standard alphabet and matrices for tests."""
        self.ab = Alphabet('ab')
        self.ab_pairs = self.ab*self.ab
        self.empty = PairMatrix([0,0,0,0], self.ab_pairs)
        self.named = PairMatrix([[1,2],[3,4]], self.ab_pairs, 'name')
        
    def test_init(self):
        """PairMatrix init requires data and alphabet"""
        #should only care about number of elements, not shape
        p = PairMatrix([1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8], RnaPairs)
        assert p.Alphabet is RnaPairs
        self.assertEqual(len(p._data), 4)
        self.assertEqual(len(p._data.flat), 16)
        self.assertEqual(p._data[0], array([1,2,3,4]))
        self.assertEqual(p._data[1], array([5,6,7,8]))

    def test_init_bad(self):
        """PairMatrix init should fail if data wrong length"""
        self.assertRaises(ValueError, PairMatrix, [1,2,3,4], RnaPairs)
        #should also require alphabet
        self.assertRaises(TypeError, PairMatrix, [1,2,3,4])

    def test_toMatlab(self):
        """PairMatrix toMatlab should return correct format string"""
        self.assertEqual(self.empty.toMatlab(), "m=[0.0 0.0;\n0.0 0.0];\n")

        self.assertEqual(self.named.toMatlab(), \
        "name=[1.0 2.0;\n3.0 4.0];\n")

    def test_str(self):
        """PairMatrix __str__ should return string correpsonding to data"""
        self.assertEqual(str(self.named), str(self.named._data))

    def test_repr(self):
        """PairMatrix __repr__ should return reconstructable string"""
        self.assertEqual(repr(self.named), \
            'PairMatrix('+ repr(self.named._data) + ',' +\
            repr(self.ab_pairs)+",'name')")

    def test_getitem(self):
        """PairMatrix __getitem__ should translate indices and get from array"""
        n = self.named
        self.assertEqual(n['a'], array([1,2]))
        self.assertEqual(n['b'], array([3,4]))
        self.assertEqual(n['a','a'], 1)
        self.assertEqual(n['a','b'], 2)
        self.assertEqual(n['b','a'], 3)
        self.assertEqual(n['b','b'], 4)

        #WARNING: m[a][b] doesn't work b/c indices not translated!
        #must access as m[a,b] instead.
        try:
            x = n['a']['b']
        except ValueError:
            pass

        #should work even if SubAlphabets not the same
        a = Alphabet('ab')
        x = Alphabet('xyz')
        j = a * x
        m = PairMatrix([1,2,3,4,5,6], j)
        self.assertEqual(m['a','x'], 1)
        self.assertEqual(m['a','y'], 2)
        self.assertEqual(m['a','z'], 3)
        self.assertEqual(m['b','x'], 4)
        self.assertEqual(m['b','y'], 5)
        self.assertEqual(m['b','z'], 6)

        #should work even if SubAlphabets are different types
        a = Alphabet([1,2,3])
        b = Alphabet(['abc', 'xyz'])
        j = a * b
        m = PairMatrix([1,2,3,4,5,6], j)
        self.assertEqual(m[1,'abc'], 1)
        self.assertEqual(m[1,'xyz'], 2)
        self.assertEqual(m[2,'abc'], 3)
        self.assertEqual(m[2,'xyz'], 4)
        self.assertEqual(m[3,'abc'], 5)
        self.assertEqual(m[3,'xyz'], 6)
        self.assertEqual(list(m[2]), [3,4])
        #gives KeyError if single item not present in first level
        self.assertRaises(KeyError, m.__getitem__, 'x')

    def test_empty(self):
        """PairMatrix empty classmethod should produce correct class"""
        p = PairMatrix.empty(self.ab_pairs)
        self.assertEqual(p._data.flat, array([0,0,0,0]))
        self.assertEqual(p._data, array([[0,0],[0,0]]))
        self.assertEqual(p._data.shape, (2,2))

    def test_eq(self):
        """Pairmatrix test for equality should check all elements"""
        p = self.ab_pairs
        a = PairMatrix.empty(p)
        b = PairMatrix.empty(p)
        assert a is not b
        self.assertEqual(a, b)
        c = PairMatrix([1,2,3,4], p)
        d = PairMatrix([1,2,3,4], p)
        assert c is not d
        self.assertEqual(c, d)
        self.assertNotEqual(a, c)
        #Note: still compare equal if alphabets are different
        x = Alphabet('xy')
        x = x*x
        y = PairMatrix([1,2,3,4], x)
        self.assertEqual(y, c)
        #should check all elements, not just first
        c = PairMatrix([1,1,1,1], p)
        d = PairMatrix([1,1,1,4], p)
        assert c is not d
        self.assertNotEqual(c, d)

    def test_ne(self):
        """PairMatrix test for inequality should check all elements"""
        p = self.ab_pairs
        a = PairMatrix.empty(p)
        b = PairMatrix.empty(p)
        c = PairMatrix([1,2,3,4], p)
        d = PairMatrix([1,2,3,4], p)
        assert a != c
        assert a == b
        assert c == d
        
        #Note: still compare equal if alphabets are different
        x = Alphabet('xy')
        x = x*x
        y = PairMatrix([1,2,3,4], x)
        assert y == c
        #should check all elements, not just first
        c = PairMatrix([1,1,1,1], p)
        d = PairMatrix([1,1,1,4], p)
        assert c != d

    def test_iter(self):
        """PairMatrix __iter__ should iterate over rows."""
        p = self.ab_pairs
        c = PairMatrix([1,2,3,4], p)
        l = list(c)
        self.assertEqual(len(l), 2)
        self.assertEqual(list(l[0]), [1,2])
        self.assertEqual(list(l[1]), [3,4])

    def test_len(self):
        """PairMatrix __len__ should return number of rows"""
        p = self.ab_pairs
        c = PairMatrix([1,2,3,4], p)
        self.assertEqual(len(c), 2)
        

class CountsTests(TestCase):
    """Tests of the Counts class, including inferring counts from sequences."""

    def test_toProbs(self):
        """Counts toProbs should return valid prob matrix."""
        c = Counts([1,2,3,4,2,2,2,2,0.2,0.4,0.6,0.8,1,0,0,0], RnaPairs)
        p = c.toProbs()
        assert isinstance(p, Probs)
        self.assertEqual(p, Probs([0.1,0.2,0.3,0.4,0.25,0.25,0.25,0.25, \
            0.1,0.2,0.3,0.4,1.0,0.0,0.0,0.0], RnaPairs))
        self.assertEqual(p['U','U'], 0.1)
        self.assertEqual(p['G','U'], 1.0)
        self.assertEqual(p['G','G'], 0.0)
        
    
    def test_fromPair(self):
        """Counts fromPair should return correct counts."""
        s = Counts.fromPair( RnaSequence('UCCGAUCGAUUAUCGGGUACGUA'), \
                             RnaSequence('GUCGAGUAUAGCGUACGGCUACG'),
                             RnaPairs)
        
        assert isinstance(s, Counts)

        vals = [
            ('U','U',0),('U','C',2.5),('U','A',1),('U','G',2.5),
            ('C','U',2.5),('C','C',1),('C','A',1),('C','G',0.5),
            ('A','U',1),('A','C',1),('A','A',1),('A','G',2),
            ('G','U',2.5),('G','C',0.5),('G','A',2),('G','G',2),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)
        #check that it works for big seqs
        s = Counts.fromPair( RnaSequence('UCAG'*1000), \
                             RnaSequence('UGAG'*1000),
                             RnaPairs)
        
        assert isinstance(s, Counts)

        vals = [
            ('U','U',1000),('U','C',0),('U','A',0),('U','G',0),
            ('C','U',0),('C','C',0),('C','A',0),('C','G',500),
            ('A','U',0),('A','C',0),('A','A',1000),('A','G',0),
            ('G','U',0),('G','C',500),('G','A',0),('G','G',1000),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)

        #check that it works for codon seqs
        s1 = ModelRnaCodonSequence('UUCGCG')
        s2 = ModelRnaCodonSequence('UUUGGG')
        c = Counts.fromPair(s1, s2, RNA.Alphabet.Triples**2)
        self.assertEqual(c._data.sum(), 2)
        self.assertEqual(c._data[0,1], 0.5)
        self.assertEqual(c._data[1,0], 0.5)
        self.assertEqual(c._data[55,63], 0.5)
        self.assertEqual(c._data[63,55], 0.5)


    def test_fromTriple(self):
        """Counts fromTriple should return correct counts."""
        cft = Counts.fromTriple
        rs = RnaSequence
        A, C, G, U = map(rs, 'ACGU')
        #counts if different from both the other groups
        s = cft(A, C, C, RnaPairs)
        assert isinstance(s, Counts)
        self.assertEqual(s['C','A'], 1)
        self.assertEqual(s['A','C'], 0)
        self.assertEqual(s['C','C'], 0)
        #try it with longer sequences
        AAA, CCC = map(rs, ['AAA', 'CCC'])
        s = cft(AAA, CCC, CCC, RnaPairs)
        self.assertEqual(s['C','A'], 3)
        self.assertEqual(s['A','C'], 0)
        #doesn't count if all three differ
        ACG, CGA, GAC = map(rs, ['ACG','CGA','GAC'])
        s = cft(ACG, CGA, GAC, RnaPairs)
        self.assertEqual(s['C','A'], 0)
        self.assertEqual(s['A','C'], 0)
        
        self.assertEqual(s, Counts.empty(RnaPairs))
        #counts as no change if same as other sequence...
        s = cft(AAA, AAA, CCC, RnaPairs)
        self.assertEqual(s['A','A'], 3)
        self.assertEqual(s['A','C'], 0)
        #...or same as the outgroup
        s = cft(AAA, CCC, AAA, RnaPairs)
        self.assertEqual(s['A','A'], 3)
        self.assertEqual(s['A','C'], 0)
        #spot-check a mixed example
        s = cft( \
            rs('AUCGCUAGCAUACGUCA'),
            rs('AAGCUGCGUAGCGCAUA'),
            rs('GCGCAUAUGACGAUAGC'),
            RnaPairs
        )
        vals = [
            ('U','U',1),('U','C',0),('U','A',0),('U','G',0),
            ('C','U',0),('C','C',0),('C','A',0),('C','G',1),
            ('A','U',1),('A','C',0),('A','A',4),('A','G',0),
            ('G','U',0),('G','C',1),('G','A',0),('G','G',1),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)
        #check a long sequence
        s = cft( \
            rs('AUCGCUAGCAUACGUCA'*1000),
            rs('AAGCUGCGUAGCGCAUA'*1000),
            rs('GCGCAUAUGACGAUAGC'*1000),
            RnaPairs
        )
        vals = [
            ('U','U',1000),('U','C',0),('U','A',0),('U','G',0),
            ('C','U',0),('C','C',0),('C','A',0),('C','G',1000),
            ('A','U',1000),('A','C',0),('A','A',4000),('A','G',0),
            ('G','U',0),('G','C',1000),('G','A',0),('G','G',1000),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)

        #check that it works when forced to use both variants of fromTriple
        s = cft( \
            rs('AUCGCUAGCAUACGUCA'*1000),
            rs('AAGCUGCGUAGCGCAUA'*1000),
            rs('GCGCAUAUGACGAUAGC'*1000),
            RnaPairs,
            threshold=0 #forces "large" method
        )
        vals = [
            ('U','U',1000),('U','C',0),('U','A',0),('U','G',0),
            ('C','U',0),('C','C',0),('C','A',0),('C','G',1000),
            ('A','U',1000),('A','C',0),('A','A',4000),('A','G',0),
            ('G','U',0),('G','C',1000),('G','A',0),('G','G',1000),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)
        
        s = cft( \
            rs('AUCGCUAGCAUACGUCA'*1000),
            rs('AAGCUGCGUAGCGCAUA'*1000),
            rs('GCGCAUAUGACGAUAGC'*1000),
            RnaPairs,
            threshold=1e12 #forces "small" method
        )
        vals = [
            ('U','U',1000),('U','C',0),('U','A',0),('U','G',0),
            ('C','U',0),('C','C',0),('C','A',0),('C','G',1000),
            ('A','U',1000),('A','C',0),('A','A',4000),('A','G',0),
            ('G','U',0),('G','C',1000),('G','A',0),('G','G',1000),
        ]
        for i, j, val in vals:
            self.assertFloatEqual(s[i,j], val)


        #check that it works for codon seqs
        s1 = ModelRnaCodonSequence('UUCGCG')
        s2 = ModelRnaCodonSequence('UUUGGG')
        s3 = s2
        c = Counts.fromTriple(s1, s2, s3, RNA.Alphabet.Triples**2)
        self.assertEqual(c._data.sum(), 2)
        self.assertEqual(c._data[0,1], 1)
        self.assertEqual(c._data[63,55], 1)

        
class ProbsTests(TestCase):
    """Tests of the Probs class."""
    def setUp(self):
        """Define an alphabet and some probs."""
        self.ab = Alphabet('ab')
        self.ab_pairs = self.ab**2
        
    def test_isValid(self):
        """Probs isValid should return True if it's a prob matrix"""
        a = self.ab_pairs
        m = Probs([0.5,0.5,1,0], a)
        self.assertEqual(m.isValid(), True)
        #fails if don't sum to 1
        m = Probs([0.5, 0, 1, 0], a)
        self.assertEqual(m.isValid(), False)
        #fails if negative elements
        m = Probs([1, -1, 0, 1], a)
        self.assertEqual(m.isValid(), False)

    def test_makeModel(self):
        """Probs makeModel should return correct substitution pattern"""
        a = Alphabet('abc')**2
        m = Probs([0.5,0.25,0.25,0.1,0.8,0.1,0.3,0.6,0.1], a)
        obs = m.makeModel(array([0,1,1,0,2,2]))
        exp = array([[0.5,0.25,0.25],[0.1,0.8,0.1],[0.1,0.8,0.1],\
            [0.5,0.25,0.25],[0.3,0.6,0.1],[0.3,0.6,0.1]])
        self.assertEqual(obs, exp)

    def test_mutate(self):
        """Probs mutate should return correct vector from input vector"""
        a = Alphabet('abc')**2
        m = Probs([0.5,0.25,0.25,0.1,0.8,0.1,0.3,0.6,0.1], a)
        #because of fp math in accumulate, can't predict boundaries exactly
        #so add/subtract eps to get the result we expect
        eps = 1e-6
        #            a b b a c c a b c
        seq = array([0,1,1,0,2,2,0,1,2])
        random_vec = array([0,.01,.8-eps,1,1,.3,.05,.9+eps,.95])
        self.assertEqual(m.mutate(seq, random_vec), \
            #      a a b c c a a c c
            array([0,0,1,2,2,0,0,2,2]))
        #check that freq. distribution is about right
        seqs = array([m.mutate(seq) for i in range(1000)])
        #WARNING: bool operators return byte arrays, whose sums wrap at 256!
        zero_count = asarray(seqs == 0, 'int32')
        sums = sum(zero_count, axis=0)
        #expect: 500, 100, 100, 500, 300, 300, 500, 100, 300
        #std dev = sqrt(npq), which is sqrt(250), sqrt(90), sqrt(210)
        means = array([500, 100, 100, 500, 300, 300, 500, 100, 300])
        var   = array([250, 90, 90,  250, 210, 210, 250, 90, 210])
        three_sd = 3 * sqrt(var)
        for obs, exp, sd in zip(sums, means, three_sd):
            assert exp - 2*sd < obs < exp + 2*sd

    def test_toCounts(self):
        """Probs toCounts should return counts object w/ right numbers"""
        a = Alphabet('abc')**2
        m = Probs([0.5,0.25,0.25,0.1,0.8,0.1,0.3,0.6,0.1], a)
        obs = m.toCounts(30)
        assert isinstance(obs, Counts)
        exp = Counts([[5.,2.5,2.5,1,8,1,3,6,1]], a)
        self.assertEqual(obs, exp)

    def test_toRates(self):
        """Probs toRates should return log of probs, optionally normalized"""
        a = Alphabet('abc')**2
        p = Probs([0.9,0.05,0.05,0.1,0.85,0.05,0.02,0.02,0.96], a)
        assert p.isValid()
        r = p.toRates()
        assert isinstance(r, Rates)
        assert r.isValid()
        assert not r.isComplex()
        self.assertEqual(r._data, logm(p._data))
        r_norm = p.toRates(normalize=True)
        self.assertFloatEqual(trace(r_norm._data), -1.0)

    def test_random_p_matrix(self):
        """Probs random should return random Probsrows that sum to 1"""
        for i in range(NUM_TESTS):
            p = Probs.random(RnaPairs)._data
            for i in p:
                self.assertFloatEqual(sum(i), 1.0)
            #length should be 4 by default
            self.assertEqual(len(p), 4)
            self.assertEqual(len(p[0]), 4)
            
    def test_random_p_matrix_diag(self):
        """Probs random should work with a scalar diagonal"""
        #if diagonal is 1, off-diagonal elements should be 0
        for i in range(NUM_TESTS):
            p = Probs.random(RnaPairs, 1)._data
            self.assertEqual(p, identity(4, 'd'))
        #if diagonal is between 0 and 1, rows should sum to 1
        for i in range(NUM_TESTS):
            p = Probs.random(RnaPairs, 0.5)._data
            for i in range(4):
                self.assertFloatEqual(sum(p[i]), 1.0)
                self.assertEqual(p[i][i], 0.5)
                assert min(p[i]) >= 0
                assert max(p[i]) <= 1
        #if diagonal > 1, rows should still sum to 1
        for i in range(NUM_TESTS):
            p = Probs.random(RnaPairs, 2)._data
            for i in range(4):
                self.assertEqual(p[i][i], 2.0)
                self.assertFloatEqual(sum(p[i]), 1.0)
                assert min(p[i]) < 0

    def test_random_p_matrix_diag_vector(self):
        """Probs random should work with a vector diagonal"""
        for i in range(NUM_TESTS):
            diag = [0, 0.2, 0.6, 1.0]
            p = Probs.random(RnaPairs, diag)._data
            for i, d, row in zip(range(4), diag, p):
                self.assertFloatEqual(sum(row), 1.0)
                self.assertEqual(row[i], diag[i])

class RatesTests(TestCase):
    """Tests of the Rates class."""
    def setUp(self):
        """Define standard alphabets."""
        self.abc = Alphabet('abc')
        self.abc_pairs = self.abc**2
        
    def test_init(self):
        """Rates init should take additional parameter to normalize"""
        r = Rates([-2,1,1,0,-1,1,0,0,0], self.abc_pairs)
        self.assertEqual(r._data, array([[-2,1,1],[0,-1,1],[0,0,0]]))
       
        r = Rates([-2.5,1,1,0,-1,1,0,0,0], self.abc_pairs)
        self.assertEqual(r._data, array([[-2.5,1.,1.],[0.,-1.,1.],[0.,0.,0.]]))
       
        r = Rates([-2,1,1,0,-1,1,2,0,-1], self.abc_pairs, normalize=True)
        self.assertEqual(r._data, \
            array([[-0.5,.25,.25],[0.,-.25,.25],[.5,0.,-.25]]))

    def test_isComplex(self):
        """Rates isComplex should return True if complex elements"""
        r = Rates([0,0,0.1j,0,0,0,0,0,0], self.abc_pairs)
        assert r.isComplex()
        
        r = Rates([0,0,0.1,0,0,0,0,0,0], self.abc_pairs)
        assert not r.isComplex()

    def test_isSignificantlyComplex(self):
        """Rates isSignificantlyComplex should be true if large imag component"""
        r = Rates([0,0,0.2j,0,0,0,0,0,0], self.abc_pairs)
        assert r.isSignificantlyComplex()
        assert r.isSignificantlyComplex(0.01)
        assert not r.isSignificantlyComplex(0.2)
        assert not r.isSignificantlyComplex(0.3)
        
        r = Rates([0,0,0.1,0,0,0,0,0,0], self.abc_pairs)
        assert not r.isSignificantlyComplex()
        assert not r.isSignificantlyComplex(1e-30)
        assert not r.isSignificantlyComplex(1e3)
        
    def test_isValid(self):
        """Rates isValid should check row sums and neg off-diags"""
        r = Rates([-2,1,1,0,-1,1,0,0,0], self.abc_pairs)
        assert r.isValid()

        r = Rates([0,0,0,0,0,0,0,0,0], self.abc_pairs)
        assert r.isValid()
        #not valid if negative off-diagonal
        r = Rates([-2,-1,3,1,-1,0,2,2,-4], self.abc_pairs)
        assert not r.isValid()
        #not valid if rows don't all sum to 0
        r = Rates([0,0.0001,0,0,0,0,0,0,0], self.abc_pairs)
        assert not r.isValid()

    def test_normalize(self):
        """Rates normalize should return normalized copy of self where trace=-1"""
        r = Rates([-2,1,1,0,-1,1,2,0,-1], self.abc_pairs)
        n = r.normalize()
        self.assertEqual(n._data, \
            array([[-0.5,.25,.25],[0.,-.25,.25],[.5,0.,-.25]]))
        #check that we didn't change the original
        assert n._data is not r._data
        self.assertEqual(r._data, \
            array([[-2,1,1,],[0,-1,1,],[2,0,-1]]))

    def test_toProbs(self):
        """Rates toProbs should return correct probability matrix"""
        a = self.abc_pairs
        p = Probs([0.75, 0.1, 0.15, 0.2, 0.7, 0.1, 0.05, 0.1, 0.85], a)
        q = p.toRates()
        self.assertEqual(q._data, logm(p._data))
        p2 = q.toProbs()
        self.assertFloatEqual(p2._data, p._data)
        
        #test a case that didn't work for DNA
        q = Rates(array(
            [[-0.64098451,  0.0217681 ,  0.35576469,  0.26345171],
             [ 0.31144238, -0.90915091,  0.25825858,  0.33944995],
             [ 0.01578521,  0.43162879, -0.99257581,  0.54516182],
             [ 0.13229986,  0.04027147,  0.05817791, -0.23074925]]),
            DnaPairs)
        self.assertFloatEqual(q.toProbs(0.5)._data, expm(q._data)(t=0.5))

    def test_timeForSimilarity(self):
        """Rates timeToSimilarity should return correct time"""
        a = self.abc_pairs
        p = Probs([0.75, 0.1, 0.15, 0.2, 0.7, 0.1, 0.05, 0.15, 0.8], a)
        q = p.toRates()
        d = 0.5
        t = q.timeForSimilarity(d)
        x = expm(q._data)(t)
        self.assertFloatEqual(average(diagonal(x), axis=0), d)
        t = q.timeForSimilarity(d, array([1/3.0]*3))
        x = expm(q._data)(t)
        self.assertFloatEqual(average(diagonal(x), axis=0), d)
        self.assertEqual(q.timeForSimilarity(1), 0)

    def test_toSimilarProbs(self):
        """Rates toSimilarProbs should match individual steps"""
        a = self.abc_pairs
        p = Probs([0.75, 0.1, 0.15, 0.2, 0.7, 0.1, 0.05, 0.15, 0.8], a)
        q = p.toRates()
        self.assertEqual(q.toSimilarProbs(0.5), \
            q.toProbs(q.timeForSimilarity(0.5)))

        #test a case that didn't work for DNA
        q = Rates(array(
            [[-0.64098451,  0.0217681 ,  0.35576469,  0.26345171],
             [ 0.31144238, -0.90915091,  0.25825858,  0.33944995],
             [ 0.01578521,  0.43162879, -0.99257581,  0.54516182],
             [ 0.13229986,  0.04027147,  0.05817791, -0.23074925]]),
            DnaPairs)
        p = q.toSimilarProbs(0.66)
        self.assertFloatEqual(average(diagonal(p._data), axis=0), 0.66)

    def test_random_q_matrix(self):
        """Rates random should return matrix of correct size"""
        for i in range(NUM_TESTS):
            q = Rates.random(RnaPairs)._data
            self.assertEqual(len(q), 4)
            self.assertEqual(len(q[0]), 4)
            for row in q:
                self.assertFloatEqual(sum(row), 0.0)
                assert min(row) < 0
                assert max(row) > 0
                l = list(row)
                l.sort()
                assert min(l[1:]) >= 0
                assert max(l[1:]) <= 1

    def test_random_q_matrix_diag(self):
        """Rates random should set diagonal correctly from scalar"""
        for i in range(NUM_TESTS):
            q = Rates.random(RnaPairs, -1)._data
            self.assertEqual(len(q), 4)
            for i, row in enumerate(q):
                self.assertFloatEqual(sum(row), 0)
                self.assertEqual(row[i], -1)
                assert max(row) <= 1
                l = list(row)
                l.sort()
                assert min(l[1:]) >= 0
                assert max(l[1:]) <= 1
        for i in range(NUM_TESTS):
            q = Rates.random(RnaPairs, -5)._data
            self.assertEqual(len(q), 4)
            for i, row in enumerate(q):
                self.assertFloatEqual(sum(row), 0)
                self.assertEqual(row[i], -5)
                assert max(row) <= 5
                l = list(row)
                l.sort()
                assert min(l[1:]) >= 0
                assert max(l[1:]) <= 5

    def test_random_q_matrix_diag_vector(self):
        """Rates random should init with vector as diagonal"""
        diag = [1, -1, 2, -2]
        for i in range(NUM_TESTS):
            q = Rates.random(RnaPairs, diag)._data
            for i, d, row in zip(range(4), diag, q):
                self.assertFloatEqual(sum(row, axis=0), 0.0)
                self.assertEqual(row[i], diag[i])

    def test_fixNegsDiag(self):
        """Rates fixNegsDiag should fix negatives by adding to diagonal"""
        q = Rates([[-6,2,2,2],[-6,-2,4,4],[2,2,-6,2],[4,4,-2,-6]], RnaPairs)
        m = q.fixNegsDiag()._data
        self.assertEqual(m,array([[-6,2,2,2],[0,-8,4,4],[2,2,-6,2],[4,4,0,-8]]))

    def test_fixNegsEven(self):
        """Rates fixNegsEven should fix negatives by adding evenly to others"""
        q = Rates([[-6,2,2,2],[-3,-2,3,2],[-2,-2,-6,2],[4,4,-6,-2]], RnaPairs)
        m = q.fixNegsEven()._data
        self.assertEqual(m,array([[-6,2,2,2],[0,-3,2,1],[0,0,-0,0],[2,2,0,-4]]))

    def test_fixNegsFmin(self):
        """Rates fixNegsFmin should fix negatives using fmin method"""
        q = Rates(array([[-0.28936029,  0.14543346, -0.02648614,  0.17041297],
            [ 0.00949624, -0.31186005,  0.17313171,  0.1292321 ],
            [ 0.10443209,  0.16134479, -0.30480186,  0.03902498],
            [ 0.01611264,  0.12999161,  0.15558259, -0.30168684]]), DnaPairs)
        r = q.fixNegsFmin()
        assert not q.isValid()
        assert r.isValid()

    def test_fixNegsConstrainedOpt(self):
        """Rates fixNegsConstrainedOpt should fix negatives w/ constrained opt"""
        q = Rates(array([[-0.28936029,  0.14543346, -0.02648614,  0.17041297],
            [ 0.00949624, -0.31186005,  0.17313171,  0.1292321 ],
            [ 0.10443209,  0.16134479, -0.30480186,  0.03902498],
            [ 0.01611264,  0.12999161,  0.15558259, -0.30168684]]), DnaPairs)
        r = q.fixNegsFmin()
        assert not q.isValid()
        assert r.isValid()

    def test_fixNegsReflect(self):
        """Rates fixNegsReflect should reflect negatives across diagonal"""
        ab = Alphabet('ab')**2
        #should leave matrix alone if no off-diagonal elements
        q = Rates([0,0,1,-1], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[0,0],[1,-1]]))
        q = Rates([-2,2,1,-1], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[-2,2],[1,-1]]))
        #should work if precisely one off-diag element in a pair is negative
        q = Rates([2,-2,1,-1], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[0,0],[3,-3]]))
        q = Rates([-1,1,-2,2], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[-3,3],[0,-0]]))
        #should work if both off-diag elements in a pair are negative
        q = Rates([1,-1,-2,2], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[-2,2],[1,-1]]))
        q = Rates([2,-2,-1,1], ab)
        self.assertEqual(q.fixNegsReflect()._data, array([[-1,1],[2,-2]]))
        
        q = Rates([[ 0,  3, -2, -1],
                   [ 2, -1,  2, -3],
                   [-1, -1,  2,  0],
                   [-3,  2,  0,  1]], RnaPairs)
        q2 = q.fixNegsReflect()
        self.assertEqual(q2._data, \
            array([[-7,  3,  1,  3],
                   [ 2, -5,  3,  0],
                   [ 2,  0, -2,  0],
                   [ 1,  5,  0, -6]]))

class GoldmanTests(TestCase):
    def setUp(self):
        pass

    def test_goldman_q_dna_pair(self):
        """Should return expected rate matrix"""
        seq1 = "ATGCATGCATGC"
        seq2 = "AAATTTGGGCCC"

        expected = array([[-(2/3.0), (1/3.0), (1/3.0), 0],
                          [(1/3.0), -(2/3.0), 0, (1/3.0)],
                          [(1/3.0), 0, -(2/3.0), (1/3.0)],
                          [0, (1/3.0), (1/3.0), -(2/3.0)]])
        observed = goldman_q_dna_pair(seq1, seq2)
        self.assertFloatEqual(observed, expected)

    def test_goldman_q_rna_pair(self):
        """Should return expected rate matrix"""
        seq1 = "AUGCAUGCAUGC"
        seq2 = "AAAUUUGGGCCC"

        expected = array([[-(2/3.0), (1/3.0), (1/3.0), 0],
                          [(1/3.0), -(2/3.0), 0, (1/3.0)],
                          [(1/3.0), 0, -(2/3.0), (1/3.0)],
                          [0, (1/3.0), (1/3.0), -(2/3.0)]])
        observed = goldman_q_rna_pair(seq1, seq2)
        self.assertFloatEqual(observed, expected)

    def test_goldman_q_dna_triple(self):
        """Should return expected rate matrix"""
        seq1 = "ATGCATGCATGC"
        seq2 = "AAATTTGGGCCC"
        outgroup = "AATTGGCCAATT"

        expected = array([[-(1/2.0), (1/2.0), 0, 0],
                          [0, 0, 0, 0],
                          [(1/3.0), 0, -(1/3.0), 0],
                          [0, 0, 0, 0]])
        observed = goldman_q_dna_triple(seq1, seq2, outgroup)
        self.assertFloatEqual(observed, expected)

    def test_goldman_q_rna_triple(self):
        """Should return expected rate matrix"""
        seq1 = "AUGCAUGCAUGC"
        seq2 = "AAAUUUGGGCCC"
        outgroup = "AAUUGGCCAAUU"

        expected = array([[-(1/2.0), (1/2.0), 0, 0],
                          [0, 0, 0, 0],
                          [(1/3.0), 0, -(1/3.0), 0],
                          [0, 0, 0, 0]])
        observed = goldman_q_rna_triple(seq1, seq2, outgroup)
        self.assertFloatEqual(observed, expected)

if __name__ == "__main__":
    main()

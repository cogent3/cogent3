#!/usr/bin/env python
"""test_sequence_generator.py: tests of the sequence_generator module.
"""
from cogent.seqsim.sequence_generators import permutations, combinations, \
    SequenceGenerator, Partition, Composition, \
    MageFrequencies, SequenceHandle, IUPAC_DNA, IUPAC_RNA, BaseFrequency, \
    PairFrequency, BasePairFrequency, RegionModel, ConstantRegion, \
    UnpairedRegion, ShuffledRegion, PairedRegion, MatchingRegion, \
    SequenceModel, Rule, Motif, Module, SequenceEmbedder
from StringIO import StringIO
from operator import mul
from sys import path
from cogent.maths.stats.util import Freqs
from cogent.util.misc import app_path
from cogent.struct.rna2d import ViennaStructure
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

#need to skip some tests if RNAfold absent
if app_path('RNAfold'):
    RNAFOLD_PRESENT = True
else:
    RNAFOLD_PRESENT = False

class FunctionTests(TestCase):
    """Tests of standalone functions"""
    def setUp(self):
        self.standards = (0, 1, 5, 30, 173, 1000, 4382)
        
    def test_permuations_negative_k(self):
        """permutations should raise IndexError if k negative"""
        self.assertRaises(IndexError, permutations, 3, -1)
    
    def test_permutations_k_more_than_n(self):
        """permutations should raise IndexError if k > n"""
        self.assertRaises(IndexError, permutations, 3, 4)
    
    def test_permutations_negative_n(self):
        """permutations should raise IndexError if n negative"""
        self.assertRaises(IndexError, permutations, -3, -2)

    def test_permutations_k_equals_1(self):
        """permutations should return n if k=1"""
        for n in self.standards[1:]:
            self.assertEqual(permutations(n,1), n)

    def test_permutations_k_equals_2(self):
        """permutations should return n*(n-1) if k=2"""
        for n in self.standards[2:]:
            self.assertEqual(permutations(n,2), n*(n-1))

    def test_permutations_k_equals_n(self):
        """permutations should return n! if k=n"""
        for n in self.standards[1:]:
            self.assertEqual(permutations(n,n), reduce(mul, range(1,n+1)))
            
    def test_combinations_k_equals_n(self):
        """combinations should return 1 if k = n"""
        for n in self.standards:
            self.assertEqual(combinations(n,n), 1)
    
    def test_combinations_k_equals_n_minus_1(self):
        """combinations should return n if k=(n-1)"""
        for n in self.standards[1:]:
            self.assertEqual(combinations(n, n-1), n)

    def test_combinations_zero_k(self):
        """combinations should return 1 if k is zero"""
        for n in self.standards:
            self.assertEqual(combinations(n, 0), 1)

    def test_combinations_symmetry(self):
        """combinations(n,k) should equal combinations(n,n-k)"""
        for n in self.standards[3:]:
            for k in (0, 1, 5, 18):
                self.assertEquals(combinations(n, k), combinations(n, n-k))

    def test_combinations_arbitrary_values(self):
        """combinations(n,k) should equal results from spreadsheet"""
        results = {
            30:{0:1, 1:30, 5:142506, 18:86493225, 29:30, 30:1},
            173:{0:1, 1:173, 5:1218218079, 18:1.204353e24, 29:7.524850e32, \
                 30:3.611928e33},
            1000:{0:1, 1:1000, 5:8.2502913e12, 18:1.339124e38,29:7.506513e55, \
                  30:2.429608e57},
            4382:{0:1, 1:4382, 5:1.343350e16, 18:5.352761e49, 29:4.184411e74, \
                  30:6.0715804e76},
                  }
        for n in self.standards[3:]:
            for k in (0, 1, 5, 18, 29, 30):
                self.assertFloatEqualRel(combinations(n,k), results[n][k], 1e-5)
                
class SequenceGeneratorTests(TestCase):
    """Tests of SequenceGenerator, which fills in degenerate bases"""
    def setUp(self):
        """Defines a few standard generators"""
        self.rna_codons = SequenceGenerator('NNN')
        self.dna_iupac_small = SequenceGenerator('RH', IUPAC_DNA)
        self.empty = SequenceGenerator('')
        self.huge = SequenceGenerator('N'*50)
        self.binary = SequenceGenerator('01??01', {'0':'0','1':'1','?':'10'})

    def test_len(self):
        """len(SequenceGenerator) should return number of possible matches"""
        lengths = ((self.rna_codons, 64), (self.dna_iupac_small, 6),
                   (self.empty, 0), (self.binary, 4))
        for item, expected in lengths:
            self.assertEqual(len(item), expected)
        try:
            len(self.huge)
        except OverflowError:
            pass
        else:
            raise AssertionError, "Failed to raise expected OverflowError"

    def test_numPossibilities(self):
        """SequenceGenerator.numPossibilities() should be robust to overflow"""
        lengths = ((self.rna_codons, 64), (self.dna_iupac_small, 6),
                   (self.empty, 0), (self.binary, 4), (self.huge, 4**50))
        for item, expected in lengths:
            self.assertEqual(item.numPossibilities(), expected)

    def test_sequences(self):
        """SequenceGenerator should produce the correct list of sequences"""
        self.assertEqual(list(self.empty), [])
        self.assertEqual(list(self.dna_iupac_small), \
            ['AT','AC','AA','GT','GC','GA'])
        codons = []
        for first in 'UCAG':
            for second in 'UCAG':
                for third in 'UCAG':
                    codons.append(''.join([first, second, third]))
        self.assertEqual(list(self.rna_codons), codons)
        #test that it still works if we call the generator a second time
        self.assertEqual(list(self.rna_codons), codons)
        
    def test_iter(self):
        """SequenceGenerator should act like a list with for..in syntax"""
        as_list = list(self.rna_codons)
        for obs, exp in zip(self.rna_codons, as_list):
            self.assertEqual(obs, exp)

    def test_getitem(self):
        """SequenceGenerator should allow __getitem__ like a list"""
        as_list = list(self.rna_codons)
        for i in range(64):
            self.assertEqual(self.rna_codons[i], as_list[i])

        for i in range(1,65):
            self.assertEqual(self.rna_codons[-i], as_list[-i])

        self.assertEqual(self.huge[-1], 'G'*50)

    def test_getitem_slices(self):
        """SequenceGenerator slicing should work the same as a list"""
        e = list(self.rna_codons)
        o = self.rna_codons
        values = ( 
                    (o[:], e[:]),
                    (o[0:], e[0:]),
                    (o[1:], e[1:]),
                    (o[5:], e[5:]),
                    (o[0:5], e[0:5]),
                    (o[1:5], e[1:5]),
                    (o[5:5], e[5:5]),
                    (o[0:1], e[0:1]),
                    (o[len(o)-1:len(o)], e[len(e)-1:len(e)]),
                    (o[len(o):len(o)], e[len(e):len(e)]),
                )
        testnum = 0
        for obs, exp in values:
            testnum += 1
            self.assertEqual(list(obs), exp)

        big = list(self.huge[1:5])
        self.assertEqual(['U'*49+'C', 'U'*49+'A', 'U'*49+'G', 'U'*48+'CU'], big)

class PartitionTests(TestCase):
    """Tests of the Paritition object."""
    def test_single_partition(self):
        """If number of objects = bins * min, only one way to partition"""
        for num_bins in range(1, 10):
            for occupancy in range(10):
                self.assertEqual(len(Partition(num_bins*occupancy, 
                num_bins, occupancy)), 1)

    def test_partitions(self):
        """Test several properties of partitions, especially start/end"""
        for num_bins in range(1, 5):
            for occupancy in range(5):
                for num_items in \
                range(num_bins*occupancy, num_bins*occupancy + 10):
                    p = Partition(num_items, num_bins, occupancy)
                    l = [i for i in p]
                    l2 = [i for i in p]
                    #check that calling it twice doesn't break it
                    self.assertEqual(l, l2)
                    #check the lengths
                    self.assertEqual(len(p), len(l))
                    #check the ranges are the same...
                    self.assertEqual(l[0][1:], l[-1][0:-1])
                    #and that they contain the right values.
                    self.assertEqual(l[0][1:], [occupancy]*(num_bins - 1))
                    #check the first and last elements
                    self.assertEqual(l[0][0], l[-1][-1])
                    self.assertEqual(l[0][0], \
                    num_items - occupancy * (num_bins - 1))
                    
    def test_values(self):
        """Partition should match precalculated values"""
        self.assertEqual(len(Partition(20, 4, 1)), 969)

    def test_str(self):
        """str(partition) should work as expected"""
        p = Partition(20,4,1)
        self.assertEqual(str(p), "Items: 20 Pieces: 4 Min Per Piece: 1")

        p.NumItems = 13
        p.NumPieces = 2
        p.MinOccupancy = 0
        self.assertEqual(len(p), len(Partition(13, 2, 0)))
        self.assertEqual(str(p), "Items: 13 Pieces: 2 Min Per Piece: 0")

class CompositionTests(TestCase):
    """Tests of the Composition class."""
    def setUp(self):
        """Define a few standard compositions."""
        self.bases_10pct = Composition(10, 0, "ACGU")
        self.bases_5pct = Composition(5, 1, "ACGU")
        self.bases_extra = Composition(10, 0, "CYGEJ")
        self.small = Composition(20, 0, "xy")
        self.unique = Composition(20, 1, "z")

    def test_lengths(self):
        """Composition should return correct number of elements"""
        self.assertEqual(len(self.bases_10pct), len(Partition(10,4,0)))
        self.assertEqual(len(self.bases_5pct), len(Partition(20,4,1)))
        self.assertEqual(len(self.bases_extra), len(Partition(10,5,0)))
        self.assertEqual(len(self.small), len(Partition(5, 2, 0)))
        self.assertEqual(len(self.unique), len(Partition(5, 1, 1)))

    def test_known_vals(self):
        """Composition should return precalculated elements for known cases"""
        self.assertEqual(len(Composition(5,1,"ACGU")), 969) 
        self.assertEqual(len(Composition(5,0,"ACGU")), 1771)
        as_list = list(Composition(5,1,"ACGU"))
        self.assertEqual(as_list[0], Freqs('A'*17+'CGU'))
        self.assertEqual(as_list[-1], Freqs('U'*17+'ACG'))

    def test_updating(self):
        """Composition updates should reset frequencies correctly."""
        exp_list = list(Composition(5, 1, "GCAUN"))
        self.bases_10pct.Spacing = 5
        self.bases_10pct.Alphabet = "GCAUN"
        self.bases_10pct.MinOccupancy = 1
        self.assertEqual(list(self.bases_10pct), exp_list)

class MageFrequenciesTests(TestCase):
    """Tests of the MageFrequencies class -- presentation for Composition."""
    def setUp(self):
        """Define a few standard compositions."""
        self.bases_10pct = Composition(10, 0, "ACGU")
    def test_str(self):
        """MageFrequencies string conversions work correctly"""
        obs_list = list(self.bases_10pct)
        self.assertEqual(str(MageFrequencies(obs_list[0])), '1.0 0.0 0.0')
        self.assertEqual(str(MageFrequencies(obs_list[-1], "last")), \
            '{last} 0.0 0.0 0.0')
        self.assertEqual(str(MageFrequencies({'C':2, 'A':3, 'T':5, 'x':17}, \
        'bases')), '{bases} 0.3 0.2 0.0')

class SequenceHandleTests(TestCase):
    """Tests of the SequenceHandle class."""
    def setUp(self):
        """Define some standard SequenceHandles."""
        self.rna = SequenceHandle('uuca', 'ucag')
        self.any = SequenceHandle(['u', 1, None])
        self.empty = SequenceHandle()
        
    def test_init_good(self):
        """SequenceHandle should init OK without alphabet"""
        self.assertEqual(SequenceHandle('abc123'), list('abc123'))
        self.assertEqual(SequenceHandle(), list())
        self.assertEqual(SequenceHandle('abcaaa', 'abcd'), list('abcaaa'))
        self.assertEqual(SequenceHandle([1,2,3]), [1,2,3])

    def test_init_bad(self):
        """SequenceHandle should raise ValueError if item not in alphabet"""
        self.assertRaises(ValueError, SequenceHandle, 'abc1', 'abc')
        self.assertRaises(ValueError, SequenceHandle, '1', [1])

    def test_setitem_good(self):
        """SequenceHandle setitem should allow items in alphabet"""
        self.rna[0] = 'c'
        self.assertEqual(self.rna, list('cuca'))
        self.rna[-1] = 'u'
        self.assertEqual(self.rna, list('cucu'))
        self.any[1] = [1, 2, 3]
        self.assertEqual(self.any, ['u', [1, 2, 3], None])

    def test_setitem_bad(self):
        """SequenceHandle setitem should reject items not in alphabet"""
        self.assertRaises(ValueError, self.rna.__setitem__, 0, 'x')

    def test_setslice_good(self):
        """SequenceHandle setslice should allow same-length slice"""
        self.rna[:] = list('aaaa')
        self.assertEqual(self.rna, list('aaaa'))
        self.rna[0:1] = ['u']
        self.assertEqual(self.rna, list('uaaa'))
        self.rna[-2:] = ['g','g']
        self.assertEqual(self.rna, list('uagg'))

    def test_setslice_bad(self):
        """SequenceHandle setslice should reject bad items or length change"""
        self.assertRaises(ValueError, self.rna.__setslice__, 0, len(self.rna), \
            ['a']*5)
        self.assertRaises(ValueError, self.any.__setslice__, 0, len(self.any), \
            ['a']*5)
        self.assertRaises(ValueError, self.rna.__setslice__, 0, 1, ['x'])

    def test_string(self):
        """SequenceHandle str should join items without spaces"""
        #use ''.join if items are strings
        self.assertEqual(str(self.rna), 'uuca')
        self.assertEqual(str(self.empty), '')
        #if some of the items raise errors, use built-in method instead
        self.assertEqual(str(self.any), str(['u', 1, None]))

    def test_naughty_methods(self):
        """SequenceHandle list mutators should raise NotImplementedError"""
        r = self.rna
        naughty = [r.__delitem__, r.__delslice__, r.__iadd__, r.__imul__, \
                   r.append, r.extend, r.insert, r.pop, r.remove]
        for n in naughty:
            self.assertRaises(NotImplementedError, n)

class BaseFrequencyTests(TestCase):
    """Tests of BaseFrequency class: wrapper for FrequencyDistibution."""
    def test_init(self):
        """BaseFrequency should init as expected"""
        self.assertEqual(BaseFrequency('UUUCCCCAG'), \
                         Freqs('UUUCCCCAG', 'UCAG'))
        self.assertEqual(BaseFrequency('TTTCAGG', RNA=False), \
                         Freqs('TTTCAGG'))
    def test_init_bad(self):
        """BaseFrequency init should disallow bad characters"""
        self.assertRaises(Exception, BaseFrequency, 'TTTCAGG')
        self.assertRaises(Exception, BaseFrequency, 'UACGUA', False)

class PairFrequencyTests(TestCase):
    """Tests of PairFrequency class: wrapper for Freqs."""
    def test_init_one_parameter(self):
        """PairFrequency should interpret single parameter as pair probs"""
        obs = PairFrequency('UCCC')
        exp = Freqs({('U','U'):0.0625, ('U','C'):0.1875, 
                      ('C','U'):0.1875, ('C','C'):0.5625})
        
        for k, v in exp.items():
            self.assertEqual(v, obs[k])
        for k, v in obs.items():
            if k not in exp:
                self.assertEqual(v, 0)

        self.assertEqual(PairFrequency('UCCC', [('U','U'),('C','C')]), \
            Freqs({('U','U'):0.1, ('C','C'):0.9}))
        #check that the alphabets are right: should not raise error on
        #incrementing characters already there, but should raise KeyError
        #on anything that's missing.
        p = PairFrequency('UCCC')
        p[('U','U')] += 1
        try:
            p[('X','U')] += 1
        except KeyError:
            pass
        else:
            raise AssertionError, "Expected KeyError."
        p = PairFrequency('UCCC', (('C','C'),))
        p[('C','C')] += 1
        try:
            p[('U','U')] += 1
        except KeyError:
            pass
        else:
            raise AssertionError, "Expected KeyError."

class BasePairFrequencyTests(TestCase):
    """Tests of the BaseFrequency class, constructed for easy initialization."""
    def test_init(self):
        """BaseFrequency init should provide correct PairFrequency"""
        WatsonCrick = [('A','U'), ('U','A'),('G','C'),('C','G')]
        Wobble = WatsonCrick + [('G','U'), ('U','G')]
        #by default, basepair should have the wobble alphabet
        bpf = BasePairFrequency('UUACG')
        pf = PairFrequency('UUACG', Wobble)
        self.assertEqual(bpf, pf)
        self.assertEqual(bpf.Constraint, pf.Constraint)
        #can turn GU off, leading to watson-crickery
        bpf = BasePairFrequency('UUACG', False)
        #make sure this gives different results...
        self.assertNotEqual(bpf, pf)
        self.assertNotEqual(bpf.Constraint, pf.Constraint)
        #...but that the results are the same when the correct alphabet is used
        pf = PairFrequency('UUACG', WatsonCrick)
        self.assertEqual(bpf, pf)
        self.assertEqual(bpf.Constraint, pf.Constraint)

class RegionModelTests(TestCase):
    """Tests of the RegionModel class. Base class just returns the template."""
    def test_init(self):
        """RegionModel base class should always return current template."""
        #test blank region model
        r = RegionModel()
        self.assertEqual(str(r.Current), '')
        self.assertEqual(len(r), 0)
        #now assign it to a template
        r.Template = ('ACGUUCGA')
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))
        #check that refresh doesn't break anything
        r.refresh()
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))
        #check composition
        self.assertEqual(r.Composition, None)
        d = {'A':3, 'U':10}
        r.Composition = Freqs(d)
        self.assertEqual(r.Composition, d)
        #check that composition doesn't break the update
        r.refresh()
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))
 
class ConstantRegionTests(TestCase):
    """Tests of the ConstantRegion class. Just returns the template."""
    def test_init(self):
        """ConstantRegion should always return current template."""
        #test blank region model
        r = ConstantRegion()
        self.assertEqual(str(r.Current), '')
        self.assertEqual(len(r), 0)
        #now assign it to a template
        r.Template = ('ACGUUCGA')
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))
        #check that refresh doesn't break anything
        r.refresh()
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))
        #check composition
        self.assertEqual(r.Composition, None)
        d = {'A':3, 'U':10}
        r.Composition = Freqs(d)
        self.assertEqual(r.Composition, d)
        #check that composition doesn't break the update
        r.refresh()
        self.assertEqual(str(r.Current), 'ACGUUCGA')
        self.assertEqual(len(r), len('ACGUUCGA'))

class UnpairedRegionTests(TestCase):
    """Tests of unpaired region: should fill in w/ single-base frequencies."""
    def test_init(self):
        """Unpaired region should generate right freqs, even after change"""
        freqs = Freqs({'C':10,'U':1, 'A':0})
        r = UnpairedRegion('NN', freqs)
        seq = r.Current
        assert seq[0] in 'CU'
        assert seq[1] in 'CU'
        self.assertEqual(len(seq), 2)
        fd = []
        for i in range(1000):
            r.refresh()
            fd.append(str(seq))
        fd = Freqs(''.join(fd))

        observed = [fd['C'], fd['U']]
        expected = [1800, 200]
        self.assertSimilarFreqs(observed, expected)
        self.assertEqual(fd['U'] + fd['C'], 2000)

        freqs2 = Freqs({'A':5, 'U':5})
        r.Composition = freqs2
        r.Template = 'NNN'  #note that changing the Template changes seq ref
        seq = r.Current
        self.assertEqual(len(seq), 3)
        assert seq[0] in 'AU'
        assert seq[1] in 'AU'
        assert seq[2] in 'AU'
        fd = []
        for i in range(1000):
            r.refresh()
            fd.append(str(seq))
        fd = Freqs(''.join(fd))
        observed = [fd['A'], fd['U']]
        expected = [1500, 1500]
        self.assertSimilarFreqs(observed, expected)
        self.assertEqual(fd['A'] + fd['U'], 3000)

class ShuffledRegionTests(TestCase):
    """Shuffled region should randomize string"""
    def test_init(self):
        """Shuffled region should init ok with string, ignoring base freqs"""
        #general strategy: seqs should be different, but sorted seqs should
        #be the same
        empty = ''
        seq = 'UUUCCCCAAAGGG'
        #check that we don't get errors on empty template
        r = ShuffledRegion(empty)
        r.refresh()
        self.assertEqual(str(r.Current), '')
        #check that changing the template changes the sequence
        r.Template = seq
        self.assertNotEqual(str(r.Current), '')
        #check that it shuffled the sequence the first time
        self.assertNotEqual(str(r.Current), seq)
        curr = str(r.Current)
        as_list = list(curr)
        #check that we have the right number of each type of base
        as_list.sort()
        exp_as_list = list(seq)
        exp_as_list.sort()
        self.assertEqual(as_list, exp_as_list)
        #check that we get something different if we refresh again
        r.refresh()
        self.assertNotEqual(str(r.Current), curr)
        as_list = list(str(r.Current))
        as_list.sort()
        self.assertEqual(as_list, exp_as_list)

class PairedRegionTests(TestCase):
    """Tests of paired region generation."""
    def test_init(self):
        """Paired region init and mutation should give expected results"""
        WatsonCrick = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
        Wobble = {'A':'U', 'U':'AG', 'C':'G', 'G':'UC'}
        #check that empty init doesn't give errors
        r = PairedRegion()
        r.refresh()
        #check that mutation works correctly
        r.Template = "N"
        self.assertEqual(len(r), 1)
        r.monomers('UCCGGA')
        upstream = r.Current[0]
        downstream = r.Current[1]
        states = {}
        num_to_do = 10000
        for i in range(num_to_do):
            r.refresh()
            curr = (upstream[0], downstream[0])
            assert upstream[0] in Wobble[downstream[0]]
            states[curr] = states.get(curr, 0) + 1
        for i in states.keys():
            assert i[1] in Wobble[i[0]]
        for i in Wobble:
            for j in Wobble[i]:
                assert (i, j) in states.keys()
        expected_dict = {('A','U'):num_to_do/14, ('U','A'):num_to_do/14,
                        ('C','G'):num_to_do/14*4, ('G','C'):num_to_do/14*4,
                        ('U','G'):num_to_do/14*2, ('G','U'):num_to_do/14*2,}
        # the following for loop was replaced with the assertSimilarFreqs
        # call below it
        #for key, val in expected.items():
            #self.assertFloatEqualAbs(val, states[key], 130) #conservative?
        expected = [val for key, val in expected_dict.items()]
        observed = [states[key] for key, val in expected_dict.items()]
        self.assertSimilarFreqs(observed, expected)

        assert ('G','U') in states
        assert ('U','G') in states
            
        r.monomers('UCGA', GU=False)
        upstream = r.Current[0]
        downstream = r.Current[1]
        states = {}
        num_to_do = 10000
        for i in range(num_to_do):
            r.refresh()
            curr = (upstream[0], downstream[0])
            assert upstream[0] in WatsonCrick[downstream[0]]
            states[curr] = states.get(curr, 0) + 1
        for i in states.keys():
            assert i[1] in WatsonCrick[i[0]]
        for i in WatsonCrick:
            for j in WatsonCrick[i]:
                assert (i, j) in states.keys()
        expected_dict = {('A','U'):num_to_do/4, ('U','A'):num_to_do/4,
                        ('C','G'):num_to_do/4, ('G','C'):num_to_do/4,}
        expected = [val for key, val in expected_dict.items()]
        observed = [states[key] for key, val in expected_dict.items()]
        self.assertSimilarFreqs(observed, expected)
        #for key, val in expected.items():
        #    self.assertFloatEqualAbs(val, states[key], 130) #3 std devs
        assert ('G','U') not in states
        assert ('U','G') not in states

class SequenceModelTests(TestCase):
    """Tests of the SequenceModel class."""
    def test_init(self):
        """SequenceModel should init OK with Isoleucine motif."""
        helices = [PairedRegion('NNN'), PairedRegion('NNNNN')]
        constants = [ConstantRegion('CUAC'), ConstantRegion('UAUUGGGG')]
        order = "H0 C0 H1 - H1 C1 H0"
        isoleucine = SequenceModel(order=order, constants=constants, \
            helices=helices)
        isoleucine.Composition = BaseFrequency('UCAG')
        #print
        #print
        for i in range(10):
            isoleucine.refresh()
            #print list(isoleucine)

        #print
        isoleucine.Composition = BaseFrequency('UCAG')
        isoleucine.GU = False
        #print
        for i in range(10):
            isoleucine.refresh()
            #print list(isoleucine)
        #print

class RuleTests(TestCase):
    """Tests of the Rule class"""
    
    def test_init_bad_params(self):
        """Rule should fail validation except with exactly 5 parameters"""
        self.assertRaises(TypeError, Rule, 1, 1, 1, 1)
        self.assertRaises(TypeError, Rule, 1, 1, 1, 1, 1, 1)
        
    def test_init_bad_length(self):
        """Rule should fail validation if helix extends past downstream start"""
        self.assertRaises(ValueError, Rule, 0, 0, 1, 0, 2)
        self.assertRaises(ValueError, Rule, 0, 0, 10, 10, 12)

    def test_init_bad_negative_params(self):
        """Rule should fail validation if any parameters are negative"""
        self.assertRaises(ValueError, Rule, -1, 0, 1, 0, 1)
        self.assertRaises(ValueError, Rule, 0, -1, 1, 1, 1)
        self.assertRaises(ValueError, Rule, 0, 0, -1, 0, 5)
        self.assertRaises(ValueError, Rule, 0, 0, 0, -1, 1)
        self.assertRaises(ValueError, Rule, 0, 0, 1, 1, -1)

    def test_init_bad_zero_length(self):
        """Rule should fail validation if length is zero"""
        self.assertRaises(ValueError, Rule, 0, 0, 1, 1, 0)

    def test_init_overlap(self):
        """Rule should fail validation if bases must pair with themselves"""
        self.assertRaises(ValueError, Rule, 0, 0, 0, 0, 1)
        self.assertRaises(ValueError, Rule, 0, 10, 0, 15, 4)

    def test_init_wrong_order(self):
        """First sequence must have lower index"""
        self.assertRaises(ValueError, Rule, 1, 0, 0, 5, 3)

    def test_init_ok_length(self):
        """Rule should init OK if helix extends to exactly downstream start"""
        x = Rule(0, 0, 1, 0, 1)
        self.assertEqual(str(x), \
        "Up Seq: 0 Up Pos: 0 Down Seq: 1 Down Pos: 0 Length: 1")
        #check adjacent bases
        x = Rule(0, 0, 0, 1, 1)
        self.assertEqual(str(x), \
        "Up Seq: 0 Up Pos: 0 Down Seq: 0 Down Pos: 1 Length: 1")
        x = Rule(1, 10, 2, 8, 7)
        #check rule that would cause overlap if motifs weren't different 
        self.assertEqual(str(x), \
        "Up Seq: 1 Up Pos: 10 Down Seq: 2 Down Pos: 8 Length: 7")

    def test_str(self):
        """Rule str method should give expected results"""
        x = Rule(1, 10, 2, 8, 7)
        self.assertEqual(str(x), \
        "Up Seq: 1 Up Pos: 10 Down Seq: 2 Down Pos: 8 Length: 7")

class RuleTests_compatibility(TestCase):
    """Tests to see whether the Rule compatibility code works"""

    def setUp(self):
        """Sets up some standard rules"""
        self.x = Rule(1, 5, 2, 10, 3)
        self.x_ok = Rule(1, 8, 2, 14, 4)
        self.x_ok_diff_sequences = Rule(3, 5, 5, 10, 3)
        self.x_bad_first = Rule(1, 0, 3, 10, 10)
        self.x_bad_first_2 = Rule(0, 0, 1, 8, 2)
        self.x_bad_second = Rule(1, 15, 2, 15, 8)
        self.x_bad_second_2 = Rule(1, 14, 2, 8, 4)

    def test_is_compatible_ok(self):
        """Rule.isCompatible should return True if rules don't overlap"""
        self.assertEqual(self.x.isCompatible(self.x_ok), True)  #no return value
        self.assertEqual(self.x.isCompatible(self.x_ok_diff_sequences), True)
        #check that it's transitive
        self.assertEqual(self.x_ok.isCompatible(self.x), True)
        self.assertEqual(self.x_ok_diff_sequences.isCompatible(self.x), True)
        
    def test_is_compatible_bad(self):
        """Rule.isComaptible should return False if rules overlap"""
        tests = [   (self.x, self.x_bad_first),
                    (self.x, self.x_bad_first_2),
                    (self.x, self.x_bad_second),
                    (self.x, self.x_bad_second_2),
                ]
        for first, second in tests:
            self.assertEqual(first.isCompatible(second), False)
            #check that it's transitive
            self.assertEqual(second.isCompatible(first), False)

    def test_fits_in_sequence(self):
        """Rule.fitsInSequence should return True if sequence long enough"""
        sequences = map('x'.__mul__, range(21))     #0 to 20 copies of 'x'
        rules = [self.x, self.x_ok, self.x_ok_diff_sequences, self.x_bad_first,
                 self.x_bad_first_2, self.x_bad_second, self.x_bad_second_2]
        #test a bunch of values for all the rules we have handy
        for s in sequences:
            for r in rules:
                if r.UpstreamPosition + r.Length > len(s):
                    self.assertEqual(r.fitsInSequence(s), False)
                else:
                    self.assertEqual(r.fitsInSequence(s), True)
        #test a couple of specific boundary cases
        #length-1 helix
        r = Rule(0, 0, 1, 0, 1)
        self.assertEqual(r.fitsInSequence(''), False)
        self.assertEqual(r.fitsInSequence('x'), True)
        self.assertEqual(r.fitsInSequence('xx'), True)
        #length-2 helix starting one base from the start
        r = Rule(1, 1, 2, 2, 2)
        self.assertEqual(r.fitsInSequence(''), False)
        self.assertEqual(r.fitsInSequence('x'), False)
        self.assertEqual(r.fitsInSequence('xx'), False)
        self.assertEqual(r.fitsInSequence('xxx'), True)
        self.assertEqual(r.fitsInSequence('xxxx'), True)

class ModuleTests(TestCase):
    """Tests of the Module class, which holds sequences and structures."""
    def test_init_bad(self):
        """Module init should fail if seq/struct missing, or mismatched lengths"""
        #test incorrect param number
        self.assertRaises(TypeError, Module, 'abc')
        self.assertRaises(TypeError, Module, 'abc', 'def', 'ghi')
        #test incorrect lengths
        self.assertRaises(ValueError, Module, 'abc', 'abcd')
        self.assertRaises(ValueError, Module, 'abcd', 'acb')

    def test_init_good(self):
        """Module init should work if seq and struct same length"""
        m = Module('U', '.')
        self.assertEqual(m.Sequence, 'U')
        self.assertEqual(m.Structure, '.')
        m.Sequence = ''
        m.Structure = ''
        self.assertEqual(m.Sequence, '')
        self.assertEqual(m.Structure, '')
        m.Sequence = 'CCUAGG'
        m.Structure = '((..))'
        self.assertEqual(m.Sequence, 'CCUAGG')
        self.assertEqual(m.Structure, '((..))')
        m.Structure = ''
        self.assertRaises(ValueError, m.__len__)

    def test_len(self):
        """Module len should work if seq and struct same length"""
        m = Module('CUAG', '....')
        self.assertEqual(len(m), 4)
        m = Module('', '')
        self.assertEqual(len(m), 0)
        m.Sequence = 'AUCGAUCGA'
        self.assertRaises(ValueError, m.__len__)
   
    def test_str(self):
       """Module str should contain sequence and structure"""
       m = Module('CUAG', '....')
       self.assertEqual(str(m), 'Sequence:  CUAG\nStructure: ....')
       m = Module('', '')
       self.assertEqual(str(m), 'Sequence:  \nStructure: ')
   
    def test_matches(self):
        """Module matches should return correct result for seq/struct match"""
        empty = Module('', '')
        short_p = Module('AC', '((')
        short_u = Module('UU', '..')
        short_up = Module('UU', '((')
        long_all = Module('GGGACGGUUGGUUGGUU', ')))((..((....((((') #struct+seq
        long_seq = Module('GGGACGGUUGGUU', ')))))))))))))') #seq but not struct
        long_struct = Module('GGGGGGGGGGGGG', ')))((..((....') #struct, not seq
        long_none = Module('GGGGGGGGGGGGG', ')))))))))))))') #not struct or seq

        #test overall matching
        for matcher in [empty, short_p, short_u, short_up]:
            self.assertEqual(matcher.matches(long_all), True)
            for longer in [long_seq, long_struct, long_none]:
                if matcher is empty:
                    self.assertEqual(matcher.matches(longer), True)
                else:
                    self.assertEqual(matcher.matches(longer), False)
        #test specific positions
        positions = {3:short_p, 11:short_u, 7:short_up, 15:short_up}
        for module in [short_p, short_u, short_up]:
            for i in range(len(long_all)):
                result = module.matches(long_all, i)
                if positions.get(i, None) is module:
                    self.assertEqual(result, True)
                else:
                    self.assertEqual(result, False)

class MotifTests(TestCase):
    """Tests of the Motif object, which has a set of Modules and Rules."""
    def setUp(self):
        """Defines a few standard motifs"""
        self.ile_mod_0 = Module('NNNCUACNNNNN', '(((((..(((((')
        self.ile_mod_1 = Module('NNNNNUAUUGGGGNNN', ')))))......)))))')
        self.ile_rule_0 = Rule(0, 0, 1, 15, 3)
        self.ile_rule_1 = Rule(0, 7, 1, 4, 5)
        self.ile = Motif([self.ile_mod_0, self.ile_mod_1], \
            [self.ile_rule_0, self.ile_rule_1])
        self.hh_mod_0 = Module('NNNNUNNNNN', '(((((.((((')
        self.hh_mod_1 = Module('NNNNCUGANGAGNNN', ')))).......((((')
        self.hh_mod_2 = Module('NNNCGAAANNNN', '))))...)))))')
        self.hh_rule_0 = Rule(0, 0, 2, 11, 5)
        self.hh_rule_1 = Rule(0, 6, 1, 3, 4)
        self.hh_rule_2 = Rule(1, 11, 2, 3, 4)
        self.hh = Motif([self.hh_mod_0, self.hh_mod_1, self.hh_mod_2], \
                        [self.hh_rule_0, self.hh_rule_1, self.hh_rule_2])
        self.simple_0 = Module('CCCCC', '(((..')
        self.simple_1 = Module('GGGGG', '..)))')
        self.simple_r = Rule(0, 0, 1, 4, 3)
        self.simple = Motif([self.simple_0, self.simple_1], [self.simple_r])

    def test_init_bad_rule_lengths(self):
        """Motif init should fail if rules don't match module lengths"""
        bad_rule = Rule(0, 0, 1, 8, 6)
        self.assertRaises(ValueError, Motif, [self.simple_0, self.simple_1], \
            [bad_rule])

    def test_init_conflicting_rules(self):
        """Motif init should fail if rules overlap"""
        interferer = Rule(0, 2, 2, 20, 4)
        self.assertRaises(ValueError, Motif, [self.ile_mod_0, self.ile_mod_1, \
            self.ile_mod_0], [self.ile_rule_0, interferer]) 
    
    def test_matches_simple(self):
        """Test of simple match should work correctly"""
        index =                    '01234567890123456789012345678901'
        seq =                      'AAACCCCCUUUGGGGGAAACCCCCUUUGGGGG'
        struct =   ViennaStructure('((..((..))....))...(((.......)))')   
        struct_2 = ViennaStructure('((((((..((())))))))).....(((.)))')
                                    #substring right, not pair

        self.assertEqual(self.simple.matches(seq, struct, [19, 27]), True)
        self.assertEqual(self.simple.matches(seq, struct_2, [19,27]), False)

        for first_pos in range(len(seq) - len(self.simple_0) + 1):
            for second_pos in range(len(seq) - len(self.simple_1) + 1):
                #should match struct only at one location
                match=self.simple.matches(seq, struct, [first_pos, second_pos])
                if (first_pos == 19) and (second_pos == 27):
                    self.assertEqual(match, True)
                else:
                    self.assertEqual(match, False)
                #should never match in struct_2
                self.assertEqual(self.simple.matches(seq, struct_2, \
                    [first_pos, second_pos]), False)

        #check that it doesn't fail if there are _two_ matches
        index =  '01234567890123456789'
        seq =    'CCCCCGGGGGCCCCCGGGGG'
        struct = '(((....)))(((....)))'
        struct = ViennaStructure(struct)
        self.assertEqual(self.simple.matches(seq, struct, [0, 5]), True)
        self.assertEqual(self.simple.matches(seq, struct, [10,15]), True)
        #not allowed to cross-pair
        self.assertEqual(self.simple.matches(seq, struct, [0, 15]), False)

    def test_matches_ile(self):
        """Test of isoleucine match should work correctly"""
        index =    '012345678901234567890123456789012345'
        seq_good = 'AAACCCCUACUUUUUCCCAAAAAUAUUGGGGGGGAA'
        seq_bad =  'AAACCCCUACUUUUUCCCAAAAAUAUUGGGCGGGAA'
        st_good =  '...(((((..(((((...)))))......)))))..'
        st_bad =   '((((((((..(((((...)))))...))))))))..'

        st_good = ViennaStructure(st_good)
        st_bad = ViennaStructure(st_bad)

        for first_pos in range(len(seq_good) - len(self.ile_mod_0) + 1):
            for second_pos in range(len(seq_good) - len(self.ile_mod_1) + 1):
                #seq_good and struct_good should match at one location
                match=self.ile.matches(seq_good,st_good,[first_pos,second_pos])
                if (first_pos == 3) and (second_pos == 18):
                    self.assertEqual(match, True)
                else:
                    self.assertEqual(match, False)
                self.assertEqual(self.ile.matches(seq_good, st_bad, \
                    [first_pos, second_pos]), False)
                self.assertEqual(self.ile.matches(seq_bad, st_good, \
                    [first_pos, second_pos]), False)
                self.assertEqual(self.ile.matches(seq_bad, st_bad, \
                    [first_pos, second_pos]), False)

    def test_matches_hh(self):
        """Test of hammerhead match should work correctly"""
        index =    '0123456789012345678901234567890123456'
        seq_good = 'CCCCUAGGGGCCCCCUGAAGAGAAAUUUCGAAAGGGG'
        seq_bad ='CCCCCAGGGGCCCCCUGAAGAGAAAUUUCGAAGGGGG'
        structure ='(((((.(((()))).......(((())))...)))))'
        struct = ViennaStructure(structure)
        self.assertEqual(self.hh.matches(seq_good, struct, [0, 10, 25]), True)
        self.assertEqual(self.hh.matches(seq_bad, struct, [0, 10, 25]), False)

    def test_structureMatches_hh(self):
        """Test of hammerhead structureMatch should work correctly"""
        index =    '0123456789012345678901234567890123456'
        seq_good = 'CCCCUAGGGGCCCCCUGAAGAGAAAUUUCGAAAGGGG'
        seq_bad ='CCCCCAGGGGCCCCCUGAAGAGAAAUUUCGAAGGGGG'
        structure ='(((((.(((()))).......(((())))...)))))'
        struct = ViennaStructure(structure)
        self.assertEqual(self.hh.structureMatches(struct, [0, 10, 25]), True)
        self.assertEqual(self.hh.structureMatches(struct, [0, 10, 25]), True)

class SequenceEmbedderTests(TestCase):
    """Tests of the SequenceEmbedder class."""
    def setUp(self):
        """Define a few standard models and motifs"""
        ile_mod_0 = Module('NNNCUACNNNNN', '(((((..(((((')
        ile_mod_1 = Module('NNNNNUAUUGGGGNNN', ')))))......)))))')
        ile_rule_0 = Rule(0, 0, 1, 15, 5)
        ile_rule_1 = Rule(0, 7, 1, 4, 5)
        ile_motif = Motif([ile_mod_0, ile_mod_1], \
            [ile_rule_0, ile_rule_1])

        helices = [PairedRegion('NNN'), PairedRegion('NNNNN')]
        constants = [ConstantRegion('CUAC'), ConstantRegion('UAUUGGGG')]
        order = "H0 C0 H1 - H1 C1 H0"
        ile_model = SequenceModel(order=order, constants=constants, \
            helices=helices, composition=BaseFrequency('UCAG'))

        self.ile_embedder = SequenceEmbedder(length=50, num_to_do=10, \
            motif=ile_motif, model=ile_model, composition=BaseFrequency('UCAG'))

        short_ile_mod_0 = Module('NCUACNN', '(((..((')
        short_ile_mod_1 = Module('NNUAUUGGGGN', '))......)))')
        short_ile_rule_0 = Rule(0, 0, 1, 10, 3)
        short_ile_rule_1 = Rule(0, 5, 1, 1, 2)
        short_ile_motif = Motif([short_ile_mod_0, short_ile_mod_1], \
            [short_ile_rule_0, short_ile_rule_1])

        short_helices = [PairedRegion('N'), PairedRegion('NN')]
        short_constants = [ConstantRegion('CUAC'), ConstantRegion('UAUUGGGG')]
        short_order = "H0 C0 H1 - H1 C1 H0"
        short_ile_model = SequenceModel(order=short_order, \
            constants=short_constants, \
            helices=short_helices, composition=BaseFrequency('UCAG'))

        self.short_ile_embedder = SequenceEmbedder(length=50, num_to_do=10, \
            motif=short_ile_motif, model=short_ile_model, \
            composition=BaseFrequency('UCAG'))


    def test_composition_change(self):
        """Changes in composition should propagate."""
        rr = str(self.ile_embedder.RandomRegion.Current)
        #for base in 'UCAG':
        #    assert base in rr
        #the above two lines should generally be true but fail stochastically
        self.ile_embedder.Composition = BaseFrequency('CG')
        self.assertEqual(self.ile_embedder.Model.Composition, \
            BaseFrequency('CG'))
        self.assertEqual(self.ile_embedder.RandomRegion.Composition, \
            BaseFrequency('CG'))
        self.ile_embedder.RandomRegion.refresh()
        self.assertEqual(len(self.ile_embedder.RandomRegion), 22)
        rr = str(self.ile_embedder.RandomRegion.Current)
        assert ('C' in rr or 'G' in rr)
        assert 'A' not in rr
        assert 'U' not in rr

    def test_choose_locations_too_short(self):
        """SequenceEmbedder _choose_locations should fail if too little space"""
        self.ile_embedder.Length = 28   #no positions left over
        self.assertRaises(ValueError, self.ile_embedder._choose_locations)
        self.ile_embedder.Length = 29   #one position left over
        self.assertRaises(ValueError, self.ile_embedder._choose_locations)

    def test_choose_locations_exact(self):
        """SequenceEmbedder _choose_locations should pick all locations"""
        self.ile_embedder.Length = 30   #two positions left: must both be filled
        for i in range(10):
            first, second = self.ile_embedder._choose_locations()
            self.assertEqual(first, 0)
            self.assertEqual(second, 1)
            
    def test_choose_locations_even(self):
        """SequenceEmbedder _choose_locations should pick locations evenly"""
        self.ile_embedder.Length = 31   #three positions left
        counts = {}
        for i in range(1000):
            key = tuple(self.ile_embedder._choose_locations())
            assert key[0] != key[1]
            curr = counts.get(key, 0)
            counts[key] = curr + 1
        expected = [333, 333, 333]
        observed = [counts[(0,1)], counts[(0,2)], counts[(1,2)]]
        self.assertSimilarFreqs(observed, expected)
        #make sure nothing else snuck in there
        self.assertEqual(counts[(0,1)]+counts[(0,2)]+counts[(1,2)], 1000)

    def test_choose_locations_with_replacement(self):
        """SequenceEmbedder _choose_locations can sample with replacement"""
        self.ile_embedder.Length = 28   #exact fit
        self.ile_embedder.WithReplacement = True
        for i in range(10):
            first, second = self.ile_embedder._choose_locations()
            self.assertEqual(first, 0)
            self.assertEqual(second, 0)

        self.ile_embedder.Length = 29 #one left over: can be 0,0 0,1 1,1
        counts = {}
        for i in range(1000):
            key = tuple(self.ile_embedder._choose_locations())
            curr = counts.get(key, 0)
            counts[key] = curr + 1
        expected = [250, 500, 250]
        observed = [counts[(0,0)], counts[(0,1)], counts[(1,1)]]
        self.assertSimilarFreqs(observed, expected)
        #make sure nothing else snuck in there
        self.assertEqual(counts[(0,0)]+counts[(0,1)]+counts[(1,1)], 1000)
  
    def test_insert_modules(self):
        """SequenceEmbedder _insert_modules should make correct sequence"""
        ile = self.ile_embedder
        ile.Length = 50
        ile.RandomRegion.Current[:] = ['A'] * 22
        modules = list(ile.Model)
        ile.Positions = [0, 0]  #try inserting at first position
        self.assertEqual(str(ile), modules[0] + modules[1] + 'A'*22)

        ile.Positions = [3, 20]
        self.assertEqual(str(ile), 'A'*3+modules[0]+'A'*17+modules[1]+'A'*2)

    def test_refresh(self):
        """SequenceEmbedder refresh should change module sequences"""
        modules_before = list(self.ile_embedder.Model)
        random_before = str(self.ile_embedder.RandomRegion.Current)
        self.ile_embedder.refresh()
        random_after = str(self.ile_embedder.RandomRegion.Current)
        self.assertNotEqual(random_before, random_after)
        modules_after = list(self.ile_embedder.Model)
        for before, after in zip(modules_before, modules_after):
            self.assertNotEqual(before, after)
        #check that it works twice
        self.ile_embedder.refresh()
        random_third = str(self.ile_embedder.RandomRegion.Current)
        modules_third = list(self.ile_embedder.Model)
        self.assertNotEqual(random_third, random_before)
        self.assertNotEqual(random_third, random_after)
        for first, second, third in \
            zip(modules_before, modules_after, modules_third):
            self.assertNotEqual(first, third)
            self.assertNotEqual(second, third)

    def test_countMatches(self):
        """Shouldn't find any Ile matches if all the pairs are GU"""
        if not RNAFOLD_PRESENT:
            return
        self.ile_embedder.NumToDo = 100
        self.ile_embedder.Composition = BaseFrequency('GGGGGGGGGU')
        self.ile_embedder.Length = 40
        good_count = self.ile_embedder.countMatches()   
        self.assertEqual(good_count, 0)

    def test_countMatches_pass(self):
        """Should find some matches against a random background"""
        if not RNAFOLD_PRESENT:
            return
        self.ile_embedder.NumToDo = 100
        self.ile_embedder.Composition = BaseFrequency('UCAG')
        self.ile_embedder.Length = 40
        good_count = self.ile_embedder.countMatches()
        self.assertNotEqual(good_count, 0)

    def test_refresh_specific_position(self):
        """Should always find the module in the same position if specified"""
        first_module = Module('AAAAA', '(((((')
        second_module = Module('UUUUU', ')))))')
        rule_1 = Rule(0, 0, 1, 4, 5)
        helix = Motif([first_module, second_module], [rule_1])
        model = SequenceModel(constants=[ConstantRegion('AAAAA'), \
            ConstantRegion('UUUUU')], order='C0 - C1', \
            composition=BaseFrequency('A'))
        embedder = SequenceEmbedder(length=30, num_to_do=100, \
            motif=helix, model=model, composition=BaseFrequency('CG'), \
            positions=[3, 6])

        last = ''
        for i in range(100):
            embedder.refresh()
            curr = str(embedder)
            self.assertEqual(curr[3:8], 'AAAAA')
            self.assertEqual(curr[11:16], 'UUUUU')
            self.assertEqual(curr.count('A'), 5)
            self.assertEqual(curr.count('U'), 5)
            self.assertNotEqual(last, curr)
            last = curr
            
    def test_refresh_primers(self):
        """Module should appear in correct location with primers"""
        first_module = Module('AAAAA', '(((((')
        second_module = Module('UUUUU', ')))))')
        rule_1 = Rule(0, 0, 1, 4, 5)
        helix = Motif([first_module, second_module], [rule_1])
        model = SequenceModel(constants=[ConstantRegion('AAAAA'), \
            ConstantRegion('UUUUU')], order='C0 - C1', \
            composition=BaseFrequency('A'))
        embedder = SequenceEmbedder(length=30, num_to_do=100, \
            motif=helix, model=model, composition=BaseFrequency('CG'), \
            positions=[3, 6], primer_5 = 'UUU', primer_3 = 'AAA')

        last = ''
        for i in range(100):
            embedder.refresh()
            curr = str(embedder)
            self.assertEqual(curr[0:3], 'UUU')
            self.assertEqual(curr[6:11], 'AAAAA')
            self.assertEqual(curr[14:19], 'UUUUU')
            self.assertEqual(curr.count('A'), 8)
            self.assertEqual(curr.count('U'), 8)
            self.assertEqual(curr[-3:], 'AAA')
            self.assertNotEqual(last, curr)
            last = curr
 
    def xxx_test_count_long(self):
        self.ile_embedder.NumToDo = 100000
        self.ile_embedder.Composition = BaseFrequency('UCAG')
        print
        print "Extended helices"
        for length in range(30, 150):
            self.ile_embedder.Length = length
            good_count = self.ile_embedder.countMatches()   
            print "Length: %s Matches: %s/100000" % (length, good_count)
        print

    def xxx_test_count_short(self):
        self.short_ile_embedder.NumToDo = 10000
        self.short_ile_embedder.Composition = BaseFrequency('UCAG')
        print
        print "Minimal motif"
        for length in range(20, 150):
            self.short_ile_embedder.Length = length
            good_count = self.short_ile_embedder.countMatches()   
            print "Length: %s Matches: %s/10000" % (length, good_count)
        print



if __name__ == '__main__':
    main()

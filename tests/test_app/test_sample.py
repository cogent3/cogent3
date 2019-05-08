from cogent3 import LoadSeqs, DNA
from unittest import TestCase, main

from cogent3.app import sample, composable

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TranslateTests(TestCase):
    def _codon_positions(self, array_align):
        """correctly return codon positions"""
        aln = LoadSeqs(data=[('a', 'ACGACGACG'),
                             ('b', 'GATGATGAT')],
                       array_align=array_align)
        one = sample.take_codon_positions(1)
        got = one(aln)
        self.assertEqual(got.todict(), {'a': 'AAA', 'b': 'GGG'})

        two = sample.take_codon_positions(2)
        got = two(aln)
        self.assertEqual(got.todict(), {'a': 'CCC', 'b': 'AAA'})
        three = sample.take_codon_positions(3)
        got = three(aln)
        self.assertEqual(got.todict(), {'a': 'GGG', 'b': 'TTT'})

        one_two = sample.take_codon_positions(1, 2)
        got = one_two(aln)
        self.assertEqual(got.todict(), {'a': 'ACACAC', 'b': 'GAGAGA'})
        one_three = sample.take_codon_positions(1, 3)
        got = one_three(aln)
        self.assertEqual(got.todict(), {'a': 'AGAGAG', 'b': 'GTGTGT'})
        two_three = sample.take_codon_positions(2, 3)
        got = two_three(aln)
        self.assertEqual(got.todict(), {'a': 'CGCGCG', 'b': 'ATATAT'})

    def test_take_codon_positions_array_align(self):
        """correctly return codon positions from ArrayAlignment"""
        self._codon_positions(array_align=True)

    def test_take_codon_positions_alignment(self):
        """correctly return codon positions from Alignment"""
        self._codon_positions(array_align=False)

    def test_filter_degen(self):
        """just_nucs correctly identifies data with only nucleotides"""
        aln = LoadSeqs(data=[('a', 'ACGA-GACG'),
                             ('b', 'GATGATGYT')])
        degen = sample.omit_degenerates(moltype='dna')
        got = degen(aln)
        self.assertEqual(got.todict(), {'a': 'ACGAGAG',
                                        'b': 'GATGTGT'})

    def test_codon_positions_4fold_degen(self):
        '''codon_positions correctly return fourfold degenerate bases'''
        #                           **4---**4---
        aln = LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                             ('b', 'GCTTTTGTCAAT')], moltype=DNA)
        expect = dict([('a', 'AT'), ('b', 'TC')])
        ffold = sample.take_codon_positions(fourfold_degenerate=True)
        got = ffold(aln)
        self.assertEqual(got.todict(), expect)
        # error if no moltype
        with self.assertRaises(AssertionError):
            _ = sample.take_codon_positions(moltype=None)

    def test_take_named(self):
        '''returns collections containing named seqs'''
        select = sample.take_named_seqs('a', 'b')
        alns = [LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT'),
                               ('c', 'GC--GCGTTTAT'),
                               ('d', 'GCAAGCNNTTAT')]),
                LoadSeqs(data=[('a', 'GGAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT'),
                               ('c', 'GC--GCGTTTAT'),
                               ('d', 'GCAAGCNNTTAT')])]
        got = [aln.todict() for aln in map(select, alns) if aln]
        expected = [dict((('a', 'GCAAGCGTTTAT'),
                          ('b', 'GCTTTTGTCAAT'))),
                    dict((('a', 'GGAAGCGTTTAT'),
                          ('b', 'GCTTTTGTCAAT')))]
        self.assertEqual(got, expected)
        # return False if a named seq absent
        aln = LoadSeqs(data=[('c', 'GC--GCGTTTAT'),
                             ('d', 'GCAAGCNNTTAT')])
        got = select(aln)
        self.assertFalse(got)
        self.assertTrue(type(got) == composable.ErrorResult)

        # using negate
        select = sample.take_named_seqs('c', negate=True)
        alns = [LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT')]),
                LoadSeqs(data=[('a', 'GGAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT'),
                               ('c', 'GC--GCGTTTAT')])]
        got = [aln.todict() for aln in map(select, alns) if aln]
        expect = [dict(d) for d in [[('a', 'GCAAGCGTTTAT'),
                                     ('b', 'GCTTTTGTCAAT')],
                                    [('a', 'GGAAGCGTTTAT'),
                                     ('b', 'GCTTTTGTCAAT')]]]
        self.assertEqual(got, expect)

    def test_minlength(self):
        '''correctly identifies data with minimal length'''
        aln = LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                             ('b', 'GCTTTTGTCAAT')])

        # if using subtract_degen, fails if incorect moltype
        ml = sample.min_length(9, subtract_degen=True)
        with self.assertRaises(ValueError):
            ml(aln)

        # but works if subtract_degen is False
        ml = sample.min_length(9, subtract_degen=False)
        aln = ml(aln)
        self.assertEqual(len(aln), 12)
        # or if moltype provided
        ml = sample.min_length(9, subtract_degen=True, moltype='dna')
        aln = ml(aln)
        self.assertEqual(len(aln), 12)

        alns = [LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT')], moltype=DNA),
                LoadSeqs(data=[('a', 'GGAAGCGT'),
                               ('b', 'GCTTT-GT')], moltype=DNA)]
        ml = sample.min_length(9)
        got = [aln.todict() for aln in map(ml, alns) if aln]
        expected = [dict((('a', 'GCAAGCGTTTAT'),
                          ('b', 'GCTTTTGTCAAT')))]
        self.assertEqual(got, expected)

        # returns ErrorResult if nothing satisifies
        got = ml(alns[1])
        self.assertTrue(type(got) == sample.ErrorResult)

        alns = [LoadSeqs(data=[('a', 'GGAAGCGT'),
                               ('b', 'GCTTNGT')],
                         aligned=False, moltype=DNA)]
        ml = sample.min_length(6)
        got = [aln.todict() for aln in map(ml, alns) if aln]
        expected = [dict((('a', 'GGAAGCGT'),
                          ('b', 'GCTTNGT')))]
        self.assertEqual(got, expected)

        ml = sample.min_length(7)
        got = [aln.todict() for aln in map(ml, alns) if aln]
        expected = []
        self.assertEqual(got, expected)

    def test_fixedlength(self):
        '''correctly returns data with specified length'''
        aln = LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                             ('b', 'GCTTTTGTCAAT')])

        fl = sample.fixed_length(4)
        got = fl(aln)
        self.assertEqual(len(got), 4)
        fl = sample.fixed_length(9, moltype='dna')
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(list(got.moltype), list(DNA))

        alns = [LoadSeqs(data=[('a', 'GCAAGCGTTTAT'),
                               ('b', 'GCTTTTGTCAAT')], moltype=DNA),
                LoadSeqs(data=[('a', 'GGAAGCGT'),
                               ('b', 'GCTTT-GT')], moltype=DNA)]
        fl = sample.fixed_length(9)
        got = [a for a in map(fl, alns) if a]
        self.assertEqual(len(got[0]), 9)
        expected = dict((('a', 'GCAAGCGTT'),
                         ('b', 'GCTTTTGTC')))
        self.assertEqual(got[0].todict(), expected)

        fl = sample.fixed_length(600)
        got = [a for a in map(fl, alns) if a]
        expected = []
        self.assertEqual(got, expected)
        # returns ErrorResult if nothing satisifies
        got = fl(alns[0])
        self.assertTrue(type(got) == sample.ErrorResult)

        fl = sample.fixed_length(9, random=True)
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(set(aln.names), set('ab'))

        # these will be just a subset as sampling one triplet
        fl = sample.fixed_length(3, random=True, motif_length=3)
        d = LoadSeqs(data=[('a', 'GCAAGCGTGTAT'),
                           ('b', 'GCTACTGTCAAT')])
        expect = d.todict()
        got = fl(d)
        self.assertEqual(len(got), 3)
        for name, seq in got.todict().items():
            self.assertIn(seq, expect[name])

        fl = sample.fixed_length(9, start=2)
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(got.todict(), aln[2:11].todict())

        fl = sample.fixed_length(4, start='random')
        expect = aln.todict()
        got = fl(aln)
        self.assertEqual(len(got), 4)
        for name, seq in got.todict().items():
            self.assertIn(seq, expect[name])

    def test_omit_bad_seqs(self):
        """correctly omit bad sequences from an alignment"""
        data = {'s1': '-ACC--TT',
                's2': '-ACC--TT',
                's3': '-ACC--TT',
                's4': 'AACCGGTT',
                's5': 'AACCGGTT',
                's6': '--------'}
        aln = LoadSeqs(data=data, moltype=DNA)
        # unless exclude_just_gap=False, which should just return self
        dropbad = sample.omit_bad_seqs(exclude_just_gap=False)
        got = dropbad(aln)
        self.assertEqual(got.todict(), data)

        # with disallowed_frac=0.6, we should drop s4&5 too
        dropbad = sample.omit_bad_seqs(disallowed_frac=0.5)
        got = dropbad(aln)
        self.assertEqual(got.todict(),
                         {'s1': '-ACC--TT',
                          's2': '-ACC--TT',
                          's3': '-ACC--TT'})

        # with disallowed_frac=0.9, we should drop s6
        dropbad = sample.omit_bad_seqs(disallowed_frac=0.9)
        got = dropbad(aln)
        self.assertEqual(got.todict(),
                         {'s1': '-ACC--TT',
                          's2': '-ACC--TT',
                          's3': '-ACC--TT',
                          's4': 'AACCGGTT',
                          's5': 'AACCGGTT'})

    def test_omit_duplicated(self):
        """correctly drop duplicated sequences"""
        # strict omit_duplicated
        data = {"a": "ACGT", "b": "ACG-",  # identical excepting -
                "c": "ACGN",  # non-strict matches above
                "d": "ACGG", "e": "ACGG", "k": "ACGG",  # strict identical
                "f": "RAAA", "g": "YAAA",  # non-strict identical
                "h": "GGGG"}  # unique!
        seqs = LoadSeqs(data=data, aligned=False, moltype=DNA)

        # mask_degen = True : [{'a', 'c', 'b'}, {'k', 'd', 'e'},
        # {'g', 'f'}] are dupe sets. Only 'h' unique
        drop = sample.omit_duplicated(
            mask_degen=True, choose=None, moltype='dna')
        got = drop(seqs)
        self.assertEqual(got.todict(), {"h": "GGGG"})
        # mask_degen = False : [{'a', 'b'}, {'k', 'd', 'e'}]
        # c, f, g, h
        drop = sample.omit_duplicated(mask_degen=False,
                                      choose=None, moltype='dna')
        got = drop(seqs)
        expect = {"a": "ACGT", "b": "ACG-", "c": "ACGN", "f": "RAAA",
                  "g": "YAAA", "h": "GGGG"}
        self.assertEqual(got.todict(), expect)

        # choose longest
        seqs = LoadSeqs(data=data, aligned=True, moltype=DNA)
        drop = sample.omit_duplicated(
            mask_degen=True, choose='longest', moltype='dna')
        got = drop(seqs)
        expect = {"a": "ACGT",
                  "k": "ACGG",
                  "g": "YAAA",
                  "h": "GGGG"}
        self.assertEqual(got.todict(), expect)

        # choose random
        drop = sample.omit_duplicated(mask_degen=True,
                                      choose='random',
                                      moltype='dna')
        got1 = drop(seqs)
        seqnames = set(got1.names)
        duplicates = [{'a', 'c', 'b'}, {'k', 'd', 'e'}, {'g', 'f'}]
        # should only be one of each group
        for dupes in duplicates:
            self.assertTrue(len(dupes & seqnames) == 1)

    def test_concat(self):
        """returns concatenated alignment"""
        alns = [LoadSeqs(data=d, moltype=DNA)
                for d in [{'seq1': 'AAA', 'seq2': 'AAA', 'seq3': 'AAA'},
                          {'seq1': 'TTT', 'seq2': 'TTT',
                           'seq3': 'TTT', 'seq4': 'TTT'},
                          {'seq1': 'CC', 'seq2': 'CC', 'seq3': 'CC'}]]
        ccat = sample.concat(intersect=True)
        got = ccat(alns)
        self.assertEqual(got.todict(), {'seq1': 'AAATTTCC',
                                        'seq2': 'AAATTTCC',
                                        'seq3': 'AAATTTCC'})

        ccat = sample.concat(intersect=False)
        got = ccat(alns)
        self.assertEqual(got.todict(), {'seq1': 'AAATTTCC',
                                        'seq2': 'AAATTTCC',
                                        'seq3': 'AAATTTCC',
                                        'seq4': '???TTT??'})


if __name__ == '__main__':
    main()

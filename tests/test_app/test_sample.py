from unittest import TestCase, main

from cogent3 import DNA, make_aligned_seqs, make_unaligned_seqs
from cogent3.app import composable, sample
from cogent3.app.composable import NotCompleted
from cogent3.core import alignment


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TranslateTests(TestCase):
    def _codon_positions(self, array_align):
        """correctly return codon positions"""
        aln = make_aligned_seqs(
            data=[("a", "ACGACGACG"), ("b", "GATGATGAT")], array_align=array_align
        )
        one = sample.take_codon_positions(1)
        got = one(aln)
        self.assertEqual(got.to_dict(), {"a": "AAA", "b": "GGG"})

        two = sample.take_codon_positions(2)
        got = two(aln)
        self.assertEqual(got.to_dict(), {"a": "CCC", "b": "AAA"})
        three = sample.take_codon_positions(3)
        got = three(aln)
        self.assertEqual(got.to_dict(), {"a": "GGG", "b": "TTT"})

        one_two = sample.take_codon_positions(1, 2)
        got = one_two(aln)
        self.assertEqual(got.to_dict(), {"a": "ACACAC", "b": "GAGAGA"})
        one_three = sample.take_codon_positions(1, 3)
        got = one_three(aln)
        self.assertEqual(got.to_dict(), {"a": "AGAGAG", "b": "GTGTGT"})
        two_three = sample.take_codon_positions(2, 3)
        got = two_three(aln)
        self.assertEqual(got.to_dict(), {"a": "CGCGCG", "b": "ATATAT"})

    def test_take_codon_positions_array_align(self):
        """correctly return codon positions from ArrayAlignment"""
        self._codon_positions(array_align=True)

    def test_take_codon_positions_alignment(self):
        """correctly return codon positions from Alignment"""
        self._codon_positions(array_align=False)

    def test_filter_degen(self):
        """just_nucs correctly identifies data with only nucleotides"""
        aln = make_aligned_seqs(data=[("a", "ACGA-GACG"), ("b", "GATGATGYT")])
        degen = sample.omit_degenerates(moltype="dna")
        got = degen(aln)
        self.assertEqual(got.to_dict(), {"a": "ACGAGAG", "b": "GATGTGT"})
        self.assertIsInstance(got, alignment.ArrayAlignment)

        # no ungapped columns
        aln = make_aligned_seqs(data=[("a", "-C-A-G-C-"), ("b", "G-T-A-G-T")])
        got = degen(aln)
        self.assertIsInstance(got, composable.NotCompleted)

        # we get back the alignment type we passed in
        aln = make_aligned_seqs(
            data=[("a", "ACGA-GACG"), ("b", "GATGATGYT")], array_align=False
        )
        got = degen(aln)
        self.assertIsInstance(got, alignment.Alignment)

    def test_omit_gapped(self):
        """omit_gap_pos correctly drops aligned columns"""
        # array alignment
        data = [("a", "ACGA-GA-CG"), ("b", "GATGATG-AT")]
        aln = make_aligned_seqs(data=data)
        nogaps = sample.omit_gap_pos(moltype="dna", allowed_frac=0)  # default
        got = nogaps(aln)
        self.assertIsInstance(got, alignment.ArrayAlignment)
        expect = dict(a="ACGAGACG", b="GATGTGAT")
        self.assertEqual(got.to_dict(), expect)
        # standard alignment
        aln = make_aligned_seqs(data=data, array_align=False)
        got = nogaps(aln)
        self.assertIsInstance(got, alignment.Alignment)
        self.assertEqual(got.to_dict(), expect)
        # non-exclusive gaps
        not_all_gaps = sample.omit_gap_pos(moltype="dna")  # default
        expect = dict(a="ACGA-GACG", b="GATGATGAT")
        aln = make_aligned_seqs(data=data)
        got = not_all_gaps(aln)
        self.assertEqual(got.to_dict(), expect)
        aln = make_aligned_seqs(data=data, array_align=False)
        got = not_all_gaps(aln)
        self.assertEqual(got.to_dict(), expect)
        # with motif length
        not_all_gaps = sample.omit_gap_pos(
            moltype="dna", allowed_frac=0, motif_length=2
        )
        aln = make_aligned_seqs(data=data)
        expect = dict(a="ACGACG", b="GATGAT")
        got = not_all_gaps(aln)
        self.assertEqual(got.to_dict(), expect)

        # no ungapped columns returns NotCompleted
        aln = make_aligned_seqs(data=[("a", "-C-A-G-C-"), ("b", "G-T-A-G-T")])
        got = nogaps(aln)
        self.assertIsInstance(got, composable.NotCompleted)

    def test_codon_positions_4fold_degen(self):
        """codon_positions correctly return fourfold degenerate bases"""
        #                           **4---**4---
        aln = make_aligned_seqs(
            data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")], moltype=DNA
        )
        expect = dict([("a", "AT"), ("b", "TC")])
        ffold = sample.take_codon_positions(fourfold_degenerate=True)
        got = ffold(aln)
        self.assertEqual(got.to_dict(), expect)
        # error if no moltype
        with self.assertRaises(AssertionError):
            _ = sample.take_codon_positions(moltype=None)

    def test_take_named(self):
        """returns collections containing named seqs"""
        select = sample.take_named_seqs("a", "b")
        alns = [
            make_aligned_seqs(
                data=[
                    ("a", "GCAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                    ("d", "GCAAGCNNTTAT"),
                ]
            ),
            make_aligned_seqs(
                data=[
                    ("a", "GGAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                    ("d", "GCAAGCNNTTAT"),
                ]
            ),
        ]
        got = [aln.to_dict() for aln in map(select, alns) if aln]
        expected = [
            dict((("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT"))),
            dict((("a", "GGAAGCGTTTAT"), ("b", "GCTTTTGTCAAT"))),
        ]
        self.assertEqual(got, expected)
        # return False if a named seq absent
        aln = make_aligned_seqs(data=[("c", "GC--GCGTTTAT"), ("d", "GCAAGCNNTTAT")])
        got = select(aln)
        self.assertFalse(got)
        self.assertTrue(type(got) == composable.NotCompleted)

        # using negate
        select = sample.take_named_seqs("c", negate=True)
        alns = [
            make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")]),
            make_aligned_seqs(
                data=[
                    ("a", "GGAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                ]
            ),
        ]
        got = [aln.to_dict() for aln in map(select, alns) if aln]
        expect = [
            dict(d)
            for d in [
                [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
                [("a", "GGAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
            ]
        ]
        self.assertEqual(got, expect)

    def test_minlength(self):
        """correctly identifies data with minimal length"""
        aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])

        # if using subtract_degen, fails if incorect moltype
        ml = sample.min_length(9, subtract_degen=True)
        got = ml(aln)
        self.assertIsInstance(got, NotCompleted)
        self.assertEqual(got.type, "ERROR")

        # but works if subtract_degen is False
        ml = sample.min_length(9, subtract_degen=False)
        aln = ml(aln)
        self.assertEqual(len(aln), 12)
        # or if moltype provided
        ml = sample.min_length(9, subtract_degen=True, moltype="dna")
        aln = ml(aln)
        self.assertEqual(len(aln), 12)

        alns = [
            make_aligned_seqs(
                data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")], moltype=DNA
            ),
            make_aligned_seqs(data=[("a", "GGAAGCGT"), ("b", "GCTTT-GT")], moltype=DNA),
        ]
        ml = sample.min_length(9)
        got = [aln.to_dict() for aln in map(ml, alns) if aln]
        expected = [dict((("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")))]
        self.assertEqual(got, expected)

        # returns NotCompletedResult if nothing satisifies
        got = ml(alns[1])
        self.assertTrue(type(got) == sample.NotCompleted)

        alns = [
            make_unaligned_seqs(data=[("a", "GGAAGCGT"), ("b", "GCTTNGT")], moltype=DNA)
        ]
        ml = sample.min_length(6)
        got = [aln.to_dict() for aln in map(ml, alns) if aln]
        expected = [dict((("a", "GGAAGCGT"), ("b", "GCTTNGT")))]
        self.assertEqual(got, expected)

        ml = sample.min_length(7)
        got = [aln.to_dict() for aln in map(ml, alns) if aln]
        expected = []
        self.assertEqual(got, expected)

    def test_fixedlength(self):
        """correctly returns data with specified length"""
        aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])

        fl = sample.fixed_length(4)
        got = fl(aln)
        self.assertEqual(len(got), 4)
        fl = sample.fixed_length(9, moltype="dna")
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(list(got.moltype), list(DNA))

        alns = [
            make_aligned_seqs(
                data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")], moltype=DNA
            ),
            make_aligned_seqs(data=[("a", "GGAAGCGT"), ("b", "GCTTT-GT")], moltype=DNA),
        ]
        fl = sample.fixed_length(9)
        got = [a for a in map(fl, alns) if a]
        self.assertEqual(len(got[0]), 9)
        expected = dict((("a", "GCAAGCGTT"), ("b", "GCTTTTGTC")))
        self.assertEqual(got[0].to_dict(), expected)

        fl = sample.fixed_length(600)
        got = [a for a in map(fl, alns) if a]
        expected = []
        self.assertEqual(got, expected)
        # returns NotCompletedResult if nothing satisifies
        got = fl(alns[0])
        self.assertTrue(type(got) == sample.NotCompleted)

        fl = sample.fixed_length(9, random=True)
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(set(aln.names), set("ab"))

        # these will be just a subset as sampling one triplet
        fl = sample.fixed_length(3, random=True, motif_length=3)
        d = make_aligned_seqs(data=[("a", "GCAAGCGTGTAT"), ("b", "GCTACTGTCAAT")])
        expect = d.to_dict()
        got = fl(d)
        self.assertEqual(len(got), 3)
        for name, seq in got.to_dict().items():
            self.assertIn(seq, expect[name])
        # as above, but with moltype defined
        fl = sample.fixed_length(3, random=True, motif_length=3, moltype="dna")
        got = fl(d)
        self.assertEqual(len(got), 3)
        for name, seq in got.to_dict().items():
            self.assertIn(seq, expect[name])

        fl = sample.fixed_length(9, start=2)
        got = fl(aln)
        self.assertEqual(len(got), 9)
        self.assertEqual(got.to_dict(), aln[2:11].to_dict())

        fl = sample.fixed_length(4, start="random")
        expect = aln.to_dict()
        got = fl(aln)
        self.assertEqual(len(got), 4)
        for name, seq in got.to_dict().items():
            self.assertIn(seq, expect[name])

    def test_omit_bad_seqs(self):
        """correctly omit bad sequences from an alignment"""
        data = {
            "s1": "---ACC---TT-",
            "s2": "---ACC---TT-",
            "s3": "---ACC---TT-",
            "s4": "--AACCG-GTT-",
            "s5": "--AACCGGGTTT",
            "s6": "AGAACCGGGTT-",
            "s7": "------------",
        }
        aln = make_aligned_seqs(data=data, moltype=DNA)
        # default just eliminates strict gap sequences
        dropbad = sample.omit_bad_seqs()
        got = dropbad(aln)
        expect = data.copy()
        del expect["s7"]
        self.assertEqual(got.to_dict(), expect)
        # providing a more stringent gap_frac
        dropbad = sample.omit_bad_seqs(gap_fraction=0.5)
        got = dropbad(aln)
        expect = data.copy()
        for n in ("s1", "s2", "s3", "s7"):
            del expect[n]
        self.assertEqual(got.to_dict(), expect)

        # setting quantile drops additional sequences
        dropbad = sample.omit_bad_seqs(quantile=6 / 7)
        got = dropbad(aln)
        expect = data.copy()
        for n in ("s6", "s7"):
            del expect[n]
        self.assertEqual(got.to_dict(), expect)

    def test_omit_duplicated(self):
        """correctly drop duplicated sequences"""
        # strict omit_duplicated
        data = {
            "a": "ACGT",
            "b": "ACG-",  # identical excepting -
            "c": "ACGN",  # non-strict matches above
            "d": "ACGG",
            "e": "ACGG",
            "k": "ACGG",  # strict identical
            "f": "RAAA",
            "g": "YAAA",  # non-strict identical
            "h": "GGGG",
        }  # unique!
        seqs = make_unaligned_seqs(data=data, moltype=DNA)

        # mask_degen = True : [{'a', 'c', 'b'}, {'k', 'd', 'e'},
        # {'g', 'f'}] are dupe sets. Only 'h' unique
        drop = sample.omit_duplicated(mask_degen=True, choose=None, moltype="dna")
        got = drop(seqs)
        self.assertEqual(got.to_dict(), {"h": "GGGG"})
        # mask_degen = False : [{'a', 'b'}, {'k', 'd', 'e'}]
        # c, f, g, h
        drop = sample.omit_duplicated(mask_degen=False, choose=None, moltype="dna")
        got = drop(seqs)
        expect = {
            "a": "ACGT",
            "b": "ACG-",
            "c": "ACGN",
            "f": "RAAA",
            "g": "YAAA",
            "h": "GGGG",
        }
        self.assertEqual(got.to_dict(), expect)

    def test_omit_duplicated_aligned(self):
        """omit_duplicated works on aligned sequences"""
        data = {
            "a": "ACGT",
            "b": "ACG-",  # identical excepting -
            "c": "ACGN",  # non-strict matches above
            "d": "ACGG",
            "e": "ACGG",
            "k": "ACGG",  # strict identical
            "f": "RAAA",
            "g": "YAAA",  # non-strict identical
            "h": "GGGG",
        }  # unique!
        # choose longest
        seqs = make_aligned_seqs(data=data, moltype=DNA)
        drop = sample.omit_duplicated(mask_degen=True, choose="longest", moltype="dna")
        got = drop(seqs)
        expect = {"a": "ACGT", "k": "ACGG", "g": "YAAA", "h": "GGGG"}
        self.assertEqual(got.to_dict(), expect)

        # choose random
        drop = sample.omit_duplicated(mask_degen=True, choose="random", moltype="dna")
        got1 = drop(seqs)
        seqnames = set(got1.names)
        duplicates = [{"a", "c", "b"}, {"k", "d", "e"}, {"g", "f"}]
        # should only be one of each group
        for dupes in duplicates:
            self.assertTrue(len(dupes & seqnames) == 1)

    def test_concat(self):
        """returns concatenated alignment"""
        alns = [
            make_aligned_seqs(data=d, moltype=DNA)
            for d in [
                {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
                {"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                {"seq1": "CC", "seq2": "CC", "seq3": "CC"},
            ]
        ]
        ccat = sample.concat(intersect=True)
        got = ccat(alns)
        self.assertEqual(
            got.to_dict(), {"seq1": "AAATTTCC", "seq2": "AAATTTCC", "seq3": "AAATTTCC"}
        )

        ccat = sample.concat(intersect=False)
        got = ccat(alns)
        self.assertEqual(
            got.to_dict(),
            {
                "seq1": "AAATTTCC",
                "seq2": "AAATTTCC",
                "seq3": "AAATTTCC",
                "seq4": "???TTT??",
            },
        )

    def test_concat_handles_moltype(self):
        """coerces to type"""
        alns = [
            make_aligned_seqs(data=d, moltype=DNA)
            for d in [
                {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
                {"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                {"seq1": "CC", "seq2": "CC", "seq3": "CC"},
            ]
        ]
        ccat = sample.concat()
        got = ccat(alns)
        self.assertIsInstance(got.moltype, type(DNA))

    def test_concat_validates_type(self):
        """raises TypeError if not known alignment type"""
        data = [
            {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
            make_aligned_seqs(
                data={"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                moltype=DNA,
            ),
        ]
        ccat = sample.concat()
        # triggered by first record
        with self.assertRaises(TypeError):
            ccat(data)

        # triggered by second record
        with self.assertRaises(TypeError):
            ccat(data[::-1])

        # triggered by no data
        with self.assertRaises(ValueError):
            ccat([])

    def test_trim_stop_codons(self):
        """trims stop codons using the specified genetic code"""
        trimmer = sample.trim_stop_codons()  # defaults to standard code
        seqs = make_unaligned_seqs(
            data={"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"}, moltype="dna"
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        self.assertEqual(got.to_dict(), expect)

        trimmer = sample.trim_stop_codons(gc=1)  # standard code
        seqs = make_unaligned_seqs(
            data={"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"}, moltype="dna"
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        self.assertEqual(got.to_dict(), expect)
        trimmer = sample.trim_stop_codons(gc=1)  # standard code
        aln = make_aligned_seqs(
            data={"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"}, moltype="dna"
        )
        got = trimmer(aln)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT---"}
        self.assertEqual(got.to_dict(), expect)

        # different genetic code
        trimmer = sample.trim_stop_codons(gc=2)  # mt code
        seqs = make_unaligned_seqs(
            data={"seq1": "AAATTTCCC", "seq2": "AAATTTAGA"}, moltype="dna"
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        self.assertEqual(got.to_dict(), expect)

    def test_take_n_seqs(self):
        """select specified number of sequences from a collection"""
        seqs1 = make_unaligned_seqs(
            data={
                "a": "ACGT",
                "b": "ACG-",
                "c": "ACGN",
                "d": "ACGG",
                "e": "ACGG",
                "k": "ACGG",
                "f": "RAAA",
                "g": "YAAA",
                "h": "GGGG",
            }
        )
        seqs2 = seqs1.take_seqs(["a", "c", "e", "g", "h"])

        # by order, fixed
        take = sample.take_n_seqs(3, fixed_choice=True)
        got = take(seqs1)
        self.assertEqual(len(got.names), 3)
        # this should return NotCompleted because it applies the names present in 1 to the next one
        got = take(seqs2)
        self.assertIsInstance(got, NotCompleted)

        take = sample.take_n_seqs(30)
        # this should fail because too few seqs
        got = take(seqs1)
        self.assertIsInstance(got, NotCompleted)

        # by order, not fixed
        take = sample.take_n_seqs(3, fixed_choice=False)
        got1 = take(seqs1)
        got2 = take(seqs2)
        self.assertNotEqual(set(got1.names), set(got2.names))

        # random choice, fixed
        take = sample.take_n_seqs(3, random=True, fixed_choice=True)
        self.assertEqual(take._fixed_choice, True)

        got1 = take(seqs2)
        got2 = take(seqs1)
        self.assertEqual(got1.names, got2.names)

        # random choice, not fixed
        take = sample.take_n_seqs(2, random=True, fixed_choice=False)
        self.assertEqual(take._fixed_choice, False)
        # testing this is hard, we simply expect the labels to differ on subsequent call
        # the probability of drawing a specific pair of names on one call is 1/(9 choose 2) = 1/36
        # at n = 11, the probability all the pairs will be identical is ~=0
        first_call = take(seqs1)
        for _ in range(11):
            got = take(seqs1)
            different = first_call.names != got.names
            if different:
                break

        self.assertTrue(different, msg="failed to generate different random sample")

        # try setting the seed
        take = sample.take_n_seqs(2, random=True, seed=123)
        got = take(seqs1)
        self.assertNotIsInstance(got, NotCompleted)


if __name__ == "__main__":
    main()

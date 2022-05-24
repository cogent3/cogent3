from unittest import TestCase, main

from cogent3 import (
    DNA,
    get_moltype,
    make_aligned_seqs,
    make_tree,
    make_unaligned_seqs,
)
from cogent3.align.align import make_generic_scoring_dict
from cogent3.app import align as align_app
from cogent3.app.align import (
    _combined_refseq_gaps,
    _gap_difference,
    _gap_union,
    _GapOffset,
    _gaps_for_injection,
    _merged_gaps,
    pairwise_to_multiple,
)
from cogent3.app.composable import NotCompleted
from cogent3.core.alignment import Aligned
from cogent3.core.location import gap_coords_to_map


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_seqs = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}

_nucleotide_models = [
    "JC69",
    "K80",
    "F81",
    "HKY85",
    "TN93",
    "GTR",
    "ssGN",
    "GN",
    "BH",
    "DT",
]

_codon_models = [
    "CNFGTR",
    "CNFHKY",
    "MG94HKY",
    "MG94GTR",
    "GY94",
    "H04G",
    "H04GK",
    "H04GGK",
    "GNC",
]


def make_pairwise(data, refseq_name, moltype="dna", array_align=False):
    """returns series of refseq, [(n, pwise aln),..]. All alignments are to ref_seq"""
    aln = make_aligned_seqs(
        data,
        array_align=array_align,
        moltype=moltype,
    )
    refseq = aln.get_seq(refseq_name)
    pwise = [
        (n, aln.take_seqs([refseq_name, n]).omit_gap_pos())
        for n in aln.names
        if n != refseq_name
    ]
    return refseq, pwise


def make_aligned(gaps_lengths, seq, name="seq1"):
    seq = seq.moltype.make_seq(seq, name=name)
    return Aligned(gap_coords_to_map(gaps_lengths, len(seq)), seq)


class RefalignmentTests(TestCase):
    seqs = make_unaligned_seqs(_seqs, moltype=DNA)

    def test_align_to_ref(self):
        """correctly aligns to a reference"""
        aligner = align_app.align_to_ref(ref_seq="Human")
        aln = aligner(self.seqs)
        expect = {
            "Bandicoot": "---NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
            "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
            "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "Rhesus": "GCCAGCTCATTACAGCATGAGAAC---AGTTTGTTACTCACT",
        }
        self.assertEqual(aln.to_dict(), expect)

    def test_align_to_ref_generic_moltype(self):
        """tests when the moltype is generic"""
        test_moltypes = ["text", "rna", "protein", "protein_with_stop", "bytes", "ab"]
        for test_moltype in test_moltypes:
            aligner = align_app.align_to_ref(moltype=test_moltype)
            self.assertEqual(aligner._moltype.label, test_moltype)
            self.assertEqual(
                aligner._kwargs["S"],
                make_generic_scoring_dict(10, get_moltype(test_moltype)),
            )

    def test_align_to_ref_result_has_moltype(self):
        """aligned object has correct moltype"""
        aligner = align_app.align_to_ref(moltype="dna")
        got = aligner(self.seqs)
        self.assertEqual(got.moltype.label, "dna")

    def test_merged_gaps(self):
        """correctly merges gaps"""
        a = dict([(2, 3), (4, 9)])
        b = dict([(2, 6), (8, 5)])
        # omitting one just returns the other
        self.assertIs(_merged_gaps(a, {}), a)
        self.assertIs(_merged_gaps({}, b), b)
        got = _merged_gaps(a, b)
        self.assertEqual(got, [(2, 6), (4, 9), (8, 5)])

    def test_aln_to_ref_known(self):
        """correctly recapitulates known case"""
        orig = make_aligned_seqs(
            {
                "Ref": "CAG---GAGAACAGAAACCCAT--TACTCACT",
                "Qu1": "CAG---GAGAACAG---CCCGTGTTACTCACT",
                "Qu2": "CAGCATGAGAACAGAAACCCGT--TA---ACT",
                "Qu3": "CAGCATGAGAACAGAAACCCGT----CTCACT",
                "Qu4": "CAGCATGAGAACAGAAACCCGTGTTACTCACT",
                "Qu5": "CAG---GAGAACAG---CCCAT--TACTCACT",
                "Qu6": "CAG---GA-AACAG---CCCAT--TACTCACT",
                "Qu7": "CAG---GA--ACAGA--CCCGT--TA---ACT",
            },
            moltype="dna",
        )
        expect = orig.to_dict()
        aligner = align_app.align_to_ref(ref_seq="Ref")
        aln = aligner(orig.degap())
        self.assertEqual(aln.to_dict(), expect)

    def test_gap_union(self):
        """correctly identifies the union of all gaps"""
        # fails if not all sequences same
        seq = DNA.make_seq("AACCCGTT")
        all_gaps = dict([(0, 3), (2, 1), (5, 3), (6, 3)])
        make_aligned(all_gaps, seq)
        gap_sets = [
            dict([(5, 1), (6, 3)]),
            dict([(2, 1), (5, 3)]),
            dict([(2, 1), (5, 1), (6, 2)]),
            dict([(0, 3)]),
        ]
        seqs = [make_aligned(gaps, seq) for gaps in gap_sets]
        got = _gap_union(seqs)
        self.assertEqual(got, dict(all_gaps))

        # must all be Aligned instances
        with self.assertRaises(TypeError):
            _gap_union(seqs + ["GGGGGGGG"])

        # must all have the same name
        with self.assertRaises(ValueError):
            _gap_union(seqs + [make_aligned({}, seq, name="blah")])

    def test_gap_difference(self):
        """correctly identifies the difference in gaps"""
        seq = DNA.make_seq("AACCCGTT")
        dict([(0, 3), (2, 1), (5, 3), (6, 3)])
        gap_sets = [
            dict([(5, 1), (6, 3)]),
            dict([(2, 1), (5, 3)]),
            dict([(2, 1), (5, 1), (6, 2)]),
            dict([(0, 3)]),
        ]
        seqs = [make_aligned(gaps, seq) for gaps in gap_sets]
        union = _gap_union(seqs)
        expects = [
            [dict([(0, 3), (2, 1)]), dict([(5, 2)])],
            [dict([(0, 3), (6, 3)]), {}],
            [dict([(0, 3)]), dict([(5, 2), (6, 1)])],
            [dict([(2, 1), (5, 3), (6, 3)]), {}],
        ]
        for seq, (plain, overlap) in zip(seqs, expects):
            seq_gaps = dict(seq.map.get_gap_coordinates())
            got_plain, got_overlap = _gap_difference(seq_gaps, union)
            self.assertEqual(got_plain, dict(plain))
            self.assertEqual(got_overlap, dict(overlap))

    def test_merged_gaps(self):
        """correctly handles gap values"""
        a_gaps = {0: 2}
        b_gaps = {2: 2}
        self.assertEqual(_merged_gaps(a_gaps, {}), a_gaps)
        self.assertEqual(_merged_gaps({}, b_gaps), b_gaps)

    def test_combined_refseq_gaps(self):
        union = dict([(0, 3), (2, 1), (5, 3), (6, 3)])
        gap_sets = [
            [(5, 1), (6, 3)],
            [(2, 1), (5, 3)],
            [(2, 1), (5, 1), (6, 2)],
            [(0, 3)],
        ]
        # for subset gaps, their alignment position is the
        # offset + their position + their gap length
        expects = [
            dict([(6, 2), (0, 3), (2, 1)]),
            dict([(0, 3), (10, 3)]),
            dict([(0, 3), (5 + 1 + 1, 2), (6 + 2 + 2, 1)]),
            dict([(2 + 3, 1), (5 + 3, 3), (6 + 3, 3)]),
        ]
        for i, gap_set in enumerate(gap_sets):
            got = _combined_refseq_gaps(dict(gap_set), union)
            self.assertEqual(got, expects[i])

        # if union gaps equals ref gaps
        got = _combined_refseq_gaps({2: 2}, {2: 2})
        self.assertEqual(got, {})

    def test_gaps_for_injection(self):
        # for gaps before any otherseq gaps, alignment coord is otherseq coord
        oseq_gaps = {2: 1, 6: 2}
        rseq_gaps = {0: 3}
        expect = {0: 3, 2: 1, 6: 2}
        seqlen = 50
        got = _gaps_for_injection(oseq_gaps, rseq_gaps, seqlen)
        self.assertEqual(got, expect)
        # for gaps after otherseq gaps seq coord is align coord minus gap
        # length totals
        got = _gaps_for_injection(oseq_gaps, {4: 3}, seqlen)
        expect = {2: 1, 3: 3, 6: 2}
        self.assertEqual(got, expect)
        got = _gaps_for_injection(oseq_gaps, {11: 3}, seqlen)
        expect = {2: 1, 6: 2, 8: 3}
        self.assertEqual(got, expect)
        # gaps beyond sequence length added to end of sequence
        got = _gaps_for_injection({2: 1, 6: 2}, {11: 3, 8: 3}, 7)
        expect = {2: 1, 6: 2, 7: 6}
        self.assertEqual(got, expect)

    def test_pairwise_to_multiple(self):
        """the standalone function constructs a multiple alignment"""
        expect = {
            "Ref": "CAG---GAGAACAGAAACCCAT--TACTCACT",
            "Qu1": "CAG---GAGAACAG---CCCGTGTTACTCACT",
            "Qu2": "CAGCATGAGAACAGAAACCCGT--TA---ACT",
            "Qu3": "CAGCATGAGAACAGAAACCCGT----CTCACT",
            "Qu7": "CAG---GA--ACAGA--CCCGT--TA---ACT",
            "Qu4": "CAGCATGAGAACAGAAACCCGTGTTACTCACT",
            "Qu5": "CAG---GAGAACAG---CCCAT--TACTCACT",
            "Qu6": "CAG---GA-AACAG---CCCAT--TACTCACT",
        }
        aln = make_aligned_seqs(expect, moltype="dna").omit_gap_pos()
        expect = aln.to_dict()
        for refseq_name in ["Qu3"]:
            refseq, pwise = make_pairwise(expect, refseq_name)
            got = pairwise_to_multiple(pwise, ref_seq=refseq, moltype=refseq.moltype)
            self.assertEqual(len(got), len(aln))
            orig = dict(pwise)
            _, pwise = make_pairwise(got.to_dict(), refseq_name)
            got = dict(pwise)
            # should be able to recover the original pairwise alignments
            for key, value in got.items():
                self.assertEqual(value.to_dict(), orig[key].to_dict(), msg=refseq_name)

            with self.assertRaises(TypeError):
                pairwise_to_multiple(pwise, "ACGG", DNA)

    def test_pairwise_to_multiple_2(self):
        """correctly handle alignments with gaps beyond end of query"""
        # cogent3.core.alignment.DataError: Not all sequences are the same length:
        # max is 425, min is 419
        def make_pwise(data, ref_name):
            result = []
            for n, seqs in data.items():
                result.append(
                    [n, make_aligned_seqs(data=seqs, moltype="dna", array_align=False)]
                )
            ref_seq = result[0][1].get_seq(ref_name)
            return result, ref_seq

        pwise = {
            "Platypus": {
                "Opossum": "-----------------GTGC------GAT-------------------------------CCAAAAACCTGTGTC--ACCGT--------GCC----CAGAGCCTCC----CTCAGGCCGCTCGGGGAG---TG-------GCCCCCCG--GC-GGAGGGCAGGGATGGGGAGT-AGGGGTGGCAGTC----GGAACTGGAAGAGCTT-TACAAACC---------GA--------------------GGCT-AGAGGGTC-TGCTTAC-------TTTTTACCTTGG------------GTTTG-CCAGGAGGTAG----------AGGATGA-----------------CTAC--ATCAAG----AGC------------TGGG-------------",
                "Platypus": "CAGGATGACTACATCAAGAGCTGGGAAGATAACCAGCAAGGAGATGAAGCTCTGGACACTACCAAAGACCCCTGCCAGAACGTGAAGTGCAGCCGACACAAGGTCTGCATCGCTCAGGGCTACCAGAGAGCCATGTGTATCAGCCGCAAGAAGCTGGAGCACAGGATCAAGCAGCCAGCCCTGAAACTCCATGGAAACAGAGAGAGCTTCTGCAAGCCTTGTCACATGACCCAGCTGGCCTCTGTCTGCGGCTCGGACGGACACACTTACAGCTCCGTGTGCAAACTGGAGCAGCAGGCCTGTCTGACCAGCAAGCAGCTGACAGTCAAGTGTGAAGGCCAGTGCCCGTGCCCCACCGATCATGTTCCAGCCTCCACCGCTGATGGAAAACAAGAGACCT",
            },
            "Wombat": {
                "Opossum": "GTGCGATCCAAAAACCTGTGTCACCGTGCCCAGAGCCTCCCTCAGGCCGCTCGG-GGAGTGGCCCCCCGGCGGAGGGCAGGGATGGGGAGTAGGGGTGGCAGTCGGAACTGGAAGAGCTTTACAAACCGAGGCTAGAGGGTCTGCTTACTTTTTACCTTGG------GTTT--GC-CAGGA---GGT----AGAGGATGACTACATCAAGAGCTGGG---------------------------",
                "Wombat": "--------CA----------TCACCGC-CCCTGCACC---------CGGCTCGGCGGAGGGGGATTCTAA-GGGGGTCAAGGATGGCGAG-ACCCCTGGCAATTTCA--TGGAGGA------CGAGCAATGGCT-----GTC-GTCCATCTCCCAGTATAGCGGCAAGATCAAGCACTGGAACCGCTTCCGAGACGATGACTACATCAAGAGCTGGGAGGACAGTCAGCAAGGAGATGAAGCGC",
            },
        }
        pwise, ref_seq = make_pwise(pwise, "Opossum")
        aln = pairwise_to_multiple(pwise, ref_seq, ref_seq.moltype)
        self.assertNotIsInstance(aln, NotCompleted)

        pwise = {
            "Platypus": {
                "Opossum": "-----------------GTGC------GAT-------------------------------CCAAAAACCTGTGTC",
                "Platypus": "CAGGATGACTACATCAAGAGCTGGGAAGATAACCAGCAAGGAGATGAAGCTCTGGACACTACCAAAGACCCCTGCC",
            },
            "Wombat": {
                "Opossum": "GTGCGATCCAAAAACCTGTGTC",
                "Wombat": "--------CA----------TC",
            },
        }
        pwise, ref_seq = make_pwise(pwise, "Opossum")
        aln = pairwise_to_multiple(pwise, ref_seq, ref_seq.moltype)
        self.assertNotIsInstance(aln, NotCompleted)


class ProgressiveAlignment(TestCase):
    seqs = make_unaligned_seqs(_seqs, moltype=DNA)
    treestring = "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus:0.06," "Human:0.0):0.04);"

    def test_progressive_align_protein_moltype(self):
        """tests guide_tree is None and moltype is protein"""
        from cogent3 import load_aligned_seqs

        seqs = load_aligned_seqs("data/nexus_aa.nxs", moltype="protein")
        seqs = seqs.degap()
        seqs = seqs.take_seqs(["Rat", "Cow", "Human", "Mouse", "Whale"])
        aligner = align_app.progressive_align(model="WG01")
        got = aligner(seqs)
        self.assertNotIsInstance(got, NotCompleted)
        aligner = align_app.progressive_align(model="protein")
        got = aligner(seqs)
        self.assertNotIsInstance(got, NotCompleted)

    def test_progressive_align_nuc(self):
        """progressive alignment with nuc models"""
        aligner = align_app.progressive_align(model="TN93", distance="TN93")
        aln = aligner(self.seqs)
        expect = {
            "Rhesus": "GCCAGCTCATTACAGCATGAGAACAG---TTTGTTACTCACT",
            "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAG---TTTATTGTCCAAC",
            "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
        }
        got = aln.to_dict()
        self.assertEqual(got, expect)

        # using default
        aligner = align_app.progressive_align(model="TN93", distance="TN93")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        self.assertEqual(aln.moltype, aligner._moltype)
        # todo the following is not robust across operating systems
        # so commenting out for now, but needs to be checked
        # expect = {'Human': 'GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
        #           'Rhesus': 'GCCAGCTCATTACAGCATGAGAA---CAGTTTGTTACTCACT',
        #           'Bandicoot': 'NACTCATTAATGCTTGAAACCAG---CAGTTTATTGTCCAAC',
        #           'FlyingFox': 'GCCAGCTCTTTACAGCATGAGAA---CAGTTTATTATACACT'}
        # got = aln.to_dict()
        # self.assertEqual(got, expect)

    def test_progressive_fails(self):
        """should return NotCompletedResult along with message"""
        # Bandicoot has an inf-frame stop codon
        seqs = make_unaligned_seqs(
            data={"Human": "GCCTCA", "Rhesus": "GCCAGCTCA", "Bandicoot": "TGATCATTA"},
            moltype="dna",
        )
        aligner = align_app.progressive_align(model="codon")
        got = aligner(seqs)
        self.assertTrue(type(got), NotCompleted)

    def test_progress_with_guide_tree(self):
        """progressive align works with provided guide tree"""
        tree = make_tree(treestring=self.treestring)
        aligner = align_app.progressive_align(
            model="nucleotide", guide_tree=self.treestring
        )
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="nucleotide", guide_tree=tree)
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        # even if it has underscores in name
        treestring = (
            "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus_macaque:0.06," "Human:0.0):0.04);"
        )
        aligner = align_app.progressive_align(model="nucleotide", guide_tree=treestring)
        data = self.seqs.to_dict()
        data["Rhesus macaque"] = data.pop("Rhesus")
        seqs = make_unaligned_seqs(data)
        aln = aligner(seqs)
        self.assertEqual(len(aln), 42)
        # guide tree with no lengths raises value error
        with self.assertRaises(ValueError):
            _ = align_app.progressive_align(
                model="nucleotide",
                guide_tree="(Bandicoot,FlyingFox,(Rhesus_macaque,Human));",
            )

    def test_progressive_align_codon(self):
        """progressive alignment with codon models"""
        aligner = align_app.progressive_align(model="GY94")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="codon")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)

    def test_pickle_progressive_align(self):
        """test progressive_align is picklable"""
        from pickle import dumps, loads

        aligner = align_app.progressive_align(model="codon")
        aln = aligner(self.seqs)
        got = loads(dumps(aln))
        self.assertTrue(got)

    def test_with_genetic_code(self):
        """handles genetic code argument"""
        aligner = align_app.progressive_align(model="GY94", gc="2")
        # the 'TGA' codon is a sense codon in vertebrate mitochondrial
        self.assertTrue("TGA" in aligner._model.get_motifs())
        aligner = align_app.progressive_align(model="codon")
        # but a stop codon in the standard nuclear
        self.assertTrue("TGA" not in aligner._model.get_motifs())
        # try using a nuclear
        with self.assertRaises(TypeError):
            aligner = align_app.progressive_align(model="nucleotide", gc="2")

    def test_progressive_align_protein(self):
        """progressive alignment with protein models"""
        seqs = self.seqs.get_translation()
        aligner = align_app.progressive_align(model="WG01", guide_tree=self.treestring)
        aln = aligner(seqs)
        self.assertEqual(len(aln), 14)
        aligner = align_app.progressive_align(
            model="protein", guide_tree=self.treestring
        )
        aln = aligner(seqs)
        self.assertEqual(len(aln), 14)


class GapOffsetTests(TestCase):
    def test_empty(self):
        """create an empty offset"""
        goff = _GapOffset({})
        for i in range(4):
            self.assertEqual(goff[i], 0)

        goff = _GapOffset({}, invert=True)
        for i in range(4):
            self.assertEqual(goff[i], 0)

    def test_repr_str(self):
        """repr and str work"""
        goff = _GapOffset({}, invert=True)
        for func in (str, repr):
            self.assertEqual(func(goff), "{}")

    def test_gap_offset(self):
        goff = _GapOffset({1: 2, 3: 4})
        self.assertEqual(goff.min_pos, 1)
        self.assertEqual(goff.max_pos, 3)
        self.assertEqual(goff.total, 6)
        self.assertEqual(goff[0], 0)
        self.assertEqual(goff[1], 0)
        self.assertEqual(goff[2], 2)
        self.assertEqual(goff[3], 2)
        self.assertEqual(goff[4], 6)

    def test_gap_offset_invert(self):
        aln2seq = _GapOffset({2: 1, 5: 2, 7: 2}, invert=True)
        self.assertEqual(aln2seq._store, {3: 1, 2: 0, 8: 3, 6: 1, 12: 5, 10: 3})
        self.assertEqual(aln2seq.max_pos, 12)
        self.assertEqual(aln2seq.min_pos, 2)
        self.assertEqual(aln2seq[11], 3)
        seq2aln = _GapOffset({2: 1, 5: 2, 7: 2})
        for seq_pos in range(20):
            aln_pos = seq_pos + seq2aln[seq_pos]
            self.assertEqual(aln_pos - aln2seq[aln_pos], seq_pos)


if __name__ == "__main__":
    main()

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
    _gap_difference,
    _gap_union,
    _GapOffset,
    _map_ref_gaps_to_seq,
    _merged_gaps,
)
from cogent3.app.composable import NotCompleted
from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.location import (
    LostSpan,
    Map,
    Span,
    _gap_insertion_data,
    _gap_pos_to_map,
    _interconvert_seq_aln_coords,
    gap_coords_to_map,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2021.5.7a"
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

    def test__gap_insertion_data(self):
        """identifies gap locations and lengths"""
        seq = DNA.make_seq("AACCCCCCGGG")
        # no gaps
        m = _gap_pos_to_map([], [], len(seq))
        got = _gap_insertion_data(Aligned(m, seq))
        expect = [], []
        self.assertEqual(got, expect)
        # one gap
        gap_positions = [2]
        gap_lengths = [3]
        m = _gap_pos_to_map(gap_positions, gap_lengths, len(seq))
        got = _gap_insertion_data(Aligned(m, seq))
        expect = list(zip(gap_positions, gap_lengths)), [0]
        self.assertEqual(got, expect)
        # two gaps
        gap_positions = [2, 8]
        gap_lengths = [3, 1]
        m = _gap_pos_to_map(gap_positions, gap_lengths, len(seq))
        got = _gap_insertion_data(Aligned(m, seq))
        expect = list(zip(gap_positions, gap_lengths)), [0, 3]
        self.assertEqual(got, expect)

    def test__merged_gaps(self):
        """correctly merges gaps"""
        with self.assertRaises(ValueError):
            _merged_gaps([], [], function="blah")

        a = [(2, 3), (4, 9)]
        b = [(2, 6), (8, 5)]
        # omitting one just returns the other
        self.assertIs(_merged_gaps(a, []), a)
        self.assertIs(_merged_gaps([], b), b)
        # specifying 'sum' means (2, 9)
        got = _merged_gaps(a, b, function="sum")
        self.assertEqual(got, [(2, 9), (4, 9), (8, 5)])
        # specifying 'max' means (2, 6)
        got = _merged_gaps(a, b, function="max")
        self.assertEqual(got, [(2, 6), (4, 9), (8, 5)])

    def test__gap_pos_to_map(self):
        """correctly converts a gap positions to a Map"""
        gap_positions = [2, 9]
        gap_lengths = [3, 1]
        map = _gap_pos_to_map(gap_positions, gap_lengths, 20)
        self.assertTrue(all(map.spans[i].lost for i in (1, 3)))
        self.assertEqual(len(map), 24)
        # no gaps
        map = _gap_pos_to_map([], [], 20)
        self.assertEqual(len(map), 20)
        self.assertEqual(len(map.spans), 1)
        # gap at start
        gap_positions = [0]
        gap_lengths = [3]
        map = _gap_pos_to_map(gap_positions, gap_lengths, 20)
        self.assertTrue(map.spans[0].lost and not map.spans[1].lost)
        self.assertEqual(len(map), 23)
        self.assertEqual(len(map.spans), 2)
        # gap at end
        gap_positions = [20]
        gap_lengths = [3]
        map = _gap_pos_to_map(gap_positions, gap_lengths, 20)
        self.assertTrue(map.spans[-1].lost and not map.spans[0].lost)
        self.assertEqual(len(map), 23)
        self.assertEqual(len(map.spans), 2)
        # fail if pos beyond sequence
        with self.assertRaises(ValueError):
            _gap_pos_to_map([40], [2], 20)

    def test__map_ref_gaps_to_seq(self):
        """correctly handle case where ref and curr ref are equal"""
        aln = make_aligned_seqs(
            {
                "Ref": "CAG---GAGAACAGAAACCCAT--TACTCACT",
                "Qu2": "CAGCATGAGAACAGAAACCCGT--TA---ACT",
            },
            array_align=False,
            moltype="dna",
        )
        ref = aln.named_seqs["Ref"]
        other = aln.named_seqs["Qu2"]
        got = _map_ref_gaps_to_seq(ref, ref, other)
        self.assertIs(got, other.map)

    def test_seq_2_aln_coords(self):
        """correctly converts sequence coordinates to alignment coordinates"""
        got = _interconvert_seq_aln_coords([], [], 20, seq_pos=True)
        assert got == 20

        # aligned seq has one gap after
        got = _interconvert_seq_aln_coords([(22, 3)], [0], 20, seq_pos=True)
        assert got == 20

        # aligned seq has one 2bp gap before
        got = _interconvert_seq_aln_coords([(12, 2)], [0], 20, seq_pos=True)
        assert got == 22

        # aligned seq has two 2bp gap before
        got = _interconvert_seq_aln_coords([(10, 2), (12, 2)], [0], 20, seq_pos=True)
        assert got == 24

    def test_aln_2_seq_coord(self):
        """correctly converts alignment coordinates to sequence coordinates"""
        got = _interconvert_seq_aln_coords([], [], 20, seq_pos=False)
        assert got == 20

        # query has one gap after
        got = _interconvert_seq_aln_coords([(22, 3)], [0], 20, seq_pos=False)
        assert got == 20

        # query has one 2bp gap before
        got = _interconvert_seq_aln_coords([(12, 2)], [0], 20, seq_pos=False)
        assert got == 18

        # query has two 2bp gap before
        got = _interconvert_seq_aln_coords([(10, 2), (12, 2)], [0], 20, seq_pos=False)
        assert got == 16

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
        final_seq = make_aligned(all_gaps, seq)
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
        all_gaps = dict([(0, 3), (2, 1), (5, 3), (6, 3)])
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

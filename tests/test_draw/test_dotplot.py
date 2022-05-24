from unittest import TestCase, main

from cogent3 import DNA, make_unaligned_seqs
from cogent3.core.alignment import Aligned, ArrayAlignment
from cogent3.draw.dotplot import (
    Dotplot,
    _convert_coords_for_scatter,
    _convert_input,
    get_align_coords,
    len_seq,
    not_gap,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TestUtilFunctions(TestCase):
    def test_converting_coords(self):
        """convert [(x1,y1), (x2,y2),..] to plotly style"""
        got = _convert_coords_for_scatter([[(0, 1), (2, 3)], [(3, 4), (2, 8)]])
        expect = [0, 2, None, 3, 2], [1, 3, None, 4, 8]
        self.assertEqual(got, expect)

    def test_len_seq(self):
        """returns length of sequence minus gaps"""
        m, seq = DNA.make_seq("ACGGT--A").parse_out_gaps()
        self.assertEqual(len_seq(m), 6)

    def test_not_gap(self):
        """distinguishes Map instances that include gap or not"""
        m, seq = DNA.make_seq("ACGGT--A").parse_out_gaps()
        self.assertTrue(not_gap(m[0]))
        self.assertFalse(not_gap(m[5]))

    def test_convert_input(self):
        """converts data for dotplotting"""
        m, seq = DNA.make_seq("ACGGT--A").parse_out_gaps()
        aligned_seq = Aligned(m, seq)
        mapped_gap, new_seq = _convert_input(aligned_seq, None)
        self.assertIs(new_seq.moltype, DNA)
        self.assertIs(mapped_gap, m)
        self.assertIs(new_seq, seq)
        mapped_gap, new_seq = _convert_input("ACGGT--A", DNA)
        self.assertEqual(str(mapped_gap), str(m))
        self.assertEqual(str(new_seq), str(seq))

    def test_get_align_coords(self):
        """correctly returns the alignment coordinates"""
        # 01234  5
        # ACGGT--A
        #   012345
        # --GGTTTA
        m1, seq1 = DNA.make_seq("ACGGT--A").parse_out_gaps()
        m2, seq2 = DNA.make_seq("--GGTTTA").parse_out_gaps()
        x, y = get_align_coords(m1, m2)
        expect = [2, 4, None, 5, 5], [0, 2, None, 5, 5]
        self.assertEqual((x, y), expect)

        # we have no gaps, so coords will be None
        m1, s1 = seq1.parse_out_gaps()
        m2, s2 = seq2.parse_out_gaps()
        self.assertEqual(get_align_coords(m1, m2), None)

        # unless we indicate the seqs came from an Alignment
        m1, seq1 = DNA.make_seq("ACGGTTTA").parse_out_gaps()
        m2, seq2 = DNA.make_seq("GGGGTTTA").parse_out_gaps()
        x, y = get_align_coords(m1, m2, aligned=True)
        self.assertEqual((x, y), ([0, len(seq1)], [0, len(seq1)]))

        # raises an exception if the Aligned seqs are different lengths
        m1, seq1 = DNA.make_seq("ACGGTTTA").parse_out_gaps()
        m2, seq2 = DNA.make_seq("GGGGTT").parse_out_gaps()
        with self.assertRaises(AssertionError):
            get_align_coords(m1, m2, aligned=True)

    def test_display2d(self):
        """correctly constructs a Display2d"""
        dp = Dotplot("-TGATGTAAGGTAGTT", "CTGG---AAG---GGT", is_aligned=True, window=5)
        expect = [0, 2, None, 6, 8, None, 12, 14], [1, 3, None, 4, 6, None, 7, 9]
        self.assertEqual(dp._aligned_coords, expect)
        dp._build_fig()
        traces = dp.traces
        self.assertEqual(len(traces), 2)  # no rev complement
        # we nudge alignment coordinates by 0.2 on x-axis
        expect = [0.2, 2.2, None, 6.2, 8.2, None, 12.2, 14.2]
        self.assertEqual(traces[-1].x, expect)
        self.assertEqual(traces[-1].name, "Alignment")
        self.assertEqual(traces[0].name, "+ strand")
        # check against hand calculated coords
        expect_x = [6, 14, None, 2, 12, None, 3, 11, None, 0, 6]
        expect_y = [0, 8, None, 0, 10, None, 2, 10, None, 2, 8]
        self.assertEqual(traces[0].x, expect_x)
        self.assertEqual(traces[0].y, expect_y)

    def test_align_without_gaps(self):
        """dotplot has alignment coordinates if no gaps"""
        aln = ArrayAlignment(
            {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
        )
        aln_plot = aln.dotplot("seq1")
        self.assertNotEqual(aln_plot._aligned_coords, None)

    def test_dotplot_seqcoll(self):
        """dotplot sequence collection, gaps are removed"""
        seqs = make_unaligned_seqs(
            {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
        )
        dp = seqs.dotplot("seq1", "seq3")
        self.assertNotEqual(dp._aligned_coords, None)
        self.assertEqual(len(dp.seq1), 4)
        self.assertEqual(len(dp.seq2), 3)

    def test_dotplot_single(self):
        """dotplot with single sequence should not fail"""
        seqs = make_unaligned_seqs({"seq1": "ACGG"}, moltype="dna")
        dp = seqs.dotplot()
        self.assertEqual(dp.seq1, dp.seq2)

    def test_dotplot_missing(self):
        """fail if a sequence name not present"""
        seqs = make_unaligned_seqs(
            {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
        )
        with self.assertRaises(ValueError):
            _ = seqs.dotplot("seq1", "seq5")
        with self.assertRaises(ValueError):
            _ = seqs.dotplot("seq5", "seq1")
        with self.assertRaises(ValueError):
            _ = seqs.dotplot("seq5", "seq6")

    def test_dotplot_title(self):
        """setting empty string title works"""
        seqs = make_unaligned_seqs(
            {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
        )
        dp = seqs.dotplot("seq1", "seq3", title="")
        self.assertEqual(dp.figure.layout.title, "")


if __name__ == "__main__":
    main()

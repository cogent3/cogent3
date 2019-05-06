from unittest import TestCase, main

from cogent3 import DNA
from cogent3.core.alignment import Aligned
from cogent3.draw.dotplot_2 import (_convert_coords_for_scatter, len_seq,
                                    not_gap, _convert_input, get_align_coords,
                                    Display2D, )

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
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
        m, seq = DNA.make_seq('ACGGT--A').parse_out_gaps()
        self.assertEqual(len_seq(m), 6)

    def test_not_gap(self):
        """distinguishes Map instances that include gap or not"""
        m, seq = DNA.make_seq('ACGGT--A').parse_out_gaps()
        self.assertTrue(not_gap(m[0]))
        self.assertFalse(not_gap(m[5]))

    def test_convert_input(self):
        """converts data for dotplotting"""
        m, seq = DNA.make_seq('ACGGT--A').parse_out_gaps()
        aligned_seq = Aligned(m, seq)
        mapped_gap, new_seq = _convert_input(aligned_seq, None)
        self.assertIs(new_seq.moltype, DNA)
        self.assertIs(mapped_gap, m)
        self.assertIs(new_seq, seq)
        mapped_gap, new_seq = _convert_input('ACGGT--A', DNA)
        self.assertEqual(str(mapped_gap), str(m))
        self.assertEqual(str(new_seq), str(seq))

    def test_get_align_coords(self):
        """correctly returns the alignment coordinates"""
        # 01234  5
        # ACGGT--A
        #   012345
        # --GGTTTA
        m1, seq1 = DNA.make_seq('ACGGT--A').parse_out_gaps()
        m2, seq2 = DNA.make_seq('--GGTTTA').parse_out_gaps()
        x, y = get_align_coords(m1, m2)
        expect = [2, 4, None, 5, 5], [0, 2, None, 5, 5]
        self.assertEqual((x, y), expect)

        # we have no gaps, so coords will be None
        m1, s1 = seq1.parse_out_gaps()
        m2, s2 = seq2.parse_out_gaps()
        self.assertEqual(get_align_coords(m1, m2), None)

    def test_display2d(self):
        """correctly constructs a Display2d"""
        dp = Display2D('-TGATGTAAGGTAGTT', 'CTGG---AAG---GGT')
        expect = [0, 2, None, 6, 8, None, 12, 14], [1, 3, None, 4, 6, None, 7,
                                                    9]
        self.assertEqual(dp._aligned_coords, expect)
        traces = dp.get_trace(window=5)
        self.assertEqual(len(traces), 2)  # no rev complement
        # we nudge alignment coordinates by 0.2 on x-axis
        expect = (0.2, 2.2, None, 6.2, 8.2, None, 12.2, 14.2)
        self.assertEqual(traces[-1].x, expect)
        self.assertEqual(traces[-1].name, 'Alignment')
        self.assertEqual(traces[0].name, '+ strand')
        # check against hand calculated coords
        expect_x = (6, 14, None, 2, 12, None, 3, 11, None, 0, 6)
        expect_y = (0, 8, None, 0, 10, None, 2, 10, None, 2, 8)
        self.assertEqual(traces[0].x, expect_x)
        self.assertEqual(traces[0].y, expect_y)


if __name__ == '__main__':
    main()

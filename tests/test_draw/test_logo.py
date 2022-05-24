from unittest import TestCase, main

from cogent3.draw.letter import get_character
from cogent3.draw.logo import _char_hts_as_lists, get_logo
from cogent3.util.dict_array import DictArrayTemplate


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class LogoTests(TestCase):
    """testing utility functions"""

    def test_get_logo(self):
        """returns Drawable"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        data = DictArrayTemplate(5, "ACGT").wrap(data)
        get_logo(data)

    def test_get_logo_missing(self):
        """copes with positions with no values"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.05, 0.8, 0.05, 0.1],
            [0, 0, 0, 0],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        data = DictArrayTemplate(5, "ACGT").wrap(data)
        get_logo(data)

    def test_get_logo_alt_input_type(self):
        """copes with positions with no values"""
        data = [
            {"A": 0.1, "C": 0.3, "G": 0.5, "T": 0.1},
            {"A": 0.05, "C": 0.8, "G": 0.05, "T": 0.1},
            {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0},
            {"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1},
            {"A": 0.6, "C": 0.15, "G": 0.05, "T": 0.2},
        ]
        get_logo(data)

        data[-2] = {}
        get_logo(data)

    def test_letter_methods(self):
        """exercising some Letter methods"""
        # shift
        l = get_character("G")
        self.assertEqual(l.x, 0)
        self.assertEqual(l.y, 0)
        l.shift(2, 2)
        self.assertEqual(l.x, 2)
        self.assertEqual(l.y, 2)
        # scale adjusts the scale attributes
        orig_width = l.scale_x
        orig_height = l.scale_y
        l.scale(x=0.5, y=2)
        self.assertEqual(l.scale_x, orig_width / 2)
        self.assertEqual(l.scale_y, orig_height * 2)
        # invert changes the degree attr
        l.rotate(180)
        self.assertEqual(l.degrees, 180)

    def test_input_conversion(self):
        """correctly convert a series of dicts or a DictArray to lists"""
        data = [dict(A=0.1, C=0.2), dict(A=0.1, C=0.2)]
        base = [("A", 0.1), ("C", 0.2)]
        expect = [base, base]
        got = _char_hts_as_lists(data)
        self.assertEqual(got, expect)
        #
        data = [dict(A=0.1, C=0.2), {}]
        base = [("A", 0.1), ("C", 0.2)]
        expect = [base, None]
        got = _char_hts_as_lists(data)
        self.assertEqual(got, expect)
        data = [dict(A=0.1, C=0.2), None]
        base = [("A", 0.1), ("C", 0.2)]
        expect = [base, None]
        got = _char_hts_as_lists(data)
        self.assertEqual(got, expect)


if __name__ == "__main__":
    main()

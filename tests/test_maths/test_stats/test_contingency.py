from unittest import TestCase, main

import numpy

from numpy.testing import assert_allclose

from cogent3.maths.stats.contingency import CategoryCounts, calc_expected
from cogent3.util.dict_array import DictArrayTemplate


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class ContingencyTests(TestCase):
    def test_chisq(self):
        """correctly compute chisq test"""
        table = CategoryCounts([[762, 327, 468], [484, 239, 477]])
        got = table.chisq_test()
        self.assertEqual(round(got.chisq, 5), 30.07015)
        self.assertEqual(got.df, 2)
        assert_allclose(got.pvalue, 2.95358918321e-07)

    def test_residuals(self):
        """correctly calculate residuals"""
        table = CategoryCounts([[762, 327], [484, 239]])
        assert_allclose(
            table.residuals.array,
            [[0.48099031, -0.71365306], [-0.59031133, 0.87585441]],
        )

    def test_chisq2(self):
        """constructed from 2D dict"""
        data = {
            "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
            "b": {"env1": 1, "env3": 1, "env2": 3},
        }
        table = CategoryCounts(data)
        got = table.chisq_test()
        assert_allclose(got.chisq, 3.02222222)
        data = {
            "AIDS": {"Males": 4, "Females": 2, "Both": 3},
            "No_AIDS": {"Males": 3, "Females": 16, "Both": 2},
        }
        table = CategoryCounts(data)
        got = table.chisq_test()
        assert_allclose(got.chisq, 7.6568405139833722)
        assert_allclose(got.pvalue, 0.0217439383468)

    def test_1D_counts(self):
        """correctly operate on a 1D count array"""
        table = CategoryCounts([762, 327])
        got = table.chisq_test()
        assert_allclose(got.chisq, 173.7603305785124)
        self.assertLess(got.pvalue, 2.2e-16)  # value from R
        _ = got._repr_html_()  # shouldn't fail
        self.assertIn("1.12e-39", str(got))  # used sci formatting

    def test_G_ind(self):
        """correctly produce G test of independence"""
        table = CategoryCounts([[762, 327, 468], [484, 239, 477]])
        got = table.G_independence(williams=True)
        self.assertEqual(got.df, 2)

    def test_G_ind_with_pseudocount(self):
        """G test of independence with pseudocount"""
        table = CategoryCounts([[762, 327, 0], [484, 239, 0]])
        got = table.G_independence(williams=True, pseudo_count=1)
        assert_allclose(table.observed.array + 1, got.observed.array)
        assert_allclose(got.expected.array, calc_expected(got.observed.array))

    def test_G_fit_with_expecteds(self):
        """compute G-fit with provided expecteds"""
        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)), expected=dict(zip(keys, exp)))

        got = table.G_fit()
        assert_allclose(got.G, 9.849234)
        assert_allclose(got.pvalue, 0.04304536)
        _ = got._repr_html_()  # shouldn't fail
        self.assertIn("0.0430", str(got))  # used normal formatting

    def test_assign_expected(self):
        """assign expected property"""
        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)))
        table.expected = dict(zip(keys, exp))
        got = table.G_fit()
        assert_allclose(got.G, 9.849234)
        table.expected = None
        _ = table.G_fit()

    def test_zero_observeds(self):
        """raises ValueError"""
        with self.assertRaises(ValueError):
            CategoryCounts(dict(a=0, b=0))

    def test_shuffling(self):
        """resampling works for G-independence"""
        table = CategoryCounts([[762, 327], [750, 340]])
        got = table.G_independence(shuffled=50)
        self.assertTrue(0 < got.pvalue < 1)  # a large interval
        got = table.chisq_test(shuffled=50)
        self.assertTrue(0 < got.pvalue < 1)  # a large interval

    def test_to_dict(self):
        """returns a dict of contents"""
        table = CategoryCounts([[762, 327], [750, 340]])
        got = table.to_dict()
        assert_allclose(got["residuals"][0][0], 0.23088925877536437)
        assert_allclose(got["observed"][1][1], 340)

        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)), expected=dict(zip(keys, exp)))
        got = table.to_dict()
        assert_allclose(got["expected"]["Marl"], 5.2)
        assert_allclose(got["observed"]["Sandstone"], 8)

    def test_str_contingency(self):
        """exercising str(CategoryCounts)"""
        table = CategoryCounts(
            {
                "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
                "b": {"env1": 1, "env3": 1, "env2": 3},
            }
        )
        str(table)
        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)), expected=dict(zip(keys, exp)))
        str(table)

    def test_repr_contingency(self):
        """exercising repr(CategoryCounts) with/without html=True"""
        table = CategoryCounts(
            {
                "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
                "b": {"env1": 1, "env3": 1, "env2": 3},
            }
        )
        str(table)
        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)), expected=dict(zip(keys, exp)))
        _ = table._get_repr_()
        _ = table._get_repr_(html=True)

    def test_accessing_elements(self):
        """successfully access elements"""
        table = CategoryCounts(
            {
                "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
                "b": {"env1": 1, "env3": 1, "env2": 3},
            }
        )
        got = table.observed["rest_of_tree"]["env1"]
        self.assertEqual(got, 2)
        obs = [2, 10, 8, 2, 4]
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)))
        got = table.expected["Clay"]
        assert_allclose(got, 5.2)

    def test_calc_expected(self):
        """expected returns new matrix with expected freqs"""
        matrix = CategoryCounts(
            dict(
                rest_of_tree=dict(env1=2, env3=1, env2=0),
                b=dict(env1=1, env3=1, env2=3),
            )
        )
        assert_allclose(matrix.expected["rest_of_tree"]["env1"], 1.125)
        assert_allclose(matrix.expected["b"]["env1"], 1.875)
        assert_allclose(
            matrix.expected.array.tolist(), [[1.875, 1.875, 1.25], [1.125, 1.125, 0.75]]
        )

    def test_validate_expecteds(self):
        """test provided expecteds total same as observed"""
        with self.assertRaises(AssertionError):
            obs = dict(a=10, b=2, c=2)
            exp = [5, 5, 5]
            CategoryCounts(obs, expected=exp)

    def test_repr_str_html(self):
        """exercising construction of different representations"""
        table = CategoryCounts(
            {
                "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
                "b": {"env1": 1, "env3": 1, "env2": 3},
            }
        )
        got_g1 = table.G_fit()
        got_g2 = table.G_independence()
        got_chisq = table.chisq_test()
        for obj in (table, got_g1, got_g2, got_chisq):
            str(obj)
            repr(obj)
            obj._repr_html_()

    def test_statistics(self):
        """returns TestResult.statistics has stats"""
        table = CategoryCounts(
            {
                "rest_of_tree": {"env1": 2, "env3": 1, "env2": 0},
                "b": {"env1": 1, "env3": 1, "env2": 3},
            }
        )
        got = table.chisq_test()
        stats = got.statistics
        self.assertEqual(stats[0, "pvalue"], got.pvalue)

    def test_calc_expected2(self):
        """handle case where expected is a single column vector"""
        nums = numpy.array([1, 2, 3]).reshape((3, 1))
        got = calc_expected(nums)
        assert_allclose(got, numpy.array([2, 2, 2]).reshape((3, 1)))

    def test_category_counts_from_non_int_arrays(self):
        """handles object and float numpy array, fails if float"""
        a = numpy.array([[31, 36], [58, 138]], dtype=object)
        darr = DictArrayTemplate(["syn", "nsyn"], ["Ts", "Tv"]).wrap(a)
        got = CategoryCounts(darr)
        assert_allclose(got.observed.array.tolist(), a.tolist())

        for dtype in (object, float):
            with self.assertRaises(TypeError):
                a = numpy.array([[31.3, 36], [58, 138]], dtype=dtype)
                darr = DictArrayTemplate(["syn", "nsyn"], ["Ts", "Tv"]).wrap(a)
                _ = CategoryCounts(darr)

        # negative values disallowed
        with self.assertRaises(ValueError):
            a = numpy.array([[31, -36], [58, 138]], dtype=int)
            darr = DictArrayTemplate(["syn", "nsyn"], ["Ts", "Tv"]).wrap(a)
            _ = CategoryCounts(darr)


if __name__ == "__main__":
    main()

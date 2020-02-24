from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3.maths.stats.contingency import CategoryCounts


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
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
        """constrtucted from 2D dict"""
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

    def test_G_ind(self):
        """correctly produce G test of independence"""
        table = CategoryCounts([[762, 327, 468], [484, 239, 477]])
        got = table.G_independence(williams=True)
        self.assertEqual(got.df, 2)

    def test_G_fit_with_expecteds(self):
        """compute G-fit with provided expecteds"""
        obs = [2, 10, 8, 2, 4]
        exp = [5.2] * 5
        keys = ["Marl", "Chalk", "Sandstone", "Clay", "Limestone"]
        table = CategoryCounts(dict(zip(keys, obs)), expected=dict(zip(keys, exp)))

        got = table.G_fit()
        assert_allclose(got.G, 9.849234)
        assert_allclose(got.pvalue, 0.04304536)

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
        for obj in (got_g1, got_g2, got_chisq):
            str(obj)
            repr(obj)
            obj._repr_html_()


if __name__ == "__main__":
    main()

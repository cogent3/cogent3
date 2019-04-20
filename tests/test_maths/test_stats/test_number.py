from unittest import TestCase, main
from cogent3.maths.stats import number
from cogent3 import LoadSeqs
import numpy


class TestNumber(TestCase):
    def test_construction(self):
        nums = number.CategoryCounter('AAAACCCGGGGT')
        self.assertEqual(nums.todict(), dict(A=4, C=3, G=4, T=1))
        self.assertEqual(nums.sum, 12)
        nums['A'] += 1

    def test_valid(self):
        """correctly identify when numbers contains numbers"""
        wrong = number.NumberCounter([0, 'a', 1, 1])
        self.assertFalse(wrong.valid)
        ints = number.NumberCounter([0, 1, 1])
        self.assertTrue(ints.valid)
        floats = number.NumberCounter([0.1, 1., 1.])
        self.assertTrue(floats.valid)
        cmplx = number.NumberCounter([1j, .2j])
        self.assertTrue(cmplx.valid)
        mixed = number.NumberCounter([0.1, 1, 1.1j])
        self.assertTrue(mixed.valid)
        for dtype in (numpy.uint8, numpy.int32, numpy.float64):
            data = numpy.arange(0, 4)
            npy = number.NumberCounter(data.astype(dtype))
            self.assertTrue(npy.valid)

    def test_number_counter_stats(self):
        """stats from NumberCounter correct"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        self.assertEqual(nums.mean, numpy.mean(data))
        self.assertEqual(nums.std, numpy.std(data, ddof=1))
        self.assertEqual(nums.median, numpy.median(data))
        self.assertEqual(nums.quantile(q=.75), numpy.quantile(data, q=.75))
        self.assertEqual(nums.mode, 4)

    def test_category_counter_stats(self):
        """stats from CategoryCounter correct"""
        data = 'TCTTTAGAGAACAGTTTATTATACACTAAA'
        values = [data.count(b) for b in 'ACGT']
        nums = number.CategoryCounter(data)
        self.assertEqual(nums.mean, numpy.mean(values))
        self.assertEqual(nums.std, numpy.std(values, ddof=1))
        self.assertEqual(nums.median, numpy.median(values))
        self.assertEqual(nums.quantile(q=.75), numpy.quantile(values, q=.75))
        self.assertEqual(nums.mode, 'A')
        data = [('T', 'C'), ('T', 'T'), ('T', 'A'), ('G', 'A'), ('G', 'A'),
                ('A', 'C'), ('A', 'G'), ('T', 'T'), ('T', 'A'), ('T', 'T'),
                ('A', 'T'), ('A', 'C'), ('A', 'C'), ('T', 'A'), ('A', 'A'),
                ('A', 'C')]
        values = [1, 3, 3, 2, 4, 1, 1, 1]
        nums = number.CategoryCounter(data)
        self.assertEqual(nums.mean, numpy.mean(values))
        self.assertEqual(nums.std, numpy.std(values, ddof=1))
        self.assertEqual(nums.median, numpy.median(values))
        self.assertEqual(nums.quantile(q=.75), numpy.quantile(values, q=.75))
        self.assertEqual(nums.mode, ('A', 'C'))


if __name__ == '__main__':
    main()

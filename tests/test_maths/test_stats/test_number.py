from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3.maths.stats import number
from cogent3 import LoadSeqs
import numpy


class TestNumber(TestCase):
    def test_construction(self):
        nums = number.CategoryCounter('AAAACCCGGGGT')
        self.assertEqual(nums.todict(), dict(A=4, C=3, G=4, T=1))
        self.assertEqual(nums.sum, 12)
        nums['A'] += 1

    def test_copy(self):
        """copy works"""
        nums = number.CategoryCounter('AAAACCCGGGGT')
        new = nums.copy()
        self.assertNotEqual(id(new), id(nums))
        self.assertEqual(new.todict(), nums.todict())

        nums = number.NumberCounter(data = [0, 0, 2, 4, 4, 4])
        new = nums.copy()
        self.assertNotEqual(id(new), id(nums))
        self.assertEqual(new.todict(), nums.todict())

    def test_construct_from_dict(self):
        """construction from dict of counts"""
        data = {'A': 20, 'Q': 30, 'X': 20}
        got = number.CategoryCounter(data)
        self.assertEqual(got['A'], 20)
        self.assertEqual(got.todict(), data)

    def test_add(self):
        """allow adding elements, or series"""
        nums = number.CategoryCounter('AAAACCCGGGGT')
        nums += 'A'
        self.assertEqual(nums['A'], 5)

    def test_sub(self):
        """allow removing elements"""
        nums = number.CategoryCounter('AAAACCCGGGGT')
        nums -= 'A'
        self.assertEqual(nums['A'], 3)

    def test_to_methods(self):
        """successfully convert to dict, list, array"""
        nums = number.CategoryCounter('AAAACCCGGGGT')
        got = nums.tolist()
        self.assertEqual(got, [4, 3, 4, 1])
        got = nums.tolist(keys='TCAG')
        self.assertEqual(got, [1, 3, 4, 4])
        got = nums.toarray(keys='TCAG')
        assert_allclose(got, numpy.array([1, 3, 4, 4], dtype=int))
        self.assertEqual(nums.todict(), dict(A=4, C=3, G=4, T=1))

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
        self.assertEqual(len(nums), 6)

    def test_category_counter_stats(self):
        """stats from CategoryCounter correct"""
        data = 'TCTTTAGAGAACAGTTTATTATACACTAAA'
        values = [data.count(b) for b in 'ACGT']
        nums = number.CategoryCounter(data)
        self.assertEqual(len(nums), len(data))
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

    def test_usage(self):
        """Alignment.counts_per_seq method correctly applies CategoryCounter"""
        data = {'DogFaced': 'TCATTA', 'FalseVamp': 'TCATTA',
                'FlyingFox': 'TCTTTA', 'FreeTaile': 'TCATTA',
                'Horse': 'TCATTG', 'LeafNose': 'TCTTTA',
                'LittleBro': 'TCATTA', 'Rhino': 'TCATTG',
                'RoundEare': 'TCATTA', 'TombBat': 'TCAGTA'}
        aln = LoadSeqs(data=data, moltype='dna')
        got = aln.counts_per_pos(motif_length=3)
        self.assertEqual(got[0, 'TCA'], 8)
        self.assertEqual(got[0, 'TCT'], 2)
        self.assertEqual(got[1, 'TTA'], 7)
        self.assertEqual(got[1, 'GTA'], 1)
        self.assertEqual(got[1, 'TTG'], 2)

    def test_numbers_update(self):
        """correctly update number counts"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        data = [2, 4, 4, 4, 6, 5]
        nums2 = number.NumberCounter(data)
        nums.update_from_counts(nums2)
        self.assertEqual(nums[2], 2)
        self.assertEqual(nums[4], 6)
        self.assertEqual(nums[1], 0)

        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        nums.update_from_series([2, 4, 4, 4, 6, 5])
        self.assertEqual(nums[2], 2)
        self.assertEqual(nums[4], 6)
        self.assertEqual(nums[1], 0)

        with self.assertRaises(TypeError):
            counts = number.CategoryCounter('AAAACCCGGGGT')
            nums.update_from_counts(counts)



if __name__ == '__main__':
    main()

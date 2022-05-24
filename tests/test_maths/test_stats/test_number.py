from unittest import TestCase, main

import numpy

from numpy.testing import assert_allclose

from cogent3 import make_aligned_seqs
from cogent3.maths.stats import number


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestNumber(TestCase):
    def test_construction(self):
        nums = number.CategoryCounter("AAAACCCGGGGT")
        self.assertEqual(nums.to_dict(), dict(A=4, C=3, G=4, T=1))
        self.assertEqual(nums.sum, 12)
        nums["A"] += 1

    def test_copy(self):
        """copy works"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        new = nums.copy()
        self.assertNotEqual(id(new), id(nums))
        self.assertEqual(new.to_dict(), nums.to_dict())

        nums = number.NumberCounter(data=[0, 0, 2, 4, 4, 4])
        new = nums.copy()
        self.assertNotEqual(id(new), id(nums))
        self.assertEqual(new.to_dict(), nums.to_dict())

    def test_construct_from_dict(self):
        """construction from dict of counts"""
        data = {"A": 20, "Q": 30, "X": 20}
        got = number.CategoryCounter(data)
        self.assertEqual(got["A"], 20)
        self.assertEqual(got.to_dict(), data)

    def test_add(self):
        """allow adding elements, or series"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        nums += "A"
        self.assertEqual(nums["A"], 5)

    def test_sub(self):
        """allow removing elements"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        nums -= "A"
        self.assertEqual(nums["A"], 3)

    def test_to_methods(self):
        """successfully convert to dict, list, array"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        got = nums.tolist()
        self.assertEqual(got, [4, 3, 4, 1])
        got = nums.tolist(keys="TCAG")
        self.assertEqual(got, [1, 3, 4, 4])
        got = nums.to_array(keys="TCAG")
        assert_allclose(got, numpy.array([1, 3, 4, 4], dtype=int))
        self.assertEqual(nums.to_dict(), dict(A=4, C=3, G=4, T=1))

    def test_to_table(self):
        """produces correct Table structure"""
        data = [
            ("Ovary-AdenoCA", "IGR"),
            ("Liver-HCC", "Intron"),
            ("Panc-AdenoCA", "Intron"),
            ("Panc-AdenoCA", "Intron"),
        ]
        nums = number.CategoryCounter(data)
        t = nums.to_table(column_names=None, title="blah")
        self.assertEqual(t.header, ("key", "count"))
        # if the key is a tuple, then the unexpanded column values are also
        self.assertIsInstance(t[0, 0], tuple)
        self.assertEqual(t.title, "blah")
        # you can use any data type as a key, but Table column is a str
        t = nums.to_table(column_names=2)
        self.assertEqual(t.header, ("2", "count"))
        t = nums.to_table(column_names="blah")
        self.assertEqual(t.header, ("blah", "count"))
        t = nums.to_table(column_names=["A", "B"])
        self.assertEqual(t.header, ("A", "B", "count"))

        with self.assertRaises(AssertionError):
            # key does not have 3 dimensions
            _ = nums.to_table(column_names=["A", "B", "C"])

        with self.assertRaises(AssertionError):
            # key does not have 1 dimension
            _ = nums.to_table(column_names=[1])

    def test_valid(self):
        """correctly identify when numbers contains numbers"""
        wrong = number.NumberCounter([0, "a", 1, 1])
        self.assertFalse(wrong.valid)
        ints = number.NumberCounter([0, 1, 1])
        self.assertTrue(ints.valid)
        floats = number.NumberCounter([0.1, 1.0, 1.0])
        self.assertTrue(floats.valid)
        cmplx = number.NumberCounter([1j, 0.2j])
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
        self.assertEqual(nums.quantile(q=0.75), numpy.quantile(data, q=0.75))
        self.assertEqual(nums.mode, 4)
        self.assertEqual(len(nums), 6)

    def test_keys_values_items(self):
        """return a list of these elements"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.CategoryCounter(data)
        self.assertEqual(nums.keys(), [0, 2, 4])
        self.assertEqual(nums.values(), [2, 1, 3])
        self.assertEqual(nums.items(), [(0, 2), (2, 1), (4, 3)])

        freqs = nums.to_freqs()
        self.assertEqual(freqs.keys(), [0, 2, 4])
        assert_allclose(freqs.values(), [0.3333333333333333, 0.16666666666666666, 0.5])
        self.assertEqual(len(freqs.items()), 3)
        self.assertEqual(freqs.items()[-1], (4, 0.5))

    def test_repr(self):
        """should precede with class name"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.CategoryCounter(data)
        got = repr(nums)
        self.assertTrue(got.startswith(nums.__class__.__name__))
        freqs = nums.to_freqs()
        got = repr(freqs)
        self.assertTrue(got.startswith(freqs.__class__.__name__))
        nums = number.NumberCounter(data)
        got = repr(nums)
        self.assertTrue(got.startswith(nums.__class__.__name__))

    def test_category_counter_stats(self):
        """stats from CategoryCounter correct"""
        data = "TCTTTAGAGAACAGTTTATTATACACTAAA"
        values = [data.count(b) for b in "ACGT"]
        nums = number.CategoryCounter(data)
        self.assertEqual(len(nums), len(data))
        self.assertEqual(nums.mean, numpy.mean(values))
        self.assertEqual(nums.std, numpy.std(values, ddof=1))
        self.assertEqual(nums.median, numpy.median(values))
        self.assertEqual(nums.quantile(q=0.75), numpy.quantile(values, q=0.75))
        self.assertEqual(nums.mode, "A")
        data = [
            ("T", "C"),
            ("T", "T"),
            ("T", "A"),
            ("G", "A"),
            ("G", "A"),
            ("A", "C"),
            ("A", "G"),
            ("T", "T"),
            ("T", "A"),
            ("T", "T"),
            ("A", "T"),
            ("A", "C"),
            ("A", "C"),
            ("T", "A"),
            ("A", "A"),
            ("A", "C"),
        ]
        values = [1, 3, 3, 2, 4, 1, 1, 1]
        nums = number.CategoryCounter(data)
        self.assertEqual(nums.mean, numpy.mean(values))
        self.assertEqual(nums.std, numpy.std(values, ddof=1))
        self.assertEqual(nums.median, numpy.median(values))
        self.assertEqual(nums.quantile(q=0.75), numpy.quantile(values, q=0.75))
        self.assertEqual(nums.mode, ("A", "C"))

    def test_usage(self):
        """Alignment.counts_per_seq method correctly applies CategoryCounter"""
        data = {
            "DogFaced": "TCATTA",
            "FalseVamp": "TCATTA",
            "FlyingFox": "TCTTTA",
            "FreeTaile": "TCATTA",
            "Horse": "TCATTG",
            "LeafNose": "TCTTTA",
            "LittleBro": "TCATTA",
            "Rhino": "TCATTG",
            "RoundEare": "TCATTA",
            "TombBat": "TCAGTA",
        }
        aln = make_aligned_seqs(data=data, moltype="dna")
        got = aln.counts_per_pos(motif_length=3)
        self.assertEqual(got[0, "TCA"], 8)
        self.assertEqual(got[0, "TCT"], 2)
        self.assertEqual(got[1, "TTA"], 7)
        self.assertEqual(got[1, "GTA"], 1)
        self.assertEqual(got[1, "TTG"], 2)

    def test_entropy(self):
        """CategoryCounter correctly calculates entropy"""
        freqs = numpy.array([4 / 12, 3 / 12, 4 / 12, 1 / 12])
        expect = -(freqs * numpy.log2(freqs)).sum()
        nums = number.CategoryCounter("AAAACCCGGGGT")
        assert_allclose(nums.entropy, expect)
        nums = number.CategoryCounter("AAAA")
        assert_allclose(nums.entropy, 0)

    def test_to_freqs(self):
        """CategoryCounter.to_freqs produces CategoryFreqs"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        freqs = nums.to_freqs()
        assert_allclose(freqs.to_array(list(freqs)), nums.to_array(list(freqs)) / 12)

    def test_expand(self):
        """correctly reconstitutes original series content"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        expanded = nums.expand()
        self.assertEqual(expanded, list("AAAACCCGGGGT"))

    def test_categoryfreqs_entropy(self):
        """correctly returns frequencies"""
        vals = numpy.array([4 / 12, 3 / 12, 4 / 12, 1 / 12])
        expect = -(vals * numpy.log2(vals)).sum()
        freqs = number.CategoryFreqs({"A": 4, "C": 3, "G": 4, "T": 1}, total=12)
        assert_allclose(freqs.entropy, expect)

    def test_to_normalized(self):
        """correctly recalibrate CategoryFreqs"""
        freqs = number.CategoryFreqs({"A": 4, "C": 2, "G": 4}, total=12)
        self.assertEqual(freqs["A"], 4 / 12)
        freqs = freqs.to_normalized()
        self.assertEqual(freqs["A"], 4 / 10)

        # from an empty dict
        freqs = number.CategoryFreqs()
        d = freqs.to_normalized()
        self.assertEqual(d.to_dict(), {})

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
            counts = number.CategoryCounter("AAAACCCGGGGT")
            nums.update_from_counts(counts)

    def test_count(self):
        """correctly counts across key dimensions"""
        data = [
            ("T", "C"),
            ("T", "T"),
            ("T", "A"),
            ("G", "A"),
            ("G", "A"),
            ("A", "C"),
            ("A", "G"),
            ("T", "T"),
            ("T", "A"),
            ("T", "T"),
            ("A", "T"),
            ("A", "C"),
            ("A", "C"),
            ("T", "A"),
            ("A", "A"),
            ("A", "C"),
        ]
        nums = number.CategoryCounter(data)
        i0 = nums.count(0)
        self.assertEqual(i0["T"], 7)
        self.assertEqual(i0["G"], 2)
        self.assertEqual(i0["A"], 7)
        self.assertEqual(i0["C"], 0)
        # works same if keys are strings
        nums = number.CategoryCounter(["".join(e) for e in data])
        i0 = nums.count(0)
        self.assertEqual(i0["T"], 7)
        self.assertEqual(i0["G"], 2)
        self.assertEqual(i0["A"], 7)
        self.assertEqual(i0["C"], 0)
        with self.assertRaises(IndexError):
            _ = nums.count([0, 3])

        i0 = nums.count([0])
        self.assertEqual(i0["G"], 2)
        with self.assertRaises(IndexError):
            _ = nums.count([0, 3])
        i1 = nums.count(1)
        self.assertEqual(i1["A"], 6)
        self.assertEqual(i1["C"], 5)
        self.assertEqual(i1["T"], 4)

        data = {
            ("A", "C", "G"): 10,
            ("A", "T", "G"): 4,
            ("C", "C", "G"): 3,
            ("G", "C", "G"): 6,
        }
        nums = number.CategoryCounter(data)
        i02 = nums.count([0, 2])
        self.assertEqual(i02[("A", "G")], 14)
        self.assertEqual(i02[("C", "G")], 3)
        self.assertEqual(i02[("G", "G")], 6)

    def test_to_dictarray(self):
        """correctly constructs dict arrays"""
        d1 = {"T": 87, "C": 81, "A": 142, "expect": [142, 81, 87]}
        d2 = {
            ("T", "G"): 87,
            ("C", "C"): 81,
            ("A", "G"): 142,
            ("T", "T"): 58,
            "expect": [[0, 142, 0], [81, 0, 0], [0, 87, 58]],
        }
        d3 = {
            ("T", "G", "A"): 87,
            ("C", "C", "C"): 81,
            ("A", "G", "A"): 142,
            ("T", "T", "C"): 58,
            "expect": [
                [[0, 0], [142, 0], [0, 0]],
                [[0, 81], [0, 0], [0, 0]],
                [[0, 0], [87, 0], [0, 58]],
            ],
        }
        for d in (d1, d2, d3):
            expect = d.pop("expect")
            cat_count = number.CategoryCounter(d)
            darr = cat_count.to_dictarray()
            assert_allclose(darr.array, expect)

    def test_to_categorical(self):
        """correctly constructs categorical data"""
        d1 = {"T": 87, "C": 81, "A": 142, "expect": [142, 81, 87]}
        d2 = {
            ("T", "G"): 87,
            ("C", "C"): 81,
            ("A", "G"): 142,
            ("T", "T"): 58,
            "expect": [[0, 142, 0], [81, 0, 0], [0, 87, 58]],
        }
        d3 = {
            ("T", "G", "A"): 87,
            ("C", "C", "C"): 81,
            ("A", "G", "A"): 142,
            ("T", "T", "C"): 58,
        }
        for d in (d1, d2):
            expect = d.pop("expect")
            cats = number.CategoryCounter(d)
            cat_count = cats.to_categorical()
            assert_allclose(cat_count.observed.array, expect, err_msg=d)

        with self.assertRaises(NotImplementedError):
            cats = number.CategoryCounter(d3)
            cats.to_categorical()


if __name__ == "__main__":
    main()

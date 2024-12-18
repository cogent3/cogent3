from unittest import TestCase

import numpy
import pytest
from numpy.testing import assert_allclose

from cogent3 import make_aligned_seqs
from cogent3.maths.stats import number


class TestNumber(TestCase):
    def test_construction(self):
        nums = number.CategoryCounter("AAAACCCGGGGT")
        assert nums.to_dict() == {"A": 4, "C": 3, "G": 4, "T": 1}
        assert nums.sum == 12
        nums["A"] += 1

    def test_copy(self):
        """copy works"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        new = nums.copy()
        assert id(new) != id(nums)
        assert new.to_dict() == nums.to_dict()

        nums = number.NumberCounter(data=[0, 0, 2, 4, 4, 4])
        new = nums.copy()
        assert id(new) != id(nums)
        assert new.to_dict() == nums.to_dict()

    def test_construct_from_dict(self):
        """construction from dict of counts"""
        data = {"A": 20, "Q": 30, "X": 20}
        got = number.CategoryCounter(data)
        assert got["A"] == 20
        assert got.to_dict() == data

    def test_add(self):
        """allow adding elements, or series"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        nums += "A"
        assert nums["A"] == 5

    def test_sub(self):
        """allow removing elements"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        nums -= "A"
        assert nums["A"] == 3

    def test_to_methods(self):
        """successfully convert to dict, list, array"""
        nums = number.CategoryCounter("AAAACCCGGGGT")
        got = nums.tolist()
        assert got == [4, 3, 4, 1]
        got = nums.tolist(keys="TCAG")
        assert got == [1, 3, 4, 4]
        got = nums.to_array(keys="TCAG")
        assert_allclose(got, numpy.array([1, 3, 4, 4], dtype=int))
        assert nums.to_dict() == {"A": 4, "C": 3, "G": 4, "T": 1}

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
        assert t.header == ("key", "count")
        # if the key is a tuple, then the unexpanded column values are also
        assert isinstance(t[0, 0], tuple)
        assert t.title == "blah"
        # you can use any data type as a key, but Table column is a str
        t = nums.to_table(column_names=2)
        assert t.header == ("2", "count")
        t = nums.to_table(column_names="blah")
        assert t.header == ("blah", "count")
        t = nums.to_table(column_names=["A", "B"])
        assert t.header == ("A", "B", "count")

        with pytest.raises(AssertionError):
            # key does not have 3 dimensions
            _ = nums.to_table(column_names=["A", "B", "C"])

        with pytest.raises(AssertionError):
            # key does not have 1 dimension
            _ = nums.to_table(column_names=[1])

    def test_valid(self):
        """correctly identify when numbers contains numbers"""
        wrong = number.NumberCounter([0, "a", 1, 1])
        assert not wrong.valid
        ints = number.NumberCounter([0, 1, 1])
        assert ints.valid
        floats = number.NumberCounter([0.1, 1.0, 1.0])
        assert floats.valid
        cmplx = number.NumberCounter([1j, 0.2j])
        assert cmplx.valid
        mixed = number.NumberCounter([0.1, 1, 1.1j])
        assert mixed.valid
        for dtype in (numpy.uint8, numpy.int32, numpy.float64):
            data = numpy.arange(0, 4)
            npy = number.NumberCounter(data.astype(dtype))
            assert npy.valid

    def test_number_counter_stats(self):
        """stats from NumberCounter correct"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        assert nums.mean == numpy.mean(data)
        assert nums.std == numpy.std(data, ddof=1)
        assert nums.median == numpy.median(data)
        assert nums.quantile(q=0.75) == numpy.quantile(data, q=0.75)
        assert nums.mode == 4
        assert len(nums) == 6

    def test_keys_values_items(self):
        """return a list of these elements"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.CategoryCounter(data)
        assert nums.keys() == [0, 2, 4]
        assert nums.values() == [2, 1, 3]
        assert nums.items() == [(0, 2), (2, 1), (4, 3)]

        freqs = nums.to_freqs()
        assert freqs.keys() == [0, 2, 4]
        assert_allclose(freqs.values(), [0.3333333333333333, 0.16666666666666666, 0.5])
        assert len(freqs.items()) == 3
        assert freqs.items()[-1] == (4, 0.5)

    def test_repr(self):
        """should precede with class name"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.CategoryCounter(data)
        got = repr(nums)
        assert got.startswith(nums.__class__.__name__)
        freqs = nums.to_freqs()
        got = repr(freqs)
        assert got.startswith(freqs.__class__.__name__)
        nums = number.NumberCounter(data)
        got = repr(nums)
        assert got.startswith(nums.__class__.__name__)

    def test_category_counter_stats(self):
        """stats from CategoryCounter correct"""
        data = "TCTTTAGAGAACAGTTTATTATACACTAAA"
        values = [data.count(b) for b in "ACGT"]
        nums = number.CategoryCounter(data)
        assert len(nums) == len(data)
        assert nums.mean == numpy.mean(values)
        assert nums.std == numpy.std(values, ddof=1)
        assert nums.median == numpy.median(values)
        assert nums.quantile(q=0.75) == numpy.quantile(values, q=0.75)
        assert nums.mode == "A"
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
        assert nums.mean == numpy.mean(values)
        assert nums.std == numpy.std(values, ddof=1)
        assert nums.median == numpy.median(values)
        assert nums.quantile(q=0.75) == numpy.quantile(values, q=0.75)
        assert nums.mode == ("A", "C")

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
        assert got[0, "TCA"] == 8
        assert got[0, "TCT"] == 2
        assert got[1, "TTA"] == 7
        assert got[1, "GTA"] == 1
        assert got[1, "TTG"] == 2

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
        assert expanded == list("AAAACCCGGGGT")

    def test_categoryfreqs_entropy(self):
        """correctly returns frequencies"""
        vals = numpy.array([4 / 12, 3 / 12, 4 / 12, 1 / 12])
        expect = -(vals * numpy.log2(vals)).sum()
        freqs = number.CategoryFreqs({"A": 4, "C": 3, "G": 4, "T": 1}, total=12)
        assert_allclose(freqs.entropy, expect)

    def test_to_normalized(self):
        """correctly recalibrate CategoryFreqs"""
        freqs = number.CategoryFreqs({"A": 4, "C": 2, "G": 4}, total=12)
        assert freqs["A"] == 4 / 12
        freqs = freqs.to_normalized()
        assert freqs["A"] == 4 / 10

        # from an empty dict
        freqs = number.CategoryFreqs()
        d = freqs.to_normalized()
        assert d.to_dict() == {}

    def test_numbers_update(self):
        """correctly update number counts"""
        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        data = [2, 4, 4, 4, 6, 5]
        nums2 = number.NumberCounter(data)
        nums.update_from_counts(nums2)
        assert nums[2] == 2
        assert nums[4] == 6
        assert nums[1] == 0

        data = [0, 0, 2, 4, 4, 4]
        nums = number.NumberCounter(data)
        nums.update_from_series([2, 4, 4, 4, 6, 5])
        assert nums[2] == 2
        assert nums[4] == 6
        assert nums[1] == 0

        with pytest.raises(TypeError):
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
        assert i0["T"] == 7
        assert i0["G"] == 2
        assert i0["A"] == 7
        assert i0["C"] == 0
        # works same if keys are strings
        nums = number.CategoryCounter(["".join(e) for e in data])
        i0 = nums.count(0)
        assert i0["T"] == 7
        assert i0["G"] == 2
        assert i0["A"] == 7
        assert i0["C"] == 0
        with pytest.raises(IndexError):
            _ = nums.count([0, 3])

        i0 = nums.count([0])
        assert i0["G"] == 2
        with pytest.raises(IndexError):
            _ = nums.count([0, 3])
        i1 = nums.count(1)
        assert i1["A"] == 6
        assert i1["C"] == 5
        assert i1["T"] == 4

        data = {
            ("A", "C", "G"): 10,
            ("A", "T", "G"): 4,
            ("C", "C", "G"): 3,
            ("G", "C", "G"): 6,
        }
        nums = number.CategoryCounter(data)
        i02 = nums.count([0, 2])
        assert i02["A", "G"] == 14
        assert i02["C", "G"] == 3
        assert i02["G", "G"] == 6

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

        with pytest.raises(NotImplementedError):
            cats = number.CategoryCounter(d3)
            cats.to_categorical()

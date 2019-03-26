from unittest import TestCase, main
import numpy
from numpy.testing import assert_allclose
from cogent3 import DNA
from cogent3.util.dict_array import DictArrayTemplate


class DictArrayTest(TestCase):
    a = numpy.identity(3, int)

    def test_construct_both_dim_str(self):
        b = DictArrayTemplate('abc', 'ABC').wrap(self.a)
        self.assertEqual(b[0].array.tolist(), [1, 0, 0])
        self.assertEqual(b['a'].array.tolist(), [1, 0, 0])
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_keye_levels(self):
        """DictArray both levels have keys."""
        b = DictArrayTemplate('abc', 'ABC').wrap(self.a)
        self.assertEqual(b.keys(), ['a', 'b', 'c'])
        self.assertEqual(b['a'].keys(), ['A', 'B', 'C'])
        self.assertEqual(list(b['a']), [1, 0, 0])
        self.assertEqual(sum(b['a']), 1)

    def test_int_labels(self):
        """DictArray with no labels."""
        b = DictArrayTemplate(3, 3).wrap(self.a)
        self.assertEqual(b.keys(), [0, 1, 2])
        self.assertEqual(b[0].keys(), [0, 1, 2])
        self.assertEqual(sum(b[0]), 1)

    def test_with_mixed_label_types(self):
        """DictArray constructed with mixed label types."""
        b = DictArrayTemplate('ABC', 3).wrap(self.a)
        self.assertEqual(b.keys(), ['A', 'B', 'C'])
        self.assertEqual(b['A'].keys(), [0, 1, 2])

    def test_numpy_ops(self):
        """DictArray should work properly in numpy operations."""
        darr = DictArrayTemplate(list(DNA), list(DNA)).wrap([[.7, .1, .1, .1],
                                                             [.1, .7, .1, .1],
                                                             [.1, .1, .7, .1],
                                                             [.1, .1, .1, .7]])
        mprobs = numpy.array([0.25, 0.25, 0.25, 0.25])
        assert_allclose(mprobs.dot(darr), [0.25, 0.25, 0.25, 0.25])
        assert_allclose(numpy.dot(mprobs, darr), [0.25, 0.25, 0.25, 0.25])

    def test_asdict(self):
        """DictArray.asdict() should convert nested DictArray instances to
        dict's too."""
        a = numpy.identity(3, int)
        b = DictArrayTemplate('abc', 'ABC').wrap(a)
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        c = DictArrayTemplate('de', 'DE').wrap([[b, b], [b, b]])
        self.assertTrue(isinstance(c.asdict()['d'], dict))


if __name__ == '__main__':
    main()

from unittest import TestCase, main
import numpy
from cogent3 import DNA
from cogent3.util.dict_array import DictArrayTemplate


class DictArrayTest(TestCase):

    def test_construct_both_dim_str(self):
        a = numpy.identity(3, int)
        b = DictArrayTemplate('abc', 'ABC').wrap(a)
        self.assertEqual(b[0].array.tolist(), [1, 0, 0])
        self.assertEqual(b['a'].array.tolist(), [1, 0, 0])
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        # checks both levels have keys.
        self.assertEqual(b.keys(), ['a', 'b', 'c'])
        self.assertEqual(b['a'].keys(), ['A', 'B', 'C'])
        self.assertEqual(list(b['a']), [1, 0, 0])
        self.assertEqual(sum(b['a']), 1)

        # tests for no labels.
        b = DictArrayTemplate(3, 3).wrap(a)
        self.assertEqual(b.keys(), [0, 1, 2])
        self.assertEqual(b[0].keys(), [0, 1, 2])
        self.assertEqual(sum(b[0]), 1)

        # tests for mixed labels.
        b = DictArrayTemplate('ABC', 3).wrap(a)
        self.assertEqual(b.keys(), ['A', 'B', 'C'])
        self.assertEqual(b['A'].keys(), [0, 1, 2])

        # ``DictArray`` should work properly in ``numpy`` operations.
        darr = DictArrayTemplate(list(DNA), list(DNA)).wrap([[.7, .1, .1, .1],
                                                             [.1, .7, .1, .1],
                                                             [.1, .1, .7, .1],
                                                             [.1, .1, .1, .7]])
        mprobs = numpy.array([0.25, 0.25, 0.25, 0.25])
        self.assertEqual(mprobs.dot(darr).tolist(), [0.25, 0.25, 0.25, 0.25])
        self.assertEqual(numpy.dot(mprobs, darr).tolist(), [0.25, 0.25, 0.25, 0.25])

        # ``DictArray.asdict()`` should convert nested ``DictArray`` instances to dict's too.
        darr = DictArrayTemplate('A', 'C')
        a = numpy.identity(3, int)
        b = DictArrayTemplate('abc', 'ABC').wrap(a)
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        c = DictArrayTemplate('de', 'DE').wrap([[b, b], [b, b]])
        self.assertTrue(isinstance(c.asdict()['d'], dict))


if __name__ == '__main__':
    main()

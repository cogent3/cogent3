from unittest import TestCase, main
import numpy
from numpy.testing import assert_allclose
from cogent3 import DNA
from cogent3.util.dict_array import (DictArrayTemplate, DictArray,
                                     convert_1D_dict, convert2Ddistance,
                                     convert2DDict, convert_dict,
                                     convert_series, convert_for_dictarray)


class DictArrayTest(TestCase):
    a = numpy.identity(3, int)

    def test_convert_series(self):
        """convert_series produces valid template input"""
        vals, row_keys, col_keys = convert_series([[4], [5]],
                                                  ['A', 'B'],
                                                  ['a'])
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4], [5]])
        data = [[245, 599]]
        vals, row_keys, col_keys = convert_series(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), data)

        vals, row_keys, col_keys = convert_series(data[0])
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), data[0])

    def test_convert_dict(self):
        """convert_dict produces valid template input"""
        twoDdict = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert_dict(twoDdict)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4], [5]])

    def test_convert2DDict(self):
        """convert2DDict produces valid template input"""
        data = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert2DDict(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4], [5]])
        # row keys, then column
        self.assertEqual(b.template.names, [['b', 'c'], ['a']])

        data = {'a': {'a': 0, 'b': 1, 'e': 0},
                'b': {'a': 1, 'b': 0, 'e': 4},
                'e': {'a': 0, 'b': 4, 'e': 0}}
        vals, row_keys, col_keys = convert2DDict(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        got = b.asdict()
        self.assertEqual(got, data)
        self.assertEqual(b.template.names,
                         [['a', 'b', 'e'], ['a', 'b', 'e']])

        data = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert2DDict(data, make_symmetric=True)
        self.assertEqual(row_keys, col_keys)
        self.assertEqual(vals, [[0, 4, 5], [4, 0, 0], [5, 0, 0]])

    def test_convert2Ddistance(self):
        """convert2Ddistance produces valid template input"""
        data = {('a', 'b'): 4, ('a', 'c'): 5}
        vals, row_keys, col_keys = convert2Ddistance(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[0, 4, 5], [4, 0, 0], [5, 0, 0]])

    def test_convert_1D_dict(self):
        """convert_1D_dict produces valid template input"""
        data = dict(a=0, b=35, c=45)
        vals, keys = convert_1D_dict(data)
        b = DictArrayTemplate(keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [0, 35, 45])

    def test_construct_both_dim_str(self):
        """correctly construct when keys for both dimensions are str"""
        b = DictArrayTemplate('abc', 'ABC').wrap(self.a)
        self.assertEqual(b[0].array.tolist(), [1, 0, 0])
        self.assertEqual(b['a'].array.tolist(), [1, 0, 0])
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_key_levels(self):
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
        b = DictArrayTemplate('abc', 'ABC')
        b = b.wrap(a)
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        c = DictArrayTemplate('de', 'DE').wrap([[b, b], [b, b]])
        self.assertTrue(isinstance(c.asdict()['d'], dict))

    def test_convert_for_dictarray(self):
        """successfully delegates when constructed from a DictArray"""
        a = numpy.identity(3, int)
        b = DictArrayTemplate('abc', 'ABC').wrap(a)
        vals, row_keys, col_keys = convert_for_dictarray(b)
        got = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(got.array.tolist(), b.array.tolist())
        # the wrap method creates a new array
        self.assertIsNot(got.array, b.array)

    def test_convert_for_dictarray(self):
        """convert_for_dictarray correctly delegates"""
        b = DictArrayTemplate('abc', 'ABC').wrap(self.a)
        data_types = ([[245, 599]],
                      dict(a=dict(b=4, c=5)),
                      {('a', 'b'): 4, ('a', 'c'): 5},
                      dict(a=0, b=35, c=45),
                      b)
        for data in data_types:
            vals, row_keys, col_keys = convert_for_dictarray(data)
            _ = DictArrayTemplate(row_keys, col_keys).wrap(vals)

    def test_direct_construction(self):
        """directly construct a dict array"""
        b = DictArrayTemplate('abc', 'ABC').wrap(self.a)
        data_types = ([[245, 599]],
                      dict(a=dict(b=4, c=5)),
                      {('a', 'b'): 4, ('a', 'c'): 5},
                      dict(a=0, b=35, c=45),
                      b)
        for data in data_types:
            g = DictArray(data)


if __name__ == '__main__':
    main()

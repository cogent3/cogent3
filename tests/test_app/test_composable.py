from unittest import TestCase, main
from cogent3.app.composable import ComposableSeq, ErrorResult

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestComposableBase(TestCase):
    def test_composable(self):
        """correctly form string"""
        aseqfunc1 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        aseqfunc2 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        comb = aseqfunc1 + aseqfunc2
        expect = ("ComposableSeq(type='sequences') + "
                  "ComposableSeq(type='sequences')")
        got = str(comb)
        self.assertEqual(got, expect)

class TestErrorResult(TestCase):
    def test_err_result(self):
        """excercise creation of ErrorResult"""
        result = ErrorResult('SKIP', 'this', 'some obj')
        self.assertFalse(result)
        self.assertEqual(result.origin, 'this')
        self.assertEqual(result.message, 'some obj')

        try:
            _ = 0
            raise ValueError("error message")
        except ValueError as err:
            result = ErrorResult('SKIP', 'this', err.args[0])

        self.assertEqual(result.message, 'error message')


if __name__ == "__main__":
    main()

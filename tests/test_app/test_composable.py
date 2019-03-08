from unittest import TestCase, main
from cogent3.app.composable import ComposableSeq, ErrorResult
from cogent3.app.translate import select_translatable
from cogent3.app.sample import omit_degenerates, min_length

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

    def test_str(self):
        """str representation correctly represents parameterisations"""
        func = select_translatable()
        got = str(func)
        self.assertEqual(got, "select_translatable(type='sequences', "
                              "moltype='dna', gc='Standard Nuclear', "
                              "allow_rc=False, trim_terminal_stop=True)")

        func = select_translatable(allow_rc=True)
        got = str(func)
        self.assertEqual(got, "select_translatable(type='sequences', "
                              "moltype='dna', gc='Standard Nuclear', "
                              "allow_rc=True, trim_terminal_stop=True)")


        nodegen = omit_degenerates()
        got = str(nodegen)
        self.assertEqual(got, "omit_degenerates(type='aligned', moltype=None, "
                              "gap_is_degen=True, motif_length=1)")
        ml = min_length(100)
        got = str(ml)
        self.assertEqual(got, "min_length(type='sequences', length=100, "
                              "motif_length=1, subtract_degen=True, "
                              "moltype=None)")

if __name__ == "__main__":
    main()

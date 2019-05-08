from unittest import TestCase, main
from tempfile import TemporaryDirectory
from cogent3.app.composable import ComposableSeq, NotCompletedResult
from cogent3.app.translate import select_translatable
from cogent3.app.sample import omit_degenerates, min_length
from cogent3.app.tree import quick_tree
from cogent3.app import io as io_app, sample as sample_app
from cogent3.core.alignment import ArrayAlignment

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestCheckpoint(TestCase):
    def test_checkpointable(self):
        """chained funcs should be be able to apply a checkpoint"""
        path = 'data/brca1.fasta'
        reader = io_app.load_aligned(moltype='dna')
        omit_degens = sample_app.omit_degenerates(moltype='dna')
        with TemporaryDirectory(dir='.') as dirname:
            writer = io_app.write_seqs(dirname)
            aln = reader(path)
            outpath = writer(aln)

            read_write = reader + writer
            got = read_write(path)  # should skip reading and return path
            self.assertEqual(got, outpath)
            read_write.disconnect()  # allows us to reuse bits
            read_write_degen = reader + writer + omit_degens
            got = read_write_degen(path)  # should return an alignment instance
            self.assertIsInstance(got, ArrayAlignment)
            self.assertTrue(len(got) > 1000)


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

    def test_composables_once(self):
        """composables can only be used in a single composition"""
        aseqfunc1 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        aseqfunc2 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        comb = aseqfunc1 + aseqfunc2
        with self.assertRaises(AssertionError):
            aseqfunc3 = ComposableSeq(
                input_type='sequences', output_type='sequences')
            comb2 = aseqfunc1 + aseqfunc3
        # the other order
        with self.assertRaises(AssertionError):
            aseqfunc3 = ComposableSeq(
                input_type='sequences', output_type='sequences')
            comb2 = aseqfunc3 + aseqfunc2

    def test_disconnect(self):
        """disconnect breaks all connections and allows parts to be reused"""
        aseqfunc1 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        aseqfunc2 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        aseqfunc3 = ComposableSeq(
            input_type='sequences', output_type='sequences')
        comb = aseqfunc1 + aseqfunc2 + aseqfunc3
        comb.disconnect()
        self.assertEqual(aseqfunc1.input, None)
        self.assertEqual(aseqfunc1.output, None)
        self.assertEqual(aseqfunc3.input, None)
        self.assertEqual(aseqfunc3.output, None)
        # should be able to compose a new one now
        comb2 = aseqfunc1 + aseqfunc3


class TestNotCompletedResult(TestCase):
    def test_err_result(self):
        """excercise creation of NotCompletedResult"""
        result = NotCompletedResult('SKIP', 'this', 'some obj')
        self.assertFalse(result)
        self.assertEqual(result.origin, 'this')
        self.assertEqual(result.message, 'some obj')

        try:
            _ = 0
            raise ValueError("error message")
        except ValueError as err:
            result = NotCompletedResult('SKIP', 'this', err.args[0])

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

        qt = quick_tree()
        self.assertEqual(str(qt), "quick_tree(type='tree', distance='TN93',"
                                  " moltype='dna')")


if __name__ == "__main__":
    main()

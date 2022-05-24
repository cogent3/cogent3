import os
import pathlib

from pickle import dumps, loads
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from unittest.mock import Mock

from scitrack import CachingLogger

from cogent3.app import io as io_app
from cogent3.app import sample as sample_app
from cogent3.app.composable import (
    SERIALISABLE_TYPE,
    ComposableSeq,
    NotCompleted,
    appify,
    user_function,
)
from cogent3.app.sample import min_length, omit_degenerates
from cogent3.app.translate import select_translatable
from cogent3.app.tree import quick_tree
from cogent3.core.alignment import ArrayAlignment


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestCheckpoint(TestCase):
    def test_checkpointable(self):
        """chained funcs should be be able to apply a checkpoint"""
        reader = io_app.load_aligned(moltype="dna")
        omit_degens = sample_app.omit_degenerates(moltype="dna")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_seqs(dirname)
            path = "data" + os.sep + "brca1.fasta"
            aln = reader(path)
            outpath = writer(aln)

            read_write = reader + writer
            got = read_write(path)  # should skip reading and return path
            self.assertTrue(got.endswith(outpath))
            read_write.disconnect()  # allows us to reuse bits
            read_write_degen = reader + writer + omit_degens
            # should return an alignment instance
            got = read_write_degen(path)
            self.assertIsInstance(got, ArrayAlignment)
            self.assertTrue(len(got) > 1000)


ComposableSeq._input_types = ComposableSeq._output_types = {None}


class TestComposableBase(TestCase):
    def test_composable(self):
        """correctly form string"""
        aseqfunc1 = ComposableSeq(input_types="sequences", output_types="sequences")
        aseqfunc2 = ComposableSeq(input_types="sequences", output_types="sequences")
        comb = aseqfunc1 + aseqfunc2
        expect = "ComposableSeq(type='sequences') + " "ComposableSeq(type='sequences')"
        got = str(comb)
        self.assertEqual(got, expect)

    def test_composables_once(self):
        """composables can only be used in a single composition"""
        aseqfunc1 = ComposableSeq(input_types="sequences", output_types="sequences")
        aseqfunc2 = ComposableSeq(input_types="sequences", output_types="sequences")
        aseqfunc1 + aseqfunc2
        with self.assertRaises(AssertionError):
            aseqfunc3 = ComposableSeq(input_types="sequences", output_types="sequences")
            aseqfunc1 + aseqfunc3
        # the other order
        with self.assertRaises(AssertionError):
            aseqfunc3 = ComposableSeq(input_types="sequences", output_types="sequences")
            aseqfunc3 + aseqfunc2

    def test_composable_to_self(self):
        """this should raise a ValueError"""
        app1 = ComposableSeq(input_types="sequences", output_types="sequences")
        with self.assertRaises(ValueError):
            _ = app1 + app1

    def test_disconnect(self):
        """disconnect breaks all connections and allows parts to be reused"""
        aseqfunc1 = ComposableSeq(input_types="sequences", output_types="sequences")
        aseqfunc2 = ComposableSeq(input_types="sequences", output_types="sequences")
        aseqfunc3 = ComposableSeq(input_types="sequences", output_types="sequences")
        comb = aseqfunc1 + aseqfunc2 + aseqfunc3
        comb.disconnect()
        self.assertEqual(aseqfunc1.input, None)
        self.assertEqual(aseqfunc1.output, None)
        self.assertEqual(aseqfunc3.input, None)
        self.assertEqual(aseqfunc3.output, None)
        # should be able to compose a new one now
        aseqfunc1 + aseqfunc3

    def test_apply_to(self):
        """correctly applies iteratively"""
        from cogent3.core.alignment import SequenceCollection

        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        reader = io_app.load_unaligned(format="fasta", moltype="dna")
        got = reader.apply_to(dstore, show_progress=False)
        self.assertEqual(len(got), len(dstore))
        # should also be able to apply the results to another composable func
        min_length = sample_app.min_length(10)
        got = min_length.apply_to(got, show_progress=False)
        self.assertEqual(len(got), len(dstore))
        # should work on a chained function
        proc = reader + min_length
        got = proc.apply_to(dstore, show_progress=False)
        self.assertEqual(len(got), len(dstore))
        # and works on a list of just strings
        got = proc.apply_to([str(m) for m in dstore], show_progress=False)
        self.assertEqual(len(got), len(dstore))
        # or a single string
        got = proc.apply_to(str(dstore[0]), show_progress=False)
        self.assertEqual(len(got), 1)
        self.assertIsInstance(got[0], SequenceCollection)
        # raises ValueError if empty list
        with self.assertRaises(ValueError):
            proc.apply_to([])

        # raises ValueError if list with empty string
        with self.assertRaises(ValueError):
            proc.apply_to(["", ""])

    def test_apply_to_strings(self):
        """apply_to handles strings as paths"""
        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        dstore = [str(m) for m in dstore]
        with TemporaryDirectory(dir=".") as dirname:
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            min_length = sample_app.min_length(10)
            outpath = os.path.join(os.getcwd(), dirname, "delme.tinydb")
            writer = io_app.write_db(outpath)
            process = reader + min_length + writer
            # create paths as strings
            r = process.apply_to(dstore, show_progress=False)
            self.assertEqual(len(process.data_store.logs), 1)
            process.data_store.close()

    def test_apply_to_non_unique_identifiers(self):
        """should fail if non-unique names"""
        dstore = [
            "brca1.bats.fasta",
            "brca1.apes.fasta",
        ]
        with TemporaryDirectory(dir=".") as dirname:
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            min_length = sample_app.min_length(10)
            process = reader + min_length
            with self.assertRaises(ValueError):
                process.apply_to(dstore)

    def test_apply_to_logging(self):
        """correctly creates log file"""
        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        with TemporaryDirectory(dir=".") as dirname:
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            min_length = sample_app.min_length(10)
            outpath = os.path.join(os.getcwd(), dirname, "delme.tinydb")
            writer = io_app.write_db(outpath)
            process = reader + min_length + writer
            r = process.apply_to(dstore, show_progress=False)
            # always creates a log
            self.assertEqual(len(process.data_store.logs), 1)
            process.data_store.close()

    def test_apply_to_logger(self):
        """correctly uses user provided logger"""
        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        with TemporaryDirectory(dir=".") as dirname:
            LOGGER = CachingLogger()
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            min_length = sample_app.min_length(10)
            outpath = os.path.join(os.getcwd(), dirname, "delme.tinydb")
            writer = io_app.write_db(outpath)
            process = reader + min_length + writer
            r = process.apply_to(dstore, show_progress=False, logger=LOGGER)
            self.assertEqual(len(process.data_store.logs), 1)
            process.data_store.close()

    def test_apply_to_invalid_logger(self):
        """incorrect logger value raises TypeError"""
        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        for logger_val in (True, "somepath.log"):
            with TemporaryDirectory(dir=".") as dirname:
                reader = io_app.load_aligned(format="fasta", moltype="dna")
                min_length = sample_app.min_length(10)
                outpath = os.path.join(os.getcwd(), dirname, "delme.tinydb")
                writer = io_app.write_db(outpath)
                process = reader + min_length + writer
                with self.assertRaises(TypeError):
                    process.apply_to(dstore, show_progress=False, logger=logger_val)
                process.data_store.close()

    def test_apply_to_not_completed(self):
        """correctly creates notcompleted"""
        dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
        with TemporaryDirectory(dir=".") as dirname:
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            # trigger creation of notcompleted
            min_length = sample_app.min_length(3000)
            outpath = os.path.join(os.getcwd(), dirname, "delme.tinydb")
            writer = io_app.write_db(outpath)
            process = reader + min_length + writer
            r = process.apply_to(dstore, show_progress=False)
            self.assertEqual(len(process.data_store.incomplete), 3)
            process.data_store.close()

    def test_apply_to_not_partially_done(self):
        """correctly applies process when result already partially done"""
        dstore = io_app.get_data_store("data", suffix="fasta")
        num_records = len(dstore)
        with TemporaryDirectory(dir=".") as dirname:
            dirname = pathlib.Path(dirname)
            reader = io_app.load_aligned(format="fasta", moltype="dna")
            outpath = dirname / "delme.tinydb"
            writer = io_app.write_db(outpath)
            _ = writer(reader(dstore[0]))
            writer.data_store.close()

            writer = io_app.write_db(outpath, if_exists="ignore")
            process = reader + writer
            _ = process.apply_to(dstore, show_progress=False)
            writer.data_store.close()
            dstore = io_app.get_data_store(outpath)
            self.assertEqual(len(dstore), num_records)
            dstore.close()


class TestNotCompletedResult(TestCase):
    def test_err_result(self):
        """excercise creation of NotCompletedResult"""
        result = NotCompleted("SKIP", "this", "some obj")
        self.assertFalse(result)
        self.assertEqual(result.origin, "this")
        self.assertEqual(result.message, "some obj")
        self.assertIs(result.source, None)

        # check source correctly deduced from provided object
        fake_source = Mock()
        fake_source.source = "blah"
        del fake_source.info
        result = NotCompleted("SKIP", "this", "err", source=fake_source)
        self.assertIs(result.source, "blah")

        fake_source = Mock()
        del fake_source.source
        fake_source.info.source = "blah"
        result = NotCompleted("SKIP", "this", "err", source=fake_source)
        self.assertIs(result.source, "blah")

        try:
            _ = 0
            raise ValueError("error message")
        except ValueError as err:
            result = NotCompleted("SKIP", "this", err.args[0])

        self.assertEqual(result.message, "error message")

    def test_str(self):
        """str representation correctly represents parameterisations"""
        func = select_translatable()
        got = str(func)
        self.assertEqual(
            got,
            "select_translatable(type='sequences', "
            "moltype='dna', gc=1, "
            "allow_rc=False,\ntrim_terminal_stop=True)",
        )

        func = select_translatable(allow_rc=True)
        got = str(func)
        self.assertEqual(
            got,
            "select_translatable(type='sequences', "
            "moltype='dna', gc=1, "
            "allow_rc=True,\ntrim_terminal_stop=True)",
        )

        nodegen = omit_degenerates()
        got = str(nodegen)
        self.assertEqual(
            got,
            "omit_degenerates(type='aligned', moltype=None, "
            "gap_is_degen=True,\nmotif_length=1)",
        )
        ml = min_length(100)
        got = str(ml)
        self.assertEqual(
            got,
            "min_length(type='sequences', length=100, "
            "motif_length=1, subtract_degen=True,\n"
            "moltype=None)",
        )

        qt = quick_tree()
        self.assertEqual(str(qt), "quick_tree(type='tree', drop_invalid=False)")


class TestPicklable(TestCase):
    def test_composite_pickleable(self):
        """composable functions should be pickleable"""

        from cogent3.app import align, evo, io, sample, translate, tree

        read = io.load_aligned(moltype="dna")
        dumps(read)
        trans = translate.select_translatable()
        dumps(trans)
        aln = align.progressive_align("nucleotide")
        dumps(aln)
        just_nucs = sample.omit_degenerates(moltype="dna")
        dumps(just_nucs)
        limit = sample.fixed_length(1000, random=True)
        dumps(limit)
        mod = evo.model("HKY85")
        dumps(mod)
        qt = tree.quick_tree()
        dumps(qt)
        proc = read + trans + aln + just_nucs + limit + mod
        dumps(proc)

    def test_not_completed_result(self):
        """should survive roundtripping pickle"""
        err = NotCompleted("FAIL", "mytest", "can we roundtrip")
        p = dumps(err)
        new = loads(p)
        self.assertEqual(err.type, new.type)
        self.assertEqual(err.message, new.message)
        self.assertEqual(err.source, new.source)
        self.assertEqual(err.origin, new.origin)

    def test_triggers_bugcatcher(self):
        """a composable that does not trap failures returns NotCompletedResult
        requesting bug report"""
        from cogent3.app import io

        read = io.load_aligned(moltype="dna")
        read.func = lambda x: None
        got = read("somepath.fasta")
        self.assertIsInstance(got, NotCompleted)
        self.assertEqual(got.type, "BUG")


def _demo(ctx, expect):
    return ctx.frame_start == expect


# for testing appify
@appify(SERIALISABLE_TYPE, SERIALISABLE_TYPE)
def slicer(val, index=2):
    """my docstring"""
    return val[:index]


class TestUserFunction(TestCase):
    def foo(self, val, *args, **kwargs):
        return val[:4]

    def bar(self, val, *args, **kwargs):
        return val.distance_matrix(calc="hamming", show_progress=False)

    def test_user_function_custom_variables(self):
        # not sure what this is meant to be testing
        demo = user_function(
            _demo, ("aligned", "serialisable"), ("aligned", "serialisable")
        )
        frame_start = 2
        demo.frame_start = frame_start
        self.assertTrue(demo(demo, 2))

    def test_user_function(self):
        """composable functions should be user definable"""
        from cogent3 import make_aligned_seqs

        u_function = user_function(self.foo, "aligned", "aligned")

        aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
        got = u_function(aln)

        self.assertEqual(got.to_dict(), {"a": "GCAA", "b": "GCTT"})

    def test_user_function_multiple(self):
        """user defined composable functions should not interfere with each other"""
        from cogent3 import make_aligned_seqs
        from cogent3.core.alignment import Alignment

        u_function_1 = user_function(self.foo, "aligned", "aligned")
        u_function_2 = user_function(self.bar, "aligned", "pairwise_distances")

        aln_1 = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
        data = dict([("s1", "ACGTACGTA"), ("s2", "GTGTACGTA")])
        aln_2 = Alignment(data=data, moltype="dna")

        got_1 = u_function_1(aln_1)
        got_2 = u_function_2(aln_2)
        self.assertEqual(got_1.to_dict(), {"a": "GCAA", "b": "GCTT"})
        self.assertEqual(got_2, {("s1", "s2"): 2.0, ("s2", "s1"): 2.0})

    def test_appify(self):
        """acts like a decorator should!"""
        self.assertEqual(slicer.__doc__, "my docstring")
        self.assertEqual(slicer.__name__, "slicer")
        app = slicer()
        self.assertTrue(SERIALISABLE_TYPE in app._input_types)
        self.assertTrue(SERIALISABLE_TYPE in app._output_types)
        self.assertEqual(app(list(range(4))), [0, 1])
        app2 = slicer(index=3)
        self.assertEqual(app2(list(range(4))), [0, 1, 2])

    def test_appify_pickle(self):
        """appified function should be pickleable"""
        app = slicer(index=6)
        dumped = dumps(app)
        loaded = loads(dumped)
        self.assertEqual(loaded(list(range(10))), list(range(6)))

    def test_user_function_repr(self):
        u_function_1 = user_function(self.foo, "aligned", "aligned")
        u_function_2 = user_function(self.bar, "aligned", "pairwise_distances")
        self.assertEqual(
            repr(u_function_1), "user_function(name='foo', module='test_composable')"
        )
        self.assertEqual(
            repr(u_function_2), "user_function(name='bar', module='test_composable')"
        )

    def test_user_function_str(self):
        u_function_1 = user_function(self.foo, "aligned", "aligned")
        u_function_2 = user_function(self.bar, "aligned", "pairwise_distances")
        self.assertEqual(
            str(u_function_1), "user_function(name='foo', module='test_composable')"
        )
        self.assertEqual(
            str(u_function_2), "user_function(name='bar', module='test_composable')"
        )
        # added into a composable func
        loader = io_app.load_aligned()
        proc = loader + u_function_1
        got = str(proc)
        self.assertTrue(got.startswith("load_aligned"))

    def test_user_function_with_args_kwargs(self):
        """correctly handles definition with args, kwargs"""
        from math import log

        def product(val, multiplier, take_log=False):
            result = val * multiplier
            if take_log:
                result = log(result)

            return result

        # without defining any args, kwargs
        ufunc = user_function(
            product,
            SERIALISABLE_TYPE,
            SERIALISABLE_TYPE,
        )
        self.assertEqual(ufunc(2, 2), 4)
        self.assertEqual(ufunc(2, 2, take_log=True), log(4))

        # defining default arg2
        ufunc = user_function(
            product,
            SERIALISABLE_TYPE,
            SERIALISABLE_TYPE,
            2,
        )
        self.assertEqual(ufunc(2), 4)
        self.assertEqual(ufunc(2, take_log=True), log(4))

        # defining default kwarg only
        ufunc = user_function(
            product, SERIALISABLE_TYPE, SERIALISABLE_TYPE, take_log=True
        )
        self.assertEqual(ufunc(2, 2), log(4))
        self.assertEqual(ufunc(2, 2, take_log=False), 4)

        # defining default arg and kwarg
        ufunc = user_function(
            product, SERIALISABLE_TYPE, SERIALISABLE_TYPE, 2, take_log=True
        )
        self.assertEqual(ufunc(2), log(4))


if __name__ == "__main__":
    main()

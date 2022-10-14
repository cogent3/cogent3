import inspect
import os
import pathlib
import pickle

from pickle import dumps, loads
from tempfile import TemporaryDirectory
from typing import Set, Tuple
from unittest import main
from unittest.mock import Mock

import pytest

from numpy import array, ndarray
from scitrack import CachingLogger

from cogent3 import make_aligned_seqs
from cogent3.app import align, evo
from cogent3.app import io as io_app
from cogent3.app import sample as sample_app
from cogent3.app import translate, tree
from cogent3.app.composable import (
    GENERIC,
    NON_COMPOSABLE,
    WRITER,
    NotCompleted,
    __app_registry,
    _add,
    _get_raw_hints,
    appify,
    define_app,
    get_object_provenance,
    is_composable,
    user_function,
)
from cogent3.app.data_store_new import DataStoreDirectory, get_data_source
from cogent3.app.sample import min_length, omit_degenerates
from cogent3.app.translate import select_translatable
from cogent3.app.tree import quick_tree
from cogent3.app.typing import (
    SERIALISABLE_TYPE,
    AlignedSeqsType,
    PairwiseDistanceType,
)
from cogent3.core.alignment import Alignment, SequenceCollection


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


@pytest.fixture(scope="function")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("datastore")


def test_composable():
    """correctly form string"""

    @define_app
    class app_dummyclass_1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    @define_app
    class app_dummyclass_2:
        def __init__(self, b):
            self.b = b

        def main(self, val: int) -> int:
            return val

    aseqfunc1 = app_dummyclass_1(1)
    aseqfunc2 = app_dummyclass_2(2)
    comb = aseqfunc1 + aseqfunc2
    expect = "app_dummyclass_1(a=1) + " "app_dummyclass_2(b=2)"
    got = str(comb)
    assert got == expect
    __app_registry.pop(get_object_provenance(app_dummyclass_1), None)
    __app_registry.pop(get_object_provenance(app_dummyclass_2), None)


def test_composables_once():
    """composables can only be used in a single composition"""

    @define_app
    class app_dummyclass_1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    @define_app
    class app_dummyclass_2:
        def __init__(self, b):
            self.b = b

        def main(self, val: int) -> int:
            return val

    @define_app
    class app_dummyclass_3:
        def __init__(self, c):
            self.c = c

        def main(self, val: int) -> int:
            return val

    one = app_dummyclass_1(1)
    two = app_dummyclass_2(2)
    three = app_dummyclass_3(3)
    one + three
    with pytest.raises(ValueError):
        two + three  # three already has an input

    __app_registry.pop(get_object_provenance(app_dummyclass_1), None)
    __app_registry.pop(get_object_provenance(app_dummyclass_2), None)
    __app_registry.pop(get_object_provenance(app_dummyclass_3), None)


def test_composable_to_self():
    """this should raise a ValueError"""

    @define_app
    class app_dummyclass_1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    app1 = app_dummyclass_1(1)
    with pytest.raises(ValueError):
        _ = app1 + app1

    __app_registry.pop(get_object_provenance(app_dummyclass_1), None)


def test_disconnect():
    """disconnect breaks all connections and allows parts to be reused"""

    @define_app
    class app_dummyclass_1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    @define_app
    class app_dummyclass_2:
        def __init__(self, b):
            self.b = b

        def main(self, val: int) -> int:
            return val

    @define_app
    class app_dummyclass_3:
        def __init__(self, c):
            self.c = c

        def main(self, val: int) -> int:
            return val

    aseqfunc1 = app_dummyclass_1(1)
    aseqfunc2 = app_dummyclass_2(2)
    aseqfunc3 = app_dummyclass_3(3)
    comb = aseqfunc1 + aseqfunc2 + aseqfunc3
    comb.disconnect()
    assert aseqfunc1.input is None
    assert aseqfunc3.input is None
    # should be able to compose a new one now
    aseqfunc1 + aseqfunc3

    __app_registry.pop(get_object_provenance(app_dummyclass_1), None)
    __app_registry.pop(get_object_provenance(app_dummyclass_2), None)
    __app_registry.pop(get_object_provenance(app_dummyclass_3), None)


def test_as_completed():
    """correctly applies iteratively"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    reader = io_app.load_unaligned(format="fasta", moltype="dna")
    got = list(reader.as_completed(dstore))
    assert len(got) == len(dstore)
    # should also be able to apply the results to another composable func
    min_length = sample_app.min_length(10)
    got = list(min_length.as_completed(got))
    assert len(got) == len(dstore)
    # should work on a chained function
    proc = reader + min_length
    got = list(proc.as_completed(dstore))
    assert len(got) == len(dstore)
    # and works on a list of just strings
    got = list(proc.as_completed([str(m) for m in dstore]))
    assert len(got) == len(dstore)
    # or a single string
    got = list(proc.as_completed(str(dstore[0])))
    assert len(got) == 1
    assert isinstance(got[0].obj, SequenceCollection)
    # raises ValueError if empty list
    with pytest.raises(ValueError):
        proc.as_completed([])

    # raises ValueError if list with empty string
    with pytest.raises(ValueError):
        list(proc.as_completed(["", ""]))


@pytest.mark.parametrize("klass", (DataStoreDirectory,))
def test_apply_to_strings(tmp_dir, klass):
    """apply_to handles strings as paths"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    dstore = [str(m) for m in dstore]
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    outpath = tmp_dir / "test_apply_to_strings"
    writer = io_app.write_seqs_new(
        klass(outpath, if_dest_exists="overwrite", suffix="fasta")
    )
    process = reader + min_length + writer
    # create paths as strings
    _ = process.apply_to(dstore, id_from_source=get_data_source)
    assert len(process.data_store.logs) == 1


def test_apply_to_non_unique_identifiers(tmp_dir):
    """should fail if non-unique names"""
    dstore = [
        "brca1.bats.fasta",
        "brca1.apes.fasta",
    ]
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    outpath = tmp_dir / "test_apply_to_non_unique_identifiers"
    writer = io_app.write_seqs_new(
        DataStoreDirectory(outpath, if_dest_exists="overwrite", suffix="fasta")
    )
    process = reader + min_length + writer
    with pytest.raises(ValueError):
        process.apply_to(dstore)


def test_apply_to_logging(tmp_dir):
    """correctly creates log file"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    outpath = os.path.join(os.getcwd(), tmp_dir, "delme.tinydb")
    writer = io_app.write_db(outpath)
    process = reader + min_length + writer
    r = process.apply_to(dstore, show_progress=False)
    # always creates a log
    assert len(process.data_store.logs) == 1
    process.data_store.close()


def test_apply_to_logger(tmp_dir):
    """correctly uses user provided logger"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    LOGGER = CachingLogger()
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    outpath = os.path.join(os.getcwd(), tmp_dir, "delme.tinydb")
    writer = io_app.write_db(outpath)
    process = reader + min_length + writer
    r = process.apply_to(dstore, show_progress=False, logger=LOGGER)
    assert len(process.data_store.logs) == 1
    process.data_store.close()


def test_apply_to_invalid_logger(tmp_dir):
    """incorrect logger value raises TypeError"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    for logger_val in (True, "somepath.log"):
        reader = io_app.load_aligned(format="fasta", moltype="dna")
        min_length = sample_app.min_length(10)
        outpath = os.path.join(os.getcwd(), tmp_dir, "delme.tinydb")
        writer = io_app.write_db(outpath)
        process = reader + min_length + writer
        with pytest.raises(TypeError):
            process.apply_to(dstore, show_progress=False, logger=logger_val)
        process.data_store.close()


def test_apply_to_not_completed(tmp_dir):
    """correctly creates notcompleted"""
    dstore = io_app.get_data_store("data", suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    # trigger creation of notcompleted
    min_length = sample_app.min_length(3000)
    outpath = os.path.join(os.getcwd(), tmp_dir, "delme.tinydb")
    writer = io_app.write_db(outpath)
    process = reader + min_length + writer
    r = process.apply_to(dstore, show_progress=False)
    assert len(process.data_store.incomplete) == 3
    process.data_store.close()


def test_apply_to_not_partially_done(tmp_dir):
    """correctly applies process when result already partially done"""
    dstore = io_app.get_data_store("data", suffix="fasta")
    num_records = len(dstore)
    dirname = pathlib.Path(tmp_dir)
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
    assert len(dstore) == num_records
    dstore.close()


def test_err_result():
    """excercise creation of NotCompletedResult"""
    result = NotCompleted("SKIP", "this", "some obj")
    assert not result
    assert result.origin == "this"
    assert result.message == "some obj"
    assert result.source is None

    # check source correctly deduced from provided object
    fake_source = Mock()
    fake_source.source = "blah"
    del fake_source.info
    result = NotCompleted("SKIP", "this", "err", source=fake_source)
    assert result.source == "blah"

    try:
        _ = 0
        raise ValueError("error message")
    except ValueError as err:
        result = NotCompleted("SKIP", "this", err.args[0])

    assert result.message == "error message"


def test_str():
    """str representation correctly represents parameterisations"""
    func = select_translatable()
    got = str(func)
    assert (
        got
        == "select_translatable(moltype='dna', gc=1, allow_rc=False,\ntrim_terminal_stop=True)"
    )

    func = select_translatable(allow_rc=True)
    got = str(func)
    assert (
        got
        == "select_translatable(moltype='dna', gc=1, allow_rc=True, trim_terminal_stop=True)"
    )

    nodegen = omit_degenerates()
    got = str(nodegen)
    assert got == "omit_degenerates(moltype=None, gap_is_degen=True, motif_length=1)"
    ml = min_length(100)
    got = str(ml)
    assert (
        got
        == "min_length(length=100, motif_length=1, subtract_degen=True, moltype=None)"
    )

    qt = quick_tree()
    assert str(qt) == "quick_tree(drop_invalid=False)"


def test_composite_pickleable():
    """composable functions should be pickleable"""

    read = io_app.load_aligned(moltype="dna")
    dumps(read)
    trans = translate.select_translatable()
    dumps(trans)
    aln = align.progressive_align("nucleotide")
    dumps(aln)
    just_nucs = sample_app.omit_degenerates(moltype="dna")
    dumps(just_nucs)
    limit = sample_app.fixed_length(1000, random=True)
    dumps(limit)
    mod = evo.model("HKY85")
    dumps(mod)
    qt = tree.quick_tree()
    dumps(qt)
    proc = read + trans + aln + just_nucs + limit + mod
    dumps(proc)


def test_not_completed_result():
    """should survive roundtripping pickle"""
    err = NotCompleted("FAIL", "mytest", "can we roundtrip")
    p = dumps(err)
    new = loads(p)
    assert err.type == new.type
    assert err.message == new.message
    assert err.source == new.source
    assert err.origin == new.origin


def test_triggers_bugcatcher():
    """a composable that does not trap failures returns NotCompletedResult
    requesting bug report"""

    read = io_app.load_aligned(moltype="dna")
    read.main = lambda x: None
    got = read("somepath.fasta")
    assert isinstance(got, NotCompleted)
    assert got.type == "BUG"


def _demo(ctx, expect):
    return ctx.frame_start == expect


# for testing appify
@appify(SERIALISABLE_TYPE, SERIALISABLE_TYPE)
def slicer(val, index=2):
    """my docstring"""
    return val[:index]


@define_app
def foo(val: AlignedSeqsType, *args, **kwargs) -> AlignedSeqsType:
    return val[:4]

@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def bar(val: AlignedSeqsType, num=3) -> PairwiseDistanceType:
    return val.distance_matrix(calc="hamming", show_progress=False)


for _app_ in (foo, bar, foo_without_arg_kwargs):
    __app_registry.pop(get_object_provenance(_app_), None)


def test_user_function():
    """composable functions should be user definable"""

    u_function = foo()

    aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
    got = u_function(aln)

    assert got.to_dict() == {"a": "GCAA", "b": "GCTT"}


def test_user_function_without_arg_kwargs():
    """composable functions should be user definable"""

    u_function = foo_without_arg_kwargs()

    aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
    got = u_function(aln)

    assert got.to_dict() == {"a": "GCAA", "b": "GCTT"}


def test_user_function_multiple():
    """user defined composable functions should not interfere with each other"""
    u_function_1 = foo()
    u_function_2 = bar()

    aln_1 = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
    data = dict([("s1", "ACGTACGTA"), ("s2", "GTGTACGTA")])
    aln_2 = Alignment(data=data, moltype="dna")

    got_1 = u_function_1(aln_1)
    got_2 = u_function_2(aln_2)
    assert got_1.to_dict() == {"a": "GCAA", "b": "GCTT"}
    assert got_2 == {("s1", "s2"): 2.0, ("s2", "s1"): 2.0}


def test_appify():
    """acts like a decorator should!"""
    assert slicer.__doc__ == "my docstring"
    assert slicer.__name__ == "slicer"
    app = slicer()
    assert SERIALISABLE_TYPE in app._input_types
    assert SERIALISABLE_TYPE in app._output_types
    assert app(list(range(4))) == [0, 1]
    app2 = slicer(index=3)
    assert app2(list(range(4))) == [0, 1, 2]


def test_appify_pickle():
    """appified function should be pickleable"""
    app = slicer(index=6)
    dumped = dumps(app)
    loaded = loads(dumped)
    assert loaded(list(range(10))) == list(range(6))


def test_user_function_repr():
    got = repr(bar(num=3))
    assert got == "bar(num=3)"


def test_user_function_str():
    got = str(bar(num=3))
    assert got == "bar(num=3)"


def test_user_function_with_args_kwargs():
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
    assert ufunc(2, 2) == 4
    assert ufunc(2, 2, take_log=True) == log(4)

    # defining default arg2
    ufunc = user_function(
        product,
        SERIALISABLE_TYPE,
        SERIALISABLE_TYPE,
        2,
    )
    assert ufunc(2) == 4
    assert ufunc(2, take_log=True) == log(4)

    # defining default kwarg only
    ufunc = user_function(product, SERIALISABLE_TYPE, SERIALISABLE_TYPE, take_log=True)
    assert ufunc(2, 2) == log(4)
    assert ufunc(2, 2, take_log=False) == 4

    # defining default arg and kwarg
    ufunc = user_function(
        product, SERIALISABLE_TYPE, SERIALISABLE_TYPE, 2, take_log=True
    )
    assert ufunc(2) == log(4)


def test_app_registry():
    """correctly registers apps"""

    @define_app
    class app_test_registry1:
        def main(self, data: int) -> int:
            return data

    assert __app_registry["test_composable.app_test_registry1"]

    # delete it to not include in app available apps
    __app_registry.pop(get_object_provenance(app_test_registry1), None)


def test_app_is_composable():
    """check is_composable for composable apps"""

    @define_app
    class app_test_iscomposable1:
        def main(self, data: int) -> int:
            return data

    assert is_composable(app_test_iscomposable1)

    # delete it to not include in app available apps
    __app_registry.pop(get_object_provenance(app_test_iscomposable1), None)


def test_app_is_not_composable():
    """check is_composable for non-composable apps"""

    class app_not_composable1:
        def main(self, data: int) -> int:
            return data

    assert not is_composable(app_not_composable1)


def test_concat_not_composable():
    from cogent3.app.sample import concat

    assert not is_composable(concat)


def test_composed_func_pickleable():

    ml = min_length(100)
    no_degen = omit_degenerates(moltype="dna")
    app = ml + no_degen

    unpickled = pickle.loads(pickle.dumps((app)))
    assert unpickled.input is not None


def test_composable_variable_positional_args():
    """correctly associate argument vals with their names when have variable
    positional args"""

    @define_app
    class pos_var_pos1:
        def __init__(self, a, b, *args):
            self.a = a
            self.b = b
            self.args = args

        def main(self, val: int) -> int:
            return val

    instance = pos_var_pos1(2, 3, 4, 5, 6)
    assert instance._init_vals == {"a": 2, "b": 3, "args": (4, 5, 6)}

    __app_registry.pop(get_object_provenance(pos_var_pos1), None)


def test_composable_minimum_parameters():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    def test_func1(arg1) -> int:
        return 1

    with pytest.raises(ValueError):
        _, _ = _get_raw_hints(test_func1, 2)


def test_composable_return_type_hint():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    def test_func1(arg1):
        return 1

    with pytest.raises(TypeError):
        _, _ = _get_raw_hints(test_func1, 1)


def test_composable_firstparam_type_hint():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    def test_func1(arg1) -> int:
        return 1

    with pytest.raises(TypeError):
        _, _ = _get_raw_hints(test_func1, 1)


def test_composable_firstparam_type_is_None():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    def test_func1(arg1: None) -> int:
        return 1

    with pytest.raises(TypeError):
        _, _ = _get_raw_hints(test_func1, 1)


def test_composable_return_type_is_None():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    def test_func1(arg1: int) -> None:
        return

    with pytest.raises(TypeError):
        _, _ = _get_raw_hints(test_func1, 1)


def test_composable_variable_positional_args_and_kwargs():
    """correctly associate argument vals with their names when have variable
    positional args and kwargs"""

    @define_app
    class pos_var_pos_kw2:
        def __init__(self, a, *args, c=False):
            self.a = a
            self.c = c
            self.args = args

        def main(self, val: int) -> int:
            return val

    instance = pos_var_pos_kw2(2, 3, 4, 5, 6, c=True)
    assert instance._init_vals == {"a": 2, "args": (3, 4, 5, 6), "c": True}

    __app_registry.pop(get_object_provenance(pos_var_pos_kw2), None)


def test_app_decoration_fails_with_slots():
    with pytest.raises(NotImplementedError):

        @define_app
        class app_not_supported_slots1:
            __slots__ = ("a",)

            def __init__(self, a):
                self.a = a

            def main(self, val: int) -> int:
                return val


def test_repeated_decoration():
    @define_app
    class app_decorated_repeated1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    with pytest.raises(TypeError):
        define_app(app_decorated_repeated1)

    __app_registry.pop(get_object_provenance(app_decorated_repeated1), None)


def test_recursive_decoration():
    @define_app
    class app_docorated_recursive1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            define_app(app_docorated_recursive1)
            return val

    with pytest.raises(TypeError):
        app_docorated_recursive1().main(1)

    __app_registry.pop(get_object_provenance(app_docorated_recursive1), None)


def test_inheritance_from_decorated_class():
    @define_app
    class app_decorated_first1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    with pytest.raises(TypeError):

        @define_app
        class app_inherits_decorated1(app_decorated_first1):
            def __init__(self, a):
                self.a = a

            def main(self, val: int) -> int:
                return val

    __app_registry.pop(get_object_provenance(app_decorated_first1), None)


# have to define this at module level for pickling to work
@define_app
def func2app(arg1: int, exponent: int) -> float:
    return arg1 ** exponent


def test_decorate_app_function():
    """works on functions now"""

    sqd = func2app(exponent=2)
    assert sqd(3) == 9
    assert inspect.isclass(func2app)

    __app_registry.pop(get_object_provenance(func2app), None)


def test_roundtrip_decorated_function():
    """decorated function can be pickled/unpickled"""

    sqd = func2app(exponent=2)
    u = pickle.loads(pickle.dumps(sqd))
    assert u(4) == 16

    __app_registry.pop(get_object_provenance(func2app), None)


def test_decorated_func_optional():
    @define_app(app_type=NON_COMPOSABLE)
    def power(val: int, pow: int = 1) -> int:
        return val ** pow

    sqd = power(2)
    assert sqd(3) == 9

    __app_registry.pop(get_object_provenance(power), None)


def test_decorated_func_repr():
    def kw(val: int = 1) -> int:
        return val ** val

    def kw_kw(val: int = 1, pow: int = 1) -> int:
        return val ** pow

    def pos(val: int) -> int:
        return val ** val

    def pos_pos(val: int, pow: int) -> int:
        return val ** pow

    def pos_kw(val: int, pow: int = 1) -> int:
        return val ** pow

    fns = {fn: func for fn, func in locals().items() if callable(func)}
    args = {"pos": 4, "kw": dict(pow=3)}
    for name, func in fns.items():
        app = define_app(func)
        if len(name.split("_")) == 1:
            instance = app()
            expect = f"{name}()"
        elif name.endswith("kw"):
            instance = app(**args["kw"])
            expect = f"{name}(pow={args['kw']['pow']})"
        else:
            instance = app(args["pos"])
            expect = f"{name}(pow={args['pos']})"

        assert repr(instance) == expect, name

        __app_registry.pop(get_object_provenance(instance), None)


def test_decorated_func_just_args():
    @define_app(app_type=NON_COMPOSABLE)
    def power(val: int, pow: int) -> int:
        return val ** pow

    sqd = power()
    assert sqd(3, 3) == 27

    __app_registry.pop(get_object_provenance(power), None)


@pytest.mark.parametrize(
    "meth",
    [
        "__call__",
        "__repr__",
        "__str__",
        "__new__",
        "__add__",
        "disconnect",
        "input",
        "apply_to",
        "_validate_data_type",
    ],
)
def test_forbidden_methods_composable_app(meth):
    class app_forbidden_methods1:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    def function1():
        pass

    setattr(app_forbidden_methods1, meth, function1)
    with pytest.raises(TypeError):
        define_app(app_type=WRITER)(app_forbidden_methods1)


@pytest.mark.parametrize(
    "meth", ["__call__", "__repr__", "__str__", "__new__", "_validate_data_type"]
)
def test_forbidden_methods_non_composable_app(meth):
    class app_forbidden_methods2:
        def __init__(self, a):
            self.a = a

        def main(self, val: int) -> int:
            return val

    def function1():
        pass

    setattr(app_forbidden_methods2, meth, function1)
    with pytest.raises(TypeError):
        define_app(app_type=NON_COMPOSABLE)(app_forbidden_methods2)


def test_add_non_composable_apps():
    @define_app(app_type=NON_COMPOSABLE)
    class app_non_composable1:
        def __init__(self):
            pass

        def main(self, val: int) -> int:
            return val

    @define_app(app_type=NON_COMPOSABLE)
    class app_non_composable2:
        def __init__(self):
            pass

        def main(self, val: int) -> int:
            return val

    setattr(app_non_composable1, "__add__", _add)
    setattr(app_non_composable2, "__add__", _add)
    app1 = app_non_composable1()
    app2 = app_non_composable2()
    with pytest.raises(TypeError):
        app1 + app2

    __app_registry.pop(get_object_provenance(app_non_composable1), None)
    __app_registry.pop(get_object_provenance(app_non_composable2), None)


_types_null = (list, []), (tuple, ())


@pytest.mark.parametrize("in_type,input", _types_null)
def test_handles_null_series_input(in_type, input):
    """apps correctly handle null output"""

    @define_app
    def null_in(val: in_type, pow: int) -> int:
        return 2

    app = null_in(pow=2)
    got = app(input)
    assert isinstance(got, NotCompleted)

    __app_registry.pop(get_object_provenance(null_in), None)


@pytest.mark.parametrize("ret_type", (0, array([]), [], {}))
def test_handles_null_output(ret_type):
    """apps correctly handle null output"""

    @define_app
    def null_out(val: ndarray, pow: int) -> int:
        return ret_type

    app = null_out(pow=2)
    d = array([3, 3])
    got = app(d)
    assert isinstance(got, type(ret_type))

    __app_registry.pop(get_object_provenance(null_out), None)


def test_handles_None():
    """apps correctly handle null output"""

    @define_app
    def none_out(val: ndarray, pow: int) -> int:
        return None

    @define_app
    def take_int(val: int) -> int:
        return val

    app = none_out(pow=2)
    d = array([3, 3])
    got = app(d)
    assert isinstance(got, NotCompleted)

    app = none_out(pow=2) + take_int()
    d = array([3, 3])
    got = app(d)
    assert isinstance(got, NotCompleted)

    __app_registry.pop(get_object_provenance(none_out), None)


def test_validate_data_type_not_completed_pass_through():
    # returns the instance of a NotCompleted created by an input
    @define_app
    def take_int1(val: int) -> int:
        return NotCompleted("ERROR", "take_int1", "external to app", source="unknown")

    @define_app
    def take_int2(val: int) -> int:
        return val

    app = take_int1() + take_int2()
    got = app(2)
    assert got.origin == "take_int1"

    __app_registry.pop(get_object_provenance(take_int1), None)
    __app_registry.pop(get_object_provenance(take_int2), None)


@pytest.mark.parametrize("first,ret", ((Tuple[Set[str]], int), (int, Tuple[Set[str]])))
def test_complex_type(first, ret):
    # disallow >2-deep nesting of types for first arg and return type
    with pytest.raises(TypeError):

        @define_app
        class x:
            def main(self, data: first) -> ret:
                return data


@pytest.mark.parametrize("hint", (Tuple[Set[str]], Tuple[Tuple[Set[str]]]))
def test_complex_type_depths(hint):
    # disallow >2-deep nesting of types for first arg and return type
    with pytest.raises(TypeError):

        @define_app
        class x:
            def main(self, data: hint) -> bool:
                return True


@pytest.mark.parametrize("hint", (int, Set[str]))
def test_complex_type_allowed_depths(hint):
    # allowed <=2-deep nesting of types
    @define_app
    class x:
        def main(self, data: hint) -> int:
            return int

    __app_registry.pop(get_object_provenance(x), None)


if __name__ == "__main__":
    main()

import inspect
import pickle
import shutil
from pathlib import Path
from pickle import dumps, loads
from unittest.mock import Mock

import pytest
from numpy import array, ndarray
from scitrack import CachingLogger

from cogent3 import get_app, make_aligned_seqs, open_data_store
from cogent3.app import align, evo, translate, tree
from cogent3.app import io as io_app
from cogent3.app import sample as sample_app
from cogent3.app import typing as c3types
from cogent3.app.composable import (
    NON_COMPOSABLE,
    WRITER,
    NotCompleted,
    _add,
    _get_raw_hints,
    define_app,
    is_app,
    is_app_composable,
    source_proxy,
)
from cogent3.app.data_store import (
    APPEND,
    OVERWRITE,
    READONLY,
    DataStoreDirectory,
    get_data_source,
)
from cogent3.app.sample import min_length, omit_degenerates
from cogent3.app.sqlite_data_store import DataStoreSqlite
from cogent3.app.translate import select_translatable
from cogent3.app.tree import quick_tree
from cogent3.app.typing import AlignedSeqsType, PairwiseDistanceType
from cogent3.util.union_dict import UnionDict


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("datastore")


@pytest.fixture
def fasta_dir(DATA_DIR, tmp_dir):
    filenames = DATA_DIR.glob("*.fasta")
    fasta_dir = tmp_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for fn in filenames:
        dest = fasta_dir / fn.name
        dest.write_text(fn.read_text())
    return fasta_dir


@pytest.fixture
def write_dir1(tmp_dir):
    write_dir1 = tmp_dir / "write1"
    write_dir1.mkdir(parents=True, exist_ok=True)
    yield write_dir1
    shutil.rmtree(write_dir1)


@pytest.fixture
def write_dir2(tmp_dir):
    write_dir2 = tmp_dir / "write2"
    write_dir2.mkdir(parents=True, exist_ok=True)
    yield write_dir2
    shutil.rmtree(write_dir2)


@pytest.fixture
def ro_dstore(fasta_dir):
    return DataStoreDirectory(fasta_dir, suffix="fasta", mode=READONLY)


@pytest.fixture
def completed_objects(ro_dstore):
    return {f"{Path(m.unique_id).stem}": m.read() for m in ro_dstore}


@pytest.fixture
def nc_objects():
    return {
        f"id_{i}": NotCompleted("ERROR", "location", "message", source=f"id_{i}")
        for i in range(3)
    }


@pytest.fixture
def log_data(DATA_DIR):
    path = DATA_DIR / "scitrack.log"
    return path.read_text()


@pytest.fixture
def full_dstore(write_dir1, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(write_dir1, suffix="fasta", mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())

    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)

    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


@pytest.fixture
def nc_dstore(tmp_dir, nc_objects):
    dstore = DataStoreDirectory(tmp_dir, suffix="fasta", mode=OVERWRITE)
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())

    return dstore


@pytest.fixture
def half_dstore1(write_dir1, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(write_dir1, suffix="fasta", mode=OVERWRITE)
    i = 0
    for id, data in nc_objects.items():
        dstore.write_not_completed(unique_id=id, data=data.to_json())
        i += 1
        if i >= len(nc_objects.items()) / 2:
            break
    i = 0
    for id, data in completed_objects.items():
        dstore.write(unique_id=id, data=data)
        i += 1
        if i >= len(completed_objects.items()) / 2:
            break
    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


@pytest.fixture
def half_dstore2(write_dir2, nc_objects, completed_objects, log_data):
    dstore = DataStoreDirectory(write_dir2, suffix="fasta", mode=OVERWRITE)
    i = -1
    for id, data in nc_objects.items():
        i += 1
        if i < len(nc_objects.items()) / 2:
            continue
        dstore.write_not_completed(unique_id=id, data=data.to_json())
    i = -1
    for id, data in completed_objects.items():
        i += 1
        if i < len(completed_objects.items()) / 2:
            continue
        dstore.write(unique_id=id, data=data)
    dstore.write_log(unique_id="scitrack.log", data=log_data)
    return dstore


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
    expect = "app_dummyclass_1(a=1) + app_dummyclass_2(b=2)"
    got = str(comb)
    assert got == expect


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


def test_as_completed(DATA_DIR):
    """correctly applies iteratively"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    reader = get_app("load_unaligned", format="fasta", moltype="dna")
    got = list(reader.as_completed(dstore, show_progress=False))
    assert len(got) == len(dstore)
    # should also be able to apply the results to another composable func
    min_length = get_app("sample.min_length", 10)
    got = list(min_length.as_completed(got, show_progress=False))
    assert len(got) == len(dstore)
    # should work on a chained function
    proc = reader + min_length
    got = list(proc.as_completed(dstore, show_progress=False))
    assert len(got) == len(dstore)
    # and works on a list of just strings
    got = list(proc.as_completed([str(m) for m in dstore], show_progress=False))
    assert len(got) == len(dstore)
    # or a single string
    path = str(Path(dstore[0].data_store.source) / dstore[0].unique_id)
    got = list(proc.as_completed(path, show_progress=False))
    assert len(got) == 1
    assert got[0].__class__.__name__.endswith("SequenceCollection")


@pytest.fixture(params=[Path, None, str])
def source_type(DATA_DIR, request):
    fname = "brca1.fasta"
    dstore = open_data_store(DATA_DIR, suffix="fasta")
    if request.param is not None:
        return request.param(DATA_DIR / fname)
    return next((m for m in dstore if m.unique_id == fname), None)


def test_composable_unwraps_source_proxy_as_completed(source_type):
    app = get_app("load_unaligned", format="fasta", moltype="dna")
    result = next(iter(app.as_completed([source_type], show_progress=False)))
    got = result.source
    assert got.endswith("brca1.fasta")
    assert not isinstance(got, source_proxy)


def test_composable_unwraps_source_proxy_call(source_type):
    app = get_app("load_unaligned", format="fasta", moltype="dna")
    result = app(source_type)
    got = result.source
    assert got.endswith("brca1.fasta")
    assert not isinstance(got, source_proxy)


@pytest.mark.parametrize("data", [(), ("", "")])
def test_as_completed_empty_data(data):
    """correctly applies iteratively"""
    reader = get_app("load_unaligned", format="fasta", moltype="dna")
    min_length = get_app("sample.min_length", 10)
    proc = reader + min_length

    # returns empty input
    got = proc.as_completed(data)
    assert got == []


@pytest.mark.parametrize(
    "data",
    [
        {"a": 2},
        UnionDict(a=2, source="blah.txt"),
        make_aligned_seqs(
            {"a": "ACGT"},
            info={"source": "blah.txt"},
            moltype="dna",
        ),
    ],
)
def test_as_completed_w_wout_source(data):
    @define_app
    def pass_through(val: dict | UnionDict | AlignedSeqsType) -> dict:
        return val

    app = pass_through()  # pylint: disable=not-callable,no-value-for-parameter
    got = list(app.as_completed([data], show_progress=False))
    assert bool(got), got


@pytest.mark.parametrize("klass", [DataStoreDirectory, DataStoreSqlite])
@pytest.mark.parametrize("cast", [str, Path])
def test_apply_to_strings(DATA_DIR, tmp_dir, klass, cast):
    """apply_to handles non-DataMember"""
    dname = "test_apply_to_strings"
    outpath = tmp_dir / dname

    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    dstore = [cast(str(m)) for m in dstore]
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    if klass == DataStoreDirectory:
        writer = io_app.write_seqs(klass(outpath, mode=OVERWRITE, suffix="fasta"))
    else:
        writer = io_app.write_seqs(klass(outpath, mode=OVERWRITE))
    process = reader + min_length + writer
    # create paths as strings
    _ = process.apply_to(dstore, id_from_source=get_data_source, show_progress=False)
    assert len(process.data_store.logs) == 1


@pytest.mark.parametrize("klass", [DataStoreDirectory, DataStoreSqlite])
@pytest.mark.parametrize("cast", [str, Path])
def test_as_completed_strings(DATA_DIR, tmp_dir, klass, cast):
    """as_completed handles non-DataMember"""
    dname = "test_apply_to_strings"
    outpath = tmp_dir / dname

    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    dstore = [cast(str(m)) for m in dstore]
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    if klass == DataStoreDirectory:
        writer = io_app.write_seqs(klass(outpath, mode=OVERWRITE, suffix="fasta"))
    else:
        writer = io_app.write_seqs(klass(outpath, mode=OVERWRITE))
    orig_length = len(writer.data_store)
    process = reader + min_length + writer
    # create paths as strings
    got = list(process.as_completed(dstore, show_progress=False))
    assert len(got) > orig_length


def test_apply_to_non_unique_identifiers(tmp_dir):
    """should fail if non-unique names"""
    dstore = [
        "brca1.fasta",
        "brca1.fasta",
    ]
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    outpath = tmp_dir / "test_apply_to_non_unique_identifiers"
    writer = io_app.write_seqs(
        DataStoreDirectory(outpath, mode=OVERWRITE, suffix="fasta"),
    )
    process = reader + min_length + writer
    with pytest.raises(ValueError):
        process.apply_to(dstore)


def test_apply_to_logging(DATA_DIR, tmp_dir):
    """correctly creates log file"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    process = reader + min_length + writer
    process.apply_to(dstore, show_progress=False)
    # always creates a log
    assert len(process.data_store.logs) == 1


def test_apply_to_logger(DATA_DIR, tmp_dir):
    """correctly uses user provided logger"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    LOGGER = CachingLogger()
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    process = reader + min_length + writer
    process.apply_to(dstore, show_progress=False, logger=LOGGER)
    assert len(process.data_store.logs) == 1


def test_apply_to_no_logger(DATA_DIR, tmp_dir):
    """correctly uses user provided logger"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    process = reader + min_length + writer
    process.apply_to(dstore, show_progress=False, logger=False)
    assert len(process.data_store.logs) == 0
    assert process.logger is None


@pytest.mark.parametrize("logger_val", [True, "somepath.log"])
def test_apply_to_invalid_logger(DATA_DIR, tmp_dir, logger_val):
    """incorrect logger value raises TypeError"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    min_length = sample_app.min_length(10)
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    process = reader + min_length + writer
    with pytest.raises(TypeError):
        process.apply_to(dstore, show_progress=False, logger=logger_val)


def test_apply_to_input_only_not_completed(DATA_DIR, nc_dstore, tmp_dir):
    """correctly skips notcompleted"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    # trigger creation of notcompleted
    outpath = tmp_dir / "delme.sqlitedb"
    out_dstore = open_data_store(outpath, mode="w")
    writer = io_app.write_db(out_dstore)
    process = (
        io_app.load_aligned(format="fasta", moltype="dna")
        + sample_app.min_length(3000)
        + writer
    )
    process.apply_to(dstore, show_progress=False)
    assert len(out_dstore.not_completed) == len(nc_dstore)


def test_apply_to_makes_not_completed(DATA_DIR, tmp_dir):
    """correctly creates notcompleted"""
    dstore = open_data_store(DATA_DIR, suffix="fasta", limit=3)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    # trigger creation of notcompleted
    min_length = sample_app.min_length(3000)
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    process = reader + min_length + writer
    process.apply_to(dstore, show_progress=False)
    assert len(out_dstore.not_completed) == 3


def test_apply_to_not_partially_done(DATA_DIR, tmp_dir):
    """correctly applies process when result already partially done"""
    dstore = open_data_store(DATA_DIR, suffix="fasta")
    num_records = len(dstore)
    reader = io_app.load_aligned(format="fasta", moltype="dna")
    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="w")
    writer = io_app.write_db(out_dstore)
    # doing the first one
    # turning off warning as apps are callable
    _ = writer(reader(dstore[0]))  # pylint: disable=not-callable
    writer.data_store.close()

    out_dstore = open_data_store(tmp_dir / "delme.sqlitedb", mode="a")
    writer = io_app.write_db(out_dstore)
    process = reader + writer
    _ = process.apply_to(dstore, show_progress=False)
    assert len(out_dstore) == num_records


@pytest.mark.xfail(reason="passes except when run in full test suite")
@pytest.mark.parametrize("show", [True, False])
def test_as_completed_progress(full_dstore, capsys, show):
    loader = get_app("load_unaligned", format="fasta", moltype="dna")
    omit = get_app("omit_degenerates")
    app = loader + omit
    list(app.as_completed(full_dstore.completed, show_progress=show))
    result = capsys.readouterr().err.splitlines()
    if show:
        assert len(result) > 0
        assert "100%" in result[-1]
    else:
        assert len(result) == 0


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
        msg = "error message"
        raise ValueError(msg)
    except ValueError as err:
        result = NotCompleted("SKIP", "this", err.args[0])

    assert result.message == "error message"


def test_str():
    """str representation correctly represents parameterisations"""
    func = select_translatable()
    got = str(func)
    assert (
        got
        == "select_translatable(moltype='dna', gc=1, allow_rc=False,\ntrim_terminal_stop=True, frame=None)"
    )

    func = select_translatable(allow_rc=True)
    got = str(func)
    assert got.startswith(
        "select_translatable(moltype='dna', gc=1, allow_rc=True, trim_terminal_stop=True",
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
    got = read("somepath.fasta")  # pylint: disable=not-callable
    assert isinstance(got, NotCompleted)
    assert got.type == "BUG"


def _demo(ctx, expect):
    return ctx.frame_start == expect


@define_app
def foo(val: AlignedSeqsType, *args, **kwargs) -> AlignedSeqsType:
    return val[:4]


@define_app
def foo_without_arg_kwargs(val: AlignedSeqsType) -> AlignedSeqsType:
    return val[:4]


@define_app
def bar(val: AlignedSeqsType, num=3) -> PairwiseDistanceType:
    return val.distance_matrix(calc="hamming")


def test_user_function():
    """composable functions should be user definable"""

    u_function = foo()

    aln = make_aligned_seqs(
        [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
        moltype="dna",
    )
    got = u_function(aln)

    assert got.to_dict() == {"a": "GCAA", "b": "GCTT"}


def test_user_function_without_arg_kwargs():
    """composable functions should be user definable"""

    u_function = foo_without_arg_kwargs()

    aln = make_aligned_seqs(
        [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
        moltype="dna",
    )
    got = u_function(aln)

    assert got.to_dict() == {"a": "GCAA", "b": "GCTT"}


def test_user_function_multiple():
    """user defined composable functions should not interfere with each other"""
    u_function_1 = foo()
    u_function_2 = bar()

    aln_1 = make_aligned_seqs(
        [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
        moltype="dna",
    )
    data = {"s1": "ACGTACGTA", "s2": "GTGTACGTA"}
    aln_2 = make_aligned_seqs(data, moltype="dna")

    got_1 = u_function_1(aln_1)
    got_2 = u_function_2(aln_2)
    assert got_1.to_dict() == {"a": "GCAA", "b": "GCTT"}
    assert got_2 == {("s1", "s2"): 2.0, ("s2", "s1"): 2.0}


def test_user_function_repr():
    got = repr(bar(num=3))
    assert got == "bar(num=3)"


def test_user_function_str():
    got = str(bar(num=3))
    assert got == "bar(num=3)"


def test_decorated_app_is_app():
    """check is_app for define_app decorated apps"""

    @define_app
    class app_test_isapp1:
        def main(self, data: int) -> int:
            return data

    assert is_app(app_test_isapp1)


def test_undecorated_app_is_not_an_app():
    """check is_app for non-decorated apps"""

    class app_not_composable1:
        def main(self, data: int) -> int:
            return data

    assert not is_app(app_not_composable1)


def test_concat_not_composable():
    from cogent3.app.sample import concat

    assert not is_app_composable(concat)


def test_composed_func_pickleable():
    ml = min_length(100)
    no_degen = omit_degenerates(moltype="dna")
    app = ml + no_degen

    unpickled = pickle.loads(pickle.dumps(app))
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


# have to define this at module level for pickling to work
@define_app
def func2app(arg1: int, exponent: int) -> float:
    return arg1**exponent


def test_decorate_app_function():
    """works on functions now"""

    sqd = func2app(exponent=2)
    assert sqd(3) == 9
    assert inspect.isclass(func2app)


def test_roundtrip_decorated_function():
    """decorated function can be pickled/unpickled"""

    sqd = func2app(exponent=2)
    u = pickle.loads(pickle.dumps(sqd))
    assert u(4) == 16


def test_decorated_func_optional():
    @define_app(app_type=NON_COMPOSABLE)
    def power(val: int, pow: int = 1) -> int:
        return val**pow

    sqd = power(2)
    assert sqd(3) == 9


def test_decorated_func_repr():
    def kw(val: int = 1) -> int:
        return val**val

    def kw_kw(val: int = 1, pow: int = 1) -> int:
        return val**pow

    def pos(val: int) -> int:
        return val**val

    def pos_pos(val: int, pow: int) -> int:
        return val**pow

    def pos_kw(val: int, pow: int = 1) -> int:
        return val**pow

    fns = {fn: func for fn, func in locals().items() if callable(func)}
    args = {"pos": 4, "kw": {"pow": 3}}
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


def test_decorated_func_just_args():
    @define_app(app_type=NON_COMPOSABLE)
    def power(val: int, pow: int) -> int:
        return val**pow

    sqd = power()
    assert sqd(3, 3) == 27


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
    "meth",
    ["__call__", "__repr__", "__str__", "__new__", "_validate_data_type"],
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

    app_non_composable1.__add__ = _add
    app_non_composable2.__add__ = _add
    app1 = app_non_composable1()
    app2 = app_non_composable2()
    with pytest.raises(TypeError):
        app1 + app2


_types_null = (list, []), (tuple, ())


@pytest.mark.parametrize(("in_type", "input"), _types_null)
def test_handles_null_series_input(in_type, input):
    """apps correctly handle null output"""

    @define_app
    def null_in(val: in_type, pow: int) -> int:
        return 2

    app = null_in(pow=2)
    got = app(input)
    assert isinstance(got, NotCompleted)


@pytest.mark.parametrize("ret_type", [0, array([]), [], {}])
def test_handles_null_output(ret_type):
    """apps correctly handle null output"""

    @define_app
    def null_out(val: ndarray, pow: int) -> int:
        return ret_type

    app = null_out(pow=2)
    d = array([3, 3])
    got = app(d)
    assert isinstance(got, type(ret_type))


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


@pytest.mark.parametrize(
    ("first", "ret"),
    [(tuple[set[str]], int), (int, tuple[set[str]])],
)
def test_complex_type(first, ret):
    # disallow >2-deep nesting of types for first arg and return type
    with pytest.raises(TypeError):

        @define_app
        class x:
            def main(self, data: first) -> ret:
                return data


@pytest.mark.parametrize("hint", [tuple[set[str]], tuple[tuple[set[str]]]])
def test_complex_type_depths(hint):
    # disallow >2-deep nesting of types for first arg and return type
    with pytest.raises(TypeError):

        @define_app
        class x:
            def main(self, data: hint) -> bool:
                return True


@pytest.mark.parametrize("hint", [int, set[str]])
def test_complex_type_allowed_depths(hint):
    # allowed <=2-deep nesting of types
    @define_app
    class x:
        def main(self, data: hint) -> int:
            return int


def test_apply_to_only_appends(half_dstore1, half_dstore2):
    half_dstore1 = open_data_store(
        half_dstore1.source,
        suffix=half_dstore1.suffix,
        mode=APPEND,
    )
    reader1 = io_app.load_aligned(format="fasta", moltype="dna")
    min_length1 = sample_app.min_length(10)
    writer1 = io_app.write_seqs(half_dstore1)
    process1 = reader1 + min_length1 + writer1

    # create paths as strings
    dstore1 = [
        str(Path(m.data_store.source) / m.unique_id) for m in half_dstore1.completed
    ]
    # check does not modify dstore when applied to same records
    orig_members = {m.unique_id for m in half_dstore1.members}
    got = process1.apply_to(dstore1, id_from_source=get_data_source)
    assert {m.unique_id for m in got.members} == orig_members

    half_dstore2 = open_data_store(
        half_dstore2.source,
        suffix=half_dstore2.suffix,
        mode=APPEND,
    )

    reader2 = io_app.load_aligned(format="fasta", moltype="dna")
    min_length2 = sample_app.min_length(10)
    writer2 = io_app.write_seqs(half_dstore2)
    process2 = reader2 + min_length2 + writer2
    # check not fail on append new records
    _ = process2.apply_to(dstore1, id_from_source=get_data_source)


def test_skip_not_completed():
    @define_app(skip_not_completed=False)
    def takes_not_completed(val: c3types.SerialisableType) -> dict:
        return val.to_rich_dict()

    app = takes_not_completed()
    nc = NotCompleted("ERROR", "test", "for tracing", source="blah")
    got = app(nc)
    assert isinstance(got, dict)
    assert got == nc.to_rich_dict()


def test_copies_doc_from_func():
    @define_app
    def delme(val: c3types.SerialisableType) -> dict:
        """my docstring"""
        return val.to_rich_dict()

    assert delme.__doc__ == "my docstring"

    @define_app
    def delme2(val: c3types.SerialisableType) -> dict:
        """my docstring
        Notes
        -----
        body
        """
        return val.to_rich_dict()

    assert delme2.__doc__ == "my docstring"
    assert delme2.__init__.__doc__.split() == ["Notes", "-----", "body"]


def test_bad_wrap():
    def foo(a: "str") -> int:
        return int(a)

    with pytest.raises(NotImplementedError):
        define_app(foo)

    def bar(a: str) -> "int":
        return int(a)

    with pytest.raises(NotImplementedError):
        define_app(bar)

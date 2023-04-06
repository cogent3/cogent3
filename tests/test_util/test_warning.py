#!/usr/bin/env python

import pickle

import pytest

from cogent3.util.warning import deprecated_args


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_function_deprecated_args():
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(a: int, b: int) -> int:
        return a + b

    expected = changed(a=5, b=3)
    got = changed(x=5, y=3)
    assert got == expected


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_function_deprecated_args_docstring():
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(a: int, b: int) -> int:
        """This is a test function"""
        return a + b

    assert changed.__doc__ == "This is a test function"


@pytest.mark.parametrize("kwargs", (dict(x=5, y=3), dict(a=5, y=3), dict(x=5, b=3)))
def test_function_deprecated_args_warn(kwargs):
    # Example target function to be decorated
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with pytest.deprecated_call():
        changed(**kwargs)


def test_function_correct_args_do_not_warn():
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with pytest.warns(None) as record:  # verify no warnings for correct parameters
        changed(a=5, b=3)
    assert not record


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_function_deprecated_args_pickled():
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(a: int, b: int) -> int:
        return a + b

    myfunc = changed(x=1, y=2)
    pickled_func = pickle.dumps(myfunc)
    assert type(pickled_func) == bytes
    unpickled_func = pickle.loads(pickled_func)
    assert unpickled_func == changed(a=1, b=2)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args():
    class foo:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def changed(self, a: int, b: int) -> int:
            return a + b

    expected = foo().changed(a=5, b=3)
    got = foo().changed(x=5, y=3)
    assert got == expected


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args_docstring():
    class foo:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def changed(self, a: int, b: int) -> int:
            """This is a test function"""
            return a + b

    assert foo().changed.__doc__ == "This is a test function"


@pytest.mark.parametrize("kwargs", (dict(x=5, y=3), dict(a=5, y=3), dict(x=5, b=3)))
def test_method_deprecated_args_warn(kwargs):
    class foo:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def changed(self, a: int, b: int) -> int:
            """This is a test function"""
            return a + b

    with pytest.deprecated_call():
        foo().changed(**kwargs)


def test_method_correct_args_do_not_warn():
    class foo:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def changed(self, a: int, b: int) -> int:
            """This is a test function"""
            return a + b

    with pytest.warns(None) as record:  # verify no warnings for correct parameters
        foo().changed(a=5, b=3)
    assert not record


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args_pickled():
    class foo:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def changed(self, a: int, b: int) -> int:
            return a + b

    myfunc = foo().changed(x=1, y=2)
    pickled_func = pickle.dumps(myfunc)
    assert type(pickled_func) == bytes
    unpickled_func = pickle.loads(pickled_func)
    assert unpickled_func == foo().changed(a=1, b=2)

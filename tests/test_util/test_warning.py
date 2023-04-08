#!/usr/bin/env python

import pickle

import pytest

from cogent3.util.warning import deprecated_args


__author__ = "Richard Morris"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Richard Morris"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


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
    assert isinstance(pickled_func, bytes)
    unpickled_func = pickle.loads(pickled_func)
    assert unpickled_func == changed(a=1, b=2)


class foo:
    def __init__(self):
        self.a = 0
        self.b = 0

    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def changed(self, a: int, b: int):
        """This is a test function"""
        self.a = a
        self.b = b


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args():
    foo_instance = foo()
    foo_instance.changed(a=5, b=3)
    assert foo_instance.a == 5
    assert foo_instance.b == 3


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args_docstring():
    assert foo.changed.__doc__ == "This is a test function"
    foo_instance = foo()
    assert foo_instance.changed.__doc__ == "This is a test function"


@pytest.mark.parametrize("kwargs", (dict(x=5, y=3), dict(a=5, y=3), dict(x=5, b=3)))
def test_method_deprecated_args_warn(kwargs):

    with pytest.deprecated_call():
        foo().changed(**kwargs)


def test_method_correct_args_do_not_warn():
    with pytest.warns(None) as record:  # verify no warnings for correct parameters
        foo().changed(a=5, b=3)
    assert not record


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_method_deprecated_args_pickled():

    foo_instance = foo()
    assert foo_instance.a == 0
    assert foo_instance.b == 0
    foo_instance.changed(x=1, y=2)
    assert foo_instance.a == 1
    assert foo_instance.b == 2

    # serialize instance with pickle
    pickled_foo = pickle.dumps(foo_instance)
    assert isinstance(pickled_foo, bytes)

    # deserialize instance with pickle
    unpickled_foo = pickle.loads(pickled_foo)
    assert isinstance(unpickled_foo, foo)
    assert unpickled_foo.a == 1
    assert unpickled_foo.b == 2
    unpickled_foo.changed(x=2, y=3)
    assert unpickled_foo.a == 2
    assert unpickled_foo.b == 3

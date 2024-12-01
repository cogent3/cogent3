import pickle
import warnings

import pytest

from cogent3.util.warning import deprecated_args, deprecated_callable


def test_function_deprecated_args():
    @deprecated_args(
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with pytest.deprecated_call():
        expected = changed(a=5, b=3)
        got = changed(x=5, y=3)
        assert got == expected


def test_function_deprecated_args_docstring():
    @deprecated_args(
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(a: int, b: int) -> int:
        """This is a test function"""
        return a + b

    assert changed.__doc__ == "This is a test function"


@pytest.mark.parametrize("kwargs", (dict(x=5, y=3), dict(a=5, y=3), dict(x=5, b=3)))
def test_function_deprecated_args_warn(kwargs):
    # Example target function to be decorated
    @deprecated_args(
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with pytest.deprecated_call():
        changed(**kwargs)


def test_function_correct_args_do_not_warn():
    @deprecated_args(
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with warnings.catch_warnings():
        # verify no warnings for correct parameters
        warnings.simplefilter("error")
        changed(a=5, b=3)


def test_function_deprecated_args_pickled():
    @deprecated_args(
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(a: int, b: int) -> int:
        return a + b

    with pytest.deprecated_call():
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
        version="a future release",
        reason="x and y are not descriptive",
        old_new=[("x", "a"), ("y", "b")],
    )
    def changed(self, a: int, b: int):
        """This is a test function"""
        self.a = a
        self.b = b


def test_method_deprecated_args():
    foo_instance = foo()
    foo_instance.changed(a=5, b=3)
    assert foo_instance.a == 5
    assert foo_instance.b == 3


def test_method_deprecated_args_docstring():
    assert foo.changed.__doc__ == "This is a test function"
    foo_instance = foo()
    assert foo_instance.changed.__doc__ == "This is a test function"


@pytest.mark.parametrize("kwargs", (dict(x=5, y=3), dict(a=5, y=3), dict(x=5, b=3)))
def test_method_deprecated_args_warn(kwargs):
    with pytest.deprecated_call():
        foo().changed(**kwargs)  # pylint: disable=no-value-for-parameter


def test_method_correct_args_do_not_warn():
    with warnings.catch_warnings():
        # verify no warnings for correct parameters
        warnings.simplefilter("error")
        foo().changed(a=5, b=3)


def test_method_deprecated_args_pickled():
    foo_instance = foo()
    assert foo_instance.a == 0
    assert foo_instance.b == 0
    with pytest.deprecated_call():
        foo_instance.changed(x=1, y=2)  # pylint: disable=no-value-for-parameter
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
    with pytest.deprecated_call():
        unpickled_foo.changed(x=2, y=3)
    assert unpickled_foo.a == 2
    assert unpickled_foo.b == 3


class foo2:
    # noting that new_meth does not exist
    @deprecated_callable("2023.9", reason="test meth", new="new_meth")
    def old_meth(self, v):
        return v**2

    @deprecated_callable("2023.9", reason="redundant", is_discontinued=True)
    def squared(self, v):
        return v * v


@deprecated_callable("2023.9", reason="test func", new="new_func")
def old_func(v):
    return v**3


@deprecated_callable("2023.9", reason="redundant", is_discontinued=True)
def cubed(v):
    return v * v


@pytest.mark.parametrize("func", (foo2().old_meth, foo2().squared, old_func, cubed))
def test_deprecated_callable_warn(func):
    with pytest.deprecated_call():
        func(2)


@pytest.mark.parametrize("func", (cubed, old_func))
def test_method_deprecated_function_pickling(recwarn, func):
    # using pytest recwarn fixture here since the filterwarning decorator
    # does not seem to play nice with parametrize
    # serialize func with pickle
    pickled_func = pickle.dumps(func)
    assert isinstance(pickled_func, bytes)

    # deserialize with pickle
    unpickled_func = pickle.loads(pickled_func)

    assert unpickled_func(20) == func(20)


def test_method_deprecated_method_pickling(recwarn):
    # serialize instance with pickle
    instance = foo2()
    pickled = pickle.dumps(instance)
    assert isinstance(pickled, bytes)

    # deserialize with pickle
    unpickled = pickle.loads(pickled)

    # the method still works after pickling
    assert unpickled.old_meth(20) == instance.old_meth(20)

    # the method still works after pickling
    assert unpickled.squared(20) == instance.squared(20)


@pytest.mark.parametrize(
    "func,_type",
    ((foo2().old_meth, "method"), (cubed, "function")),
)
def test_deprecated_callable_resolves_type(recwarn, func, _type):
    func(2)
    assert _type in recwarn.list[0].message.args[0]


def test_function_deprecated_args_deprecated_callable_chained_decorators(recwarn):
    @deprecated_args("2023.6", "x is not descriptive", old_new=[("x", "a")])
    @deprecated_args("2023.6", "b is no longer required", discontinued=["b"])
    @deprecated_callable(
        "2023.6",
        "Improved change function",
        new="changed2",
        is_discontinued=True,
    )
    def changed(a: int, b: int) -> int:
        return a + b

    got = changed(x=5, b=3)  # pylint: disable=no-value-for-parameter
    assert got == 8
    warnings = [warning.message.args[0] for warning in recwarn.list]
    assert any("argument x which will be removed" in warning for warning in warnings)
    assert any("argument b is discontinued" in warning for warning in warnings)
    assert any("function changed is discontinued" in warning for warning in warnings)


def test_class_deprecated(recwarn):
    class fooclass:
        @deprecated_callable("2023.6", "Improved change function", is_discontinued=True)
        def __init__(self): ...

    fooclass()
    warnings = [warning.message.args[0] for warning in recwarn.list]
    assert any("fooclass is discontinued" in warning for warning in warnings)

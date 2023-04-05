#!/usr/bin/env python

"""Unit tests for union_dict.
"""
from unittest import TestCase, main

import pytest

from cogent3.util.warning import deprecated_args


def test_deprecated_args():
    # Example target function to be decorated
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def test_function(a: int, b: int) -> int:
        return a + b

    expected = test_function(a=5, b=3)
    got = test_function(x=5, y=3)
    assert got == expected


def test_deprecated_args_emits_warning():
    # Example target function to be decorated
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def test_function(a: int, b: int) -> int:
        return a + b

    with pytest.deprecated_call():
        test_function(x=5, y=3)

    with pytest.deprecated_call():
        test_function(a=5, y=3)
    with pytest.deprecated_call():
        test_function(x=5, b=3)

    with pytest.warns(None) as record:  # verify no warnings for correct parameters
        test_function(a=5, b=3)
    assert not record


def test_function_deprecated_args_maintains_docstring():
    @deprecated_args(
        [("x", "a"), ("y", "b")],
        version="a future release",
        reason="x and y are not descriptive",
    )
    def test_function(a: int, b: int) -> int:
        """This is a test function"""
        return a + b

    assert test_function.__doc__ == "This is a test function"


def test_method_deprecated_args():
    class test_class:
        @deprecated_args(
            [("x", "a"), ("y", "b")],
            version="a future release",
            reason="x and y are not descriptive",
        )
        def test_method(self, a: int, b: int) -> int:
            return a + b

    test_object = test_class()
    expected = test_object.test_method(a=5, b=3)
    got = test_object.test_method(x=5, y=3)
    assert got == expected

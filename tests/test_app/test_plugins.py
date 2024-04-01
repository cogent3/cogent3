import json
import sys
import time

from importlib.metadata import EntryPoint
from inspect import getsourcelines, isclass
from os import path
from shutil import rmtree
from subprocess import check_call
from tempfile import mkdtemp
from unittest.mock import Mock, patch

import pkg_resources
import pytest

from stevedore import extension
from stevedore.extension import ExtensionManager

import cogent3

from cogent3.app import (
    _make_apphelp_docstring,
    app_help,
    apps,
    available_apps,
    get_app,
)
from cogent3.app import typing as c3types
from cogent3.app.composable import define_app


@pytest.fixture
def extension_manager_factory():
    """Fixture to create mocked ExtensionManager instances with given extensions."""

    def _factory(extensions):
        with patch("stevedore.ExtensionManager") as mock_manager_constructor:
            mock_manager = ExtensionManager.make_test_instance(
                extensions=extensions, namespace="TESTING"
            )
            mock_manager_constructor.return_value = mock_manager
            return mock_manager

    return _factory

@pytest.fixture
def mock_extension_manager(extension_manager_factory, monkeypatch):
    """Fixture to mock the ExtensionManager with the given extensions."""
    def _mock_extension_manager(extensions):
        # Create a mocked ExtensionManager with the specified mock extensions
        mocked_manager = extension_manager_factory(extensions)
        # Patch the __apps variable in cogent3.apps module to use the mocked_manager
        monkeypatch.setattr(cogent3.app, "__apps", mocked_manager)
        return mocked_manager
    return _mock_extension_manager


def create_extension(plugin: object, name: str = None,  module_name:str = "module1") -> extension.Extension:
    if name is None:
        name = plugin.__name__
    return extension.Extension(
        name=name,
        entry_point=EntryPoint(name=name, value=f"{module_name}:{name}", group="TESTING"),
        plugin=plugin,
        obj=None,
    )

def test_Install_app(mock_extension_manager):
    @define_app
    class uppercase:
        """Convert a string to uppercase."""

        def main(self, data: str) -> str:
            return data.upper()
    mock_extension_manager([create_extension(uppercase)])

    appercase = get_app("uppercase")
    assert appercase("hello") == "HELLO"
    assert appercase.__doc__ in _make_apphelp_docstring(appercase.__class__)


@pytest.mark.parametrize("app_doc", [None, "text"])
@pytest.mark.parametrize("init_doc", [None, "text"])
def test_app_docs(mock_extension_manager, app_doc, init_doc):
    @define_app
    class documented_app:
        """This is a test app that has a __init__"""

        def __init__(self):
            self.constant = 2

        def main(self, val: int) -> int:
            return val + self.constant
    mock_extension_manager([create_extension(documented_app)])

    assert cogent3.app.apps().names() == ["documented_app"]
    app = get_app("documented_app")
    app.__class__.__doc__ = app_doc
    app.__class__.__init__.__doc__ = init_doc
    app_help("documented_app")
    got = _make_apphelp_docstring(app.__class__)
    assert "Options" in got


def test_namespace_collision(mock_extension_manager):
    @define_app
    class app1:
        def main(self, data: str) -> str:
            return data.upper()
    @define_app
    class app2:
        def main(self, data: str) -> str:
            return data.lower()
    # create two apps with the same name and different modules
    mock_extension_manager([create_extension(app1, module_name="module1"), create_extension(app2, name="app1",module_name="module2")])

    assert cogent3.app.apps().names() == ["app1", "app1"]

    with pytest.raises(NameError):
        app_by_name = get_app("app1")

    app_by_module_name_1 = get_app("module1.app1")
    app_by_module_name_2 = get_app("module2.app1")
    assert app_by_module_name_1("Hello") == "HELLO"
    assert app_by_module_name_2("Hello") == "hello"

    composition = app_by_module_name_1 + app_by_module_name_2
    assert composition("Hello") == "hello"

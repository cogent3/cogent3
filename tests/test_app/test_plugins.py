from importlib.metadata import EntryPoint
import random
import string
from unittest.mock import patch

import pytest

from stevedore import extension
from stevedore.extension import ExtensionManager

import cogent3

from cogent3.app import (
    _get_extension_attr,
    _make_apphelp_docstring,
    app_help,
    available_apps,
    get_app,
    get_app_manager,
)
from cogent3.app.composable import define_app
from cogent3.util.table import Table


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


def create_extension(
    plugin: object, name: str = None, module_name: str = "module1"
) -> extension.Extension:
    if name is None:
        name = plugin.__name__
    return extension.Extension(
        name=name,
        entry_point=EntryPoint(
            name=name, value=f"{module_name}:{name}", group="TESTING"
        ),
        plugin=plugin,
        obj=None,
    )


def test_install_app_class(mock_extension_manager):
    @define_app
    class uppercase:
        """Test app that converts a string to uppercase"""

        def main(self, data: str) -> str:
            return data.upper()

    mock_extension_manager([create_extension(uppercase)])

    appercase = get_app("uppercase")
    assert appercase("hello") == "HELLO"
    assert appercase.__doc__ in _make_apphelp_docstring(appercase.__class__)


def test_install_app_function(mock_extension_manager):
    @define_app
    def uppercase(data: str) -> str:
        """Test function that converts a string to uppercase"""
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
        """This is a test app that has a __init__, and a docstring"""

        def __init__(self):
            self.constant = 2

        def main(self, val: int) -> int:
            return val + self.constant

    mock_extension_manager([create_extension(documented_app)])

    assert cogent3.app.get_app_manager().names() == ["documented_app"]
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
    mock_extension_manager(
        [
            create_extension(app1, module_name="module1"),
            create_extension(app2, name="app1", module_name="module2"),
        ]
    )

    assert cogent3.app.get_app_manager().names() == ["app1", "app1"]

    with pytest.raises(NameError):
        _ = get_app(
            "app1"
        )  # request app by name only, when there are multiple apps with the same name

    app_by_module_name_1 = get_app("module1.app1")  # request app by name and module
    app_by_module_name_2 = get_app("module2.app1")
    assert app_by_module_name_1("Hello") == "HELLO"
    assert app_by_module_name_2("Hello") == "hello"

    composition = app_by_module_name_1 + app_by_module_name_2
    assert composition("Hello") == "hello"


def test_available_apps_local(mock_extension_manager):
    """available_apps robust to local scope apps"""

    @define_app
    def dummy(val: int) -> int:
        return val

    mock_extension_manager([create_extension(dummy)])
    apps = available_apps()
    assert isinstance(apps, Table)
    apps = apps.filtered(lambda x: dummy.__name__ == x, columns="name")
    assert apps.shape[0] == 1

def test_stevedore_finds_non_apps(mock_extension_manager):
    """available_apps should return plugins that fail is_app() but emit a warning"""

    def not_an_app(val: int) -> int:
        return val

    not_an_app_extension = create_extension(not_an_app)
    with pytest.warns(UserWarning, match=r".* is not a valid cogent3 app, skipping"):
        mock_extension_manager([create_extension(not_an_app)])
        apps = available_apps()
        assert isinstance(apps, Table)
        apps = apps.filtered(lambda x: not_an_app.__name__ == x, columns="name")
        assert apps.shape[0] == 1
        
def test_apps_can_be_deprecated(mock_extension_manager):
    """if an app is in __deprecated[] then it should not appear in available_apps"""
    
    @define_app
    def deprecated_app(val: int) -> int:
        return val

    mock_extension_manager([create_extension(deprecated_app)])
    assert deprecated_app.__name__ in cogent3.app.get_app_manager().names()

    with patch("cogent3.app.__deprecated", new=[deprecated_app.__name__]):
        apps = available_apps()
        assert deprecated_app.__name__ not in apps.columns["name"]
        assert deprecated_app.__name__ in cogent3.app.get_app_manager().names()

    # check that the app is still available if we ask for it by name
    app = get_app(deprecated_app.__name__)
    assert app(5) == 5  
    
def test_unknown_app_name(mock_extension_manager):
    """get_app should raise a ValueError if the app name is not found"""

    mock_extension_manager([])
    unknown_app_name = "".join(random.choices(string.ascii_lowercase, k=10))

    with pytest.raises(ValueError, match=f"App '{unknown_app_name}' not found. Please check for typos."):
        _ = get_app(unknown_app_name)

def test_unknown_module_name(mock_extension_manager):
    """get_app should raise a ValueError if the app name is not found"""

    @define_app
    def dummy(val: int) -> int:
        return val

    mock_extension_manager([create_extension(dummy,module_name="module1")])

    assert dummy.__name__ in cogent3.app.get_app_manager().names()
    assert get_app(dummy.__name__)(5) == 5
    with pytest.raises(ValueError, match=".* not found. Please check for typos."):
        _ = get_app(''.join(['module_2',dummy.__name__]))
